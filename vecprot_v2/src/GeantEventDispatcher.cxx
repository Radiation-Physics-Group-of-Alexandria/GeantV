#include "GeantEventDispatcher.h"
#include "GeantEventServer.h"
#include "GeantEventReceiver.h"

#include "globals.h"
#include "Geant/Error.h"

#include "GeantTrack.h"
#include "GeantEvent.h"
#include "GeantRunManager.h"
#include "LocalityManager.h"
#include "PrimaryGenerator.h"
#include "GeantTaskData.h"
#include "GeantBasket.h"
#include "MCTruthMgr.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/SimpleNavigator.h"
#include "volumes/PlacedVolume.h"
#include "management/GeoManager.h"
#else
#ifdef USE_ROOT
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#endif
#endif

#ifdef USE_HPC
#include "zmq.hpp"
#include <json.hpp>
using nlohmann::json;

#include "GeantJobPool.h"
#endif

#include <thread>
#include <functional>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <fstream>
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GeantEventDispatcher::GeantEventDispatcher(GeantConfig *config)
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REP), fConfig(config), zmqSocketOut(zmqContext,ZMQ_REQ){

  jobPool = config->jobPool;
}

json GeantEventDispatcher::NewWorker(json& req){
  workerCounter++;
  GeantHPCWorker worker;
  worker.id = workerCounter;
  worker.reqSocket = req["adr"].get<std::string>();
  worker.lastContact = std::chrono::system_clock::now();
  workers[worker.id] = worker;
  json reply;
  reply["type"] = "new_wrk_rep";
  reply["wrk_id"] = worker.id;
  return reply;
}

json GeantEventDispatcher::JobReq(json& req){
  int n = req["events"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  if(workers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = workers[wrk_id];
  worker.lastContact = std::chrono::system_clock::now();
  GeantHPCJob job = jobPool->GetJob(n,worker);
  job.dispatchTime = std::chrono::system_clock::now();
  if (job.type == JobType::HepMC) {
    if (job.hepMCJobs.size() == 0) {
      return "{\"type\": \"wait\"}"_json;
    }
    job.workerId = worker.id;
    pendingJobs[job.uid]=job;
    json reply;
    reply["type"] = "job_rep";
    reply["wrk_id"] = wrk_id;
    reply["job_id"] = job.uid;
    reply["job_type"] = "hepmc";
    reply["files"] = job.hepMCJobs;
    return reply;
  } else if (job.type == JobType::Generator){
    if (job.events == 0) {
      return "{\"type\": \"wait\"}"_json;
    }
    job.workerId = worker.id;
    pendingJobs[job.uid]=job;
    json reply;
    reply["type"] = "job_rep";
    reply["wrk_id"] = wrk_id;
    reply["job_id"] = job.uid;
    reply["job_type"] = "generate";
    reply["event"] = job.events;
    return reply;
  }
}

json GeantEventDispatcher::JobConfirm(json& req){
  int job_id = req["job_id"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  if(workers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = workers[wrk_id];

  auto timeNow = std::chrono::system_clock::now();
  worker.lastContact = timeNow;

  if(pendingJobs.count(job_id) == 0)
    return "{\"type\": \"error\"}"_json;
  auto& job = pendingJobs[job_id];

  worker.expectedTimeForEvent = (worker.transportedEvents)*(worker.expectedTimeForEvent) +
                                 std::chrono::duration_cast<std::chrono::seconds>(timeNow - job.dispatchTime);
  worker.transportedEvents += (job.events + job.hepMCJobs.size());
  worker.expectedTimeForEvent /= worker.transportedEvents;


  pendingJobs.erase(job_id);
  json rep;
  rep["type"] = "job_done_rep";
  return rep;
}

json GeantEventDispatcher::HeartBeat(json &req) {
  int wrk_id = req["wrk_id"].get<int>();
  if(workers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = workers[wrk_id];
  worker.lastContact = std::chrono::system_clock::now();
  return "{\"type\": \"hb_rep\"}"_json;
}

void GeantEventDispatcher::Initialize()
{
  zmqSocket.bind("tcp://*:"+std::to_string(fConfig->fMasterPort));
  std::cout << "HPC: Event dispatcher bounded to tcp://*:"+std::to_string(fConfig->fMasterPort) << std::endl;

  std::ifstream inputFile(fConfig->fHostnameFile);
  std::string line;
  while(inputFile >> line) {
    Host host;
    host.hostname = line;
    fHosts.push_back(host);
  }

  for(auto& h : fHosts){
    std::cout << "Starting worker on remote host: " + h.hostname << std::endl;
    system((fConfig->fRemoteStartScript + " " + h.hostname + " " + fHosts[0].hostname + " &").c_str());
  }

}

void GeantEventDispatcher::RunReqReplyLoop()
{
  zmq_pollitem_t items [] = { { zmqSocket, 0, ZMQ_POLLIN, 0 } };
  while(!jobPool->IsEmpty() || !pendingJobs.empty()){

    zmq::poll(items,1,std::chrono::milliseconds(30));
    if (items[0].revents & ZMQ_POLLIN) {
      zmq::message_t request;
      zmqSocket.recv(&request);
      std::string requestMsg((char* )request.data(),(char* )request.data()+request.size());
      std::cout << "master req: " << requestMsg << std::endl;
      json rep;
      json req = json::parse(requestMsg);
      auto req_type = req["type"].get<std::string>();
      if (req_type == "new_wrk_req") {
        rep = NewWorker(req);
      } else if (req_type == "job_req") {
        rep = JobReq(req);
      } else if (req_type == "job_done_req") {
        rep = JobConfirm(req);
      } else if (req_type == "hb") {
        rep = HeartBeat(req);
      }

      std::string responce = rep.dump();
      std::cout << "master rep: " << responce << std::endl;

      zmq::message_t reply(responce.size());
      memcpy(reply.data(), responce.c_str(), responce.size());
      zmqSocket.send(reply);
    }
    CleanDeadWorkers();
    for(auto& j : pendingJobs){
      CancelJob(j.second);
    }
  }

  FinishWorkers();


}

void GeantEventDispatcher::CleanDeadWorkers(){
  auto now = std::chrono::system_clock::now();
  for(auto w_it = workers.begin(); w_it!=workers.end();){
    auto w = &w_it->second;
    auto since = std::chrono::duration_cast<std::chrono::seconds>(now - w->lastContact);
    if(since > std::chrono::seconds(45)){
      for(auto job_it = pendingJobs.begin(); job_it != pendingJobs.end();){
        if (job_it->second.workerId == w->id){
          jobPool->ReturnToPool(job_it->second);
          job_it = pendingJobs.erase(job_it);
        } else{
          ++job_it;
        }
      }
      std::cout << "mast: worker is dead id: " << w->id << std::endl;
      w_it = workers.erase(w_it);
    }else{
      ++w_it;
    }
  }
}

bool GeantEventDispatcher::SendMsgWorker(const std::string &adr, const std::string &req, std::string &rep) {
  zmqSocketOut.setsockopt(ZMQ_RCVTIMEO,5*1000);
  zmqSocketOut.connect(adr);
  zmq::message_t request(req.size());
  memcpy(request.data(),req.c_str(),req.size());
  std::cout << "mast to wrk: msg req " << req << std::endl;

  zmqSocketOut.send(request);
  zmq::message_t reply;
  bool ok = zmqSocketOut.recv(&reply);
  if (!ok) return false;

  std::string rep_msg((char* )reply.data(),(char* )reply.data() + reply.size());
  rep = rep_msg;
  std::cout << "mast to wrk: msg rep" << rep << std::endl;

  zmqSocketOut.disconnect(adr);
  return true;
}

void GeantEventDispatcher::FinishWorkers() {
  for(auto & worker_pair: workers){
    auto& worker = worker_pair.second;
    std::string msg = FinishMsg(worker).dump();
    std::string rep;
    SendMsgWorker(worker.reqSocket,msg,rep);
  }
}

json GeantEventDispatcher::FinishMsg(const GeantHPCWorker &worker) {
  json j;
  j["type"] = "finish";
  j["wrk_id"] = worker.id;
  return j;
}

json GeantEventDispatcher::CancelJobMsg(const GeantHPCWorker &worker, const GeantHPCJob &job) {
  json j;
  j["type"] = "cancel";
  j["wrk_id"] = worker.id;
  j["job_id"] = job.uid;
  return j;
}

void GeantEventDispatcher::CancelJob(const GeantHPCJob &job) {
  if (workers.count(job.workerId) != 0){
    std::string rep;
    bool ok = SendMsgWorker(workers[job.workerId].reqSocket,CancelJobMsg(workers[job.workerId],job).dump(),rep);
    if(!ok) ReplaceOutSocket();
  }
  if(pendingJobs.count(job.uid) != 0){
    jobPool->ReturnToPool(pendingJobs[job.uid]);
    pendingJobs.erase(job.uid);
  }
}

void GeantEventDispatcher::ReplaceOutSocket() {
  zmqSocketOut.close();
  zmqSocketOut = zmq::socket_t(zmqContext,ZMQ_REQ);
}


} // GEANT_IMPL_NAMESPACE
} // Geant
