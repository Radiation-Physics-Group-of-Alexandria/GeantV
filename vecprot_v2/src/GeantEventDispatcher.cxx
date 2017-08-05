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
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REP), fConfig(config) {

  jobPool = config->jobPool;
}

json GeantEventDispatcher::NewWorker(json& req){
  workerCounter++;
  GeantHPCWorker worker;
  worker.id = workerCounter;
  worker.reqSocket = req["adr"].get<std::string>();
  worker.lastContact = std::chrono::system_clock::now();
  workers.push_back(worker);
  json reply;
  reply["type"] = "new_wrk_rep";
  reply["wrk_id"] = worker.id;
  return reply;
}

json GeantEventDispatcher::JobReq(json& req){
  int n = req["events"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  auto worker = std::find_if(workers.begin(),workers.end(),[&](const GeantHPCWorker& w){return w.id == wrk_id;});
  if(worker == workers.end()){
    return "{\"type\": \"error\"}"_json;
  }
  worker->lastContact = std::chrono::system_clock::now();
  GeantHPCJob job = jobPool->GetJob(n,*worker);
  if (job.type == JobType::HepMC) {
    if (job.hepMCJobs.size() == 0) {
      return "{\"type\": \"wait\"}"_json;
    }
    job.workerId = worker->id;
    pendingJobs.push_back(job);
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
    job.workerId = worker->id;
    pendingJobs.push_back(job);
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
  int worker_id = req["wrk_id"].get<int>();
  auto worker = std::find_if(workers.begin(),workers.end(),[&](const GeantHPCWorker& w){return w.id == worker_id;});
  if(worker == workers.end())
    return "{\"type\": \"error\"}"_json;

  worker->lastContact = std::chrono::system_clock::now();

  auto job = std::find_if(pendingJobs.begin(),pendingJobs.end(),[&](const GeantHPCJob& job){ return job.uid == job_id; });
  if(job == pendingJobs.end())
    return "{\"type\": \"error\"}"_json;
  pendingJobs.erase(job);

  json rep;
  rep["type"] = "job_done_rep";
  return rep;
}

json GeantEventDispatcher::HeartBeat(json &req) {
  int wrk_id = req["wrk_id"].get<int>();
  auto worker = std::find_if(workers.begin(),workers.end(),[&](const GeantHPCWorker& w){return w.id == wrk_id;});
  if(worker == workers.end())
    return "{\"type\": \"error\"}"_json;
  worker->lastContact = std::chrono::system_clock::now();
  return "{\"type\": \"hb_rep\"}"_json;
}

void GeantEventDispatcher::Initialize()
{
  zmqSocket.bind("tcp://*:5678");
  std::cout << "HPC: Event dispatcher bounded to \"tcp://*:5678\"" << std::endl;


}

void GeantEventDispatcher::RunReqReplyLoop()
{
  zmq_pollitem_t items [] = { { zmqSocket, 0, ZMQ_POLLIN, 0 } };
  while(!jobPool->IsEmpty() || !pendingJobs.empty()){

    zmq::poll(items,1,std::chrono::milliseconds(30));
    if (items[0].revents & ZMQ_POLLIN) {
      zmq::message_t request(1024);
      zmqSocket.recv(&request);
      char requestMsg[1024];
      memcpy(requestMsg, request.data(), 1024);
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

      zmq::message_t reply(1024);
      memcpy(reply.data(), responce.c_str(), responce.size() + 1);
      zmqSocket.send(reply);
    }
    CleanDeadWorkers();
  }


//TODO: Here we send "kill yourself" to all workers
}

void GeantEventDispatcher::CleanDeadWorkers(){
  auto now = std::chrono::system_clock::now();
  for(auto w = workers.begin(); w!=workers.end();){
    auto since = std::chrono::duration_cast<std::chrono::seconds>(now - w->lastContact);
    if(since > std::chrono::seconds(45)){
      for(auto it = pendingJobs.begin(); it != pendingJobs.end();){
        if (it->workerId == w->id){
          jobPool->ReturnToPool(*it);
          it = pendingJobs.erase(it);
        } else{
          ++it;
        }
      }
      std::cout << "mast: worker is dead id: " << w->id << std::endl;
      w = workers.erase(w);
    }else{
      ++w;
    }
  }
}




} // GEANT_IMPL_NAMESPACE
} // Geant
