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
    : zmqContext(1), fSocket(zmqContext, ZMQ_ROUTER), fConfig(config){

  jobPool = config->jobPool;
}

json GeantEventDispatcher::HandleNewWorkerMsg(json &req){
  workerCounter++;
  GeantHPCWorker worker;
  worker.id = workerCounter;
  worker.zmqID = req["zmq_id"].get<std::string>();
  worker.lastContact.Update();
  workers[worker.id] = worker;
  workerId[worker.zmqID] = worker.id;
  json reply;
  reply["type"] = "new_wrk_rep";
  reply["wrk_id"] = worker.id;
  return reply;
}

json GeantEventDispatcher::HandleNewJobReq(json &req){
  int n = req["events"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  if(workers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = workers[wrk_id];
  worker.lastContact.Update();
  GeantHPCJob job = jobPool->GetJob(n,worker);
  if (job.type == JobType::HepMC) {
    if (job.hepMCJobs.size() == 0) {
      bool foundDupeJob = false;
      for(auto pendIt = pendingJobs.rbegin(); pendIt != pendingJobs.rend(); pendIt++){
        auto& pendJob = pendIt->second;
        if(IsWorkerDoingJob(worker,pendJob)) continue;
        if(n < pendJob.hepMCJobs.size()) continue;
        job = jobPool->GetDublicateJob(pendJob);
        foundDupeJob = true;
        std::cout << "Duplicating job: " << pendJob.uid << '\n';
        std::cout << "Dupe id: ";
        for(auto dum : *job.dublicateUIDs){
          std::cout << dum << ' ';
        }
        std::cout << '\n';
        break;
      }
      if(!foundDupeJob) return "{\"type\": \"wait\"}"_json;
    }
    job.workerId = worker.id;
    job.dispatchTime.Update();
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
      bool foundDupeJob = false;
      for(auto pendIt = pendingJobs.rbegin(); pendIt != pendingJobs.rend(); pendIt++){
        auto& pendJob = pendIt->second;
        if(IsWorkerDoingJob(worker,pendJob)) continue;
        if(n < pendJob.events) continue;
        job = jobPool->GetDublicateJob(pendJob);
        foundDupeJob = true;
        std::cout << "Duplicating job: " << pendJob.uid << '\n';
        std::cout << "Dupe id: ";
        for(auto dum : *job.dublicateUIDs){
          std::cout << dum << ' ';
        }
        std::cout << '\n';
        break;
      }
      if(!foundDupeJob) return "{\"type\": \"wait\"}"_json;
    }
    job.workerId = worker.id;
    job.dispatchTime.Update();
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

bool GeantEventDispatcher::IsWorkerDoingJob(GeantHPCWorker& worker, GeantHPCJob& job){
  if(job.workerId == worker.id) return true;
  for(int jobUid : *job.dublicateUIDs){
    if(pendingJobs.count(jobUid)>0) {
      if (pendingJobs[jobUid].workerId == worker.id)
        return true;
    }
  }
  return false;
}

json GeantEventDispatcher::HandleJobConfirm(json &req){
  int job_id = req["job_id"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  if(workers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = workers[wrk_id];

  worker.lastContact.Update();

  if(pendingJobs.count(job_id) == 0)
    return "{\"type\": \"error\"}"_json;
  auto& job = pendingJobs[job_id];

  worker.expectedTimeForEvent = (worker.transportedEvents)*(worker.expectedTimeForEvent) +
                                 job.dispatchTime.Since();
  worker.transportedEvents += (job.events + job.hepMCJobs.size());
  worker.expectedTimeForEvent /= worker.transportedEvents;


  job.dublicateUIDs->erase(job_id);
  for(auto dupeId : *job.dublicateUIDs){
    SendJobCancelMsg(pendingJobs[dupeId], false);
  }
  pendingJobs.erase(job_id);
  json rep;
  rep["type"] = "job_done_rep";
  return rep;
}

json GeantEventDispatcher::HandleHeartbeat(json &req) {
  int wrk_id = req["wrk_id"].get<int>();
  if(workers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = workers[wrk_id];
  worker.lastContact.Update();
  return "{\"type\": \"hb_rep\"}"_json;
}

void GeantEventDispatcher::Initialize()
{
  BindSocket();

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
  while(!jobPool->IsEmpty() || !pendingJobs.empty()){
    PollForMsg();
    ResendMsg();
    CleanDeadWorkers();
  }

  FinishWorkers();
}

void GeantEventDispatcher::CleanDeadWorkers(){
  for(auto w_it = workers.begin(); w_it!=workers.end();){
    auto w = &w_it->second;
    auto since = w->lastContact.Since();
    if(since > std::chrono::seconds(45)){
      for(auto job_it = pendingJobs.begin(); job_it != pendingJobs.end();){
        if (job_it->second.workerId == w->id) {
          job_it->second.dublicateUIDs->erase(job_it->second.uid);
          if (job_it->second.dublicateUIDs->size() == 0){
            jobPool->ReturnToPool(job_it->second);
          }
          job_it = pendingJobs.erase(job_it);
        } else{
          ++job_it;
        }
      }
      std::cout << "mast: worker is dead id: " << w->id << std::endl;
      if(workerId[w->zmqID] == w->id)
        workerId.erase(w->zmqID);
      w_it = workers.erase(w_it);
    }else{
      ++w_it;
    }
  }
}

void GeantEventDispatcher::BindSocket() {
  zmqPollItem =  { fSocket, 0, ZMQ_POLLIN, 0 };
  fSocket.bind("tcp://*:"+std::to_string(fConfig->fMasterPort));
  std::cout << "HPC: Event dispatcher bounded to tcp://*:"+std::to_string(fConfig->fMasterPort) << std::endl;
}

void GeantEventDispatcher::SendMessage(const std::string &msg, const std::string &type, size_t uid,
                                       const std::string &address) {

  zmq::message_t m_addr(address.size());
  memcpy(m_addr.data(),address.c_str(),address.size());
  zmq::message_t m_type(type.size());
  memcpy(m_type.data(),type.c_str(),type.size());
  zmq::message_t m_uid(sizeof(size_t));
  memcpy(m_uid.data(),&uid,sizeof(size_t));
  zmq::message_t message(msg.size());
  memcpy(message.data(),msg.c_str(),msg.size());

  fSocket.send(m_addr,ZMQ_SNDMORE);
  fSocket.send(m_type,ZMQ_SNDMORE);
  fSocket.send(m_uid,ZMQ_SNDMORE);
  fSocket.send(message);

  std::cout << "send: " << address;
  std::cout << " " <<  type << '\n';
  std::cout << msg << '\n';
}


void GeantEventDispatcher::SendRep(const std::string &msg, size_t uid, const std::string &address) {
  SendMessage(msg,"REP",uid,address);
}

string GeantEventDispatcher::RecvReq(const std::string &msg) {
  json req = json::parse(msg);
  json rep;
  rep["type"] = "error";

  auto req_type = req["type"].get<std::string>();
  if (req_type == "new_wrk_req") {
    rep = HandleNewWorkerMsg(req);
  } else if (req_type == "job_req") {
    rep = HandleNewJobReq(req);
  } else if (req_type == "job_done_req") {
    rep = HandleJobConfirm(req);
  } else if (req_type == "hb") {
    rep = HandleHeartbeat(req);
  }
  return rep.dump();
}

void GeantEventDispatcher::SendReq(const std::string &msg, GeantHPCWorker &worker) {
  if(worker.pendingMsg.size() >= 3)
    return;
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  size_t msgUid = std::hash<std::string>{}(msg+std::to_string(time));
  pendingRequests[msgUid] = MasterPendingMessage(msg,worker.id);
  worker.pendingMsg.push_back(msgUid);
  SendMessage(msg,"REQ",msgUid,worker.zmqID);
}

void GeantEventDispatcher::RecvRep(const string &msg) {
  return;
}

void GeantEventDispatcher::PollForMsg() {
  zmq::poll(&zmqPollItem,1,std::chrono::milliseconds(10));
  if (zmqPollItem.revents & ZMQ_POLLIN) {
    zmq::message_t addr;
    zmq::message_t type;
    zmq::message_t muid;
    zmq::message_t message;
    fSocket.recv(&addr);
    fSocket.recv(&type);
    fSocket.recv(&muid);
    fSocket.recv(&message);
    std::string messageAddr((char*)addr.data(),(char*)addr.data() + addr.size());
    std::string messageType((char*)type.data(),(char*)type.data() + type.size());
    size_t messageUid;
    memcpy(&messageUid,muid.data(),muid.size());
    std::string messageContent((char*)message.data(),(char*)message.data() + message.size());
    std::cout << "recv: " << messageAddr;
    std::cout << " " <<  messageType << '\n';
    std::cout << messageContent << '\n';
    if(messageType == "REQ"){
      std::string replyContent = RecvReq(messageContent);
      SendRep(replyContent,messageUid,messageAddr);
    }else { //messageType == REP
      if(pendingRequests.count(messageUid) > 0){
        pendingRequests.erase(messageUid);
        if(workerId.count(messageAddr) > 0 && workers.count(workerId[messageAddr]) > 0){
          auto& workerPending = workers[workerId[messageAddr]].pendingMsg;
          workerPending.erase(std::find(workerPending.begin(),workerPending.end(),messageUid));
        }
        RecvRep(messageContent);
      }
    }
  }
}

bool GeantEventDispatcher::ResendMsg() {
  bool ok = true;
  for(auto it = pendingRequests.begin(); it != pendingRequests.end();){
    auto& pend = it->second;
    bool remove = false;
    if(workers.count(pend.workerId) > 0) {
      if (pend.lastRetry.Since() > std::chrono::seconds(5)) {
        if (pend.retries >= 1) {
          remove = true;
        } else {
          std::cout << "Resending msg: " << pend.msg << std::endl;
          SendMessage(pend.msg, "REQ", it->first, workers[pend.workerId].zmqID);
          pend.retries++;
          pend.lastRetry.Update();
        }
      }
    }else{
      remove = true; //Worker died - removing request
    }
    if(remove) {
      std::cout << "Removing msg: " << pend.msg << std::endl;
      it = pendingRequests.erase(it);
      workers[pend.workerId].discardedMsg++;
      ok = false;
    }
    else {
      it++;
    }
  }
  return ok;
}

void GeantEventDispatcher::FinishWorkers() {
  for(auto& it : workers){
    auto& w = it.second;
    SendFinishMsg(w);
  }
  while(pendingRequests.size() > 0){
    PollForMsg();
    ResendMsg();
  }
}

void GeantEventDispatcher::SendFinishMsg(GeantHPCWorker &worker) {
  SendReq("{\"type\":\"finish_req\"}",worker);
}

void GeantEventDispatcher::SendJobCancelMsg(GeantHPCJob &job, bool retToPool) {
  json req;
  req["type"] = "job_cancel_req";
  req["wrk_id"] = job.workerId;
  req["job_id"] = job.uid;
  if(pendingJobs.count(job.uid)) {
    if(retToPool) {
      jobPool->ReturnToPool(pendingJobs[job.uid]);
    }
    pendingJobs.erase(job.uid);
  }
  SendReq(req.dump(),workers[job.workerId]);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
