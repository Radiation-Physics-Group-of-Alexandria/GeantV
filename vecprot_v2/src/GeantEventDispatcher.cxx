#include "GeantEventDispatcher.h"

#include <iostream>
#include <fstream>

#include "globals.h"
#include "Geant/Error.h"

#include "GeantJobPool.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GeantEventDispatcher::GeantEventDispatcher(GeantConfig *config)
    : fZMQContext(1), fSocket(fZMQContext, ZMQ_ROUTER), fConfig(config), fJobPool(config->jobPool),
      fWorkerCounter(0)
{

}

json GeantEventDispatcher::HandleNewWorkerMsg(json &req){
  fWorkerCounter++;
  GeantHPCWorker worker;
  worker.fID = fWorkerCounter;
  worker.fZMQID = req["zmq_id"].get<std::string>();
  worker.fLastContact.Update();
  fHPCWorkers[worker.fID] = worker;
  fHPCWorkerZMQIDtoUID[worker.fZMQID] = worker.fID;
  json reply;
  reply["type"] = "new_wrk_rep";
  reply["wrk_id"] = worker.fID;
  return reply;
}

json GeantEventDispatcher::HandleNewJobReq(json &req){
  int n = req["events"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  if(fHPCWorkers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = fHPCWorkers[wrk_id];
  worker.fLastContact.Update();
  GeantHPCJob job = fJobPool->GetJob(n,worker);
  if (job.fType == JobType::HepMC) {
    if (job.fHepMCJobs.size() == 0) {
      bool foundDupeJob = false;
      for(auto pendIt = fPendingJobs.rbegin(); pendIt != fPendingJobs.rend(); pendIt++){
        auto& pendJob = pendIt->second;
        if(IsWorkerDoingJob(worker,pendJob)) continue;
        if(n < pendJob.fHepMCJobs.size()) continue;
        job = fJobPool->GetDublicateJob(pendJob);
        foundDupeJob = true;
        std::cout << "Duplicating job: " << pendJob.fUID << '\n';
        std::cout << "Dupe id: ";
        for(auto dum : *job.fDublicateUIDs){
          std::cout << dum << ' ';
        }
        std::cout << '\n';
        break;
      }
      if(!foundDupeJob) return "{\"type\": \"wait\"}"_json;
    }
    job.fWorkerID = worker.fID;
    job.fDispatchTime.Update();
    fPendingJobs[job.fUID]=job;
    json reply;
    reply["type"] = "job_rep";
    reply["wrk_id"] = wrk_id;
    reply["job_id"] = job.fUID;
    reply["job_type"] = "hepmc";
    reply["files"] = job.fHepMCJobs;
    return reply;
  } else if (job.fType == JobType::Generator){
    if (job.fEvents == 0) {
      bool foundDupeJob = false;
      for(auto pendIt = fPendingJobs.rbegin(); pendIt != fPendingJobs.rend(); pendIt++){
        auto& pendJob = pendIt->second;
        if(IsWorkerDoingJob(worker,pendJob)) continue;
        if(n < pendJob.fEvents) continue;
        job = fJobPool->GetDublicateJob(pendJob);
        foundDupeJob = true;
        std::cout << "Duplicating job: " << pendJob.fUID << '\n';
        std::cout << "Dupe id: ";
        for(auto dum : *job.fDublicateUIDs){
          std::cout << dum << ' ';
        }
        std::cout << '\n';
        break;
      }
      if(!foundDupeJob) return "{\"type\": \"wait\"}"_json;
    }
    job.fWorkerID = worker.fID;
    job.fDispatchTime.Update();
    fPendingJobs[job.fUID]=job;
    json reply;
    reply["type"] = "job_rep";
    reply["wrk_id"] = wrk_id;
    reply["job_id"] = job.fUID;
    reply["job_type"] = "generate";
    reply["event"] = job.fEvents;
    return reply;
  }
  return "{\"type\": \"error\"}"_json;
}

bool GeantEventDispatcher::IsWorkerDoingJob(GeantHPCWorker& worker, GeantHPCJob& job){
  if(job.fWorkerID == worker.fID) return true;
  for(int jobUid : *job.fDublicateUIDs){
    if(fPendingJobs.count(jobUid)>0) {
      if (fPendingJobs[jobUid].fWorkerID == worker.fID)
        return true;
    }
  }
  return false;
}

json GeantEventDispatcher::HandleJobConfirm(json &req){
  int job_id = req["job_id"].get<int>();
  int wrk_id = req["wrk_id"].get<int>();
  if(fHPCWorkers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = fHPCWorkers[wrk_id];

  worker.fLastContact.Update();

  if(fPendingJobs.count(job_id) == 0)
    return "{\"type\": \"error\"}"_json;
  auto& job = fPendingJobs[job_id];

  worker.fExpectedTimeForEvent = (worker.fTransportedEvents)*(worker.fExpectedTimeForEvent) +
                                 job.fDispatchTime.Since();
  worker.fTransportedEvents += (job.fEvents + job.fHepMCJobs.size());
  worker.fExpectedTimeForEvent /= worker.fTransportedEvents;


  job.fDublicateUIDs->erase(job_id);
  for(auto dupeId : *job.fDublicateUIDs){
    SendJobCancelMsg(fPendingJobs[dupeId], false);
  }
  fPendingJobs.erase(job_id);
  json rep;
  rep["type"] = "job_done_rep";
  return rep;
}

json GeantEventDispatcher::HandleHeartbeat(json &req) {
  int wrk_id = req["wrk_id"].get<int>();
  if(fHPCWorkers.count(wrk_id) == 0){
    return "{\"type\": \"error\"}"_json;
  }
  auto& worker = fHPCWorkers[wrk_id];
  worker.fLastContact.Update();
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
  while(!fJobPool->IsEmpty() || !fPendingJobs.empty()){
    PollForMsg();
    ResendMsg();
    CleanDeadWorkers();
  }

  FinishWorkers();
}

void GeantEventDispatcher::CleanDeadWorkers(){
  for(auto w_it = fHPCWorkers.begin(); w_it!=fHPCWorkers.end();){
    auto w = &w_it->second;
    auto since = w->fLastContact.Since();
    if(since > std::chrono::seconds(45)){
      for(auto job_it = fPendingJobs.begin(); job_it != fPendingJobs.end();){
        if (job_it->second.fWorkerID == w->fID) {
          job_it->second.fDublicateUIDs->erase(job_it->second.fUID);
          if (job_it->second.fDublicateUIDs->size() == 0){
            fJobPool->ReturnToPool(job_it->second);
          }
          job_it = fPendingJobs.erase(job_it);
        } else{
          ++job_it;
        }
      }
      std::cout << "mast: worker is dead id: " << w->fID << std::endl;
      if(fHPCWorkerZMQIDtoUID[w->fZMQID] == w->fID)
        fHPCWorkerZMQIDtoUID.erase(w->fZMQID);
      w_it = fHPCWorkers.erase(w_it);
    }else{
      ++w_it;
    }
  }
}

void GeantEventDispatcher::BindSocket() {
  fZMQSocketPollItem =  { fSocket, 0, ZMQ_POLLIN, 0 };
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

std::string GeantEventDispatcher::RecvReq(const std::string &msg) {
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
  if(worker.fPendingRequests.size() >= 3)
    return;
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  size_t msgUid = std::hash<std::string>{}(msg+std::to_string(time));
  fPendingRequests[msgUid] = MasterPendingMessage(msg,worker.fID);
  worker.fPendingRequests.push_back(msgUid);
  SendMessage(msg,"REQ",msgUid,worker.fZMQID);
}

void GeantEventDispatcher::RecvRep(const std::string &msg) {
  return;
}

void GeantEventDispatcher::PollForMsg() {
  zmq::poll(&fZMQSocketPollItem,1,std::chrono::milliseconds(10));
  if (fZMQSocketPollItem.revents & ZMQ_POLLIN) {
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
      if(fPendingRequests.count(messageUid) > 0){
        fPendingRequests.erase(messageUid);
        if(fHPCWorkerZMQIDtoUID.count(messageAddr) > 0 && fHPCWorkers.count(fHPCWorkerZMQIDtoUID[messageAddr]) > 0){
          auto& workerPending = fHPCWorkers[fHPCWorkerZMQIDtoUID[messageAddr]].fPendingRequests;
          workerPending.erase(std::find(workerPending.begin(),workerPending.end(),messageUid));
        }
        RecvRep(messageContent);
      }
    }
  }
}

bool GeantEventDispatcher::ResendMsg() {
  bool ok = true;
  for(auto it = fPendingRequests.begin(); it != fPendingRequests.end();){
    auto& pend = it->second;
    bool remove = false;
    if(fHPCWorkers.count(pend.workerId) > 0) {
      if (pend.lastRetry.Since() > std::chrono::seconds(5)) {
        if (pend.retries >= 1) {
          remove = true;
        } else {
          std::cout << "Resending msg: " << pend.msg << std::endl;
          SendMessage(pend.msg, "REQ", it->first, fHPCWorkers[pend.workerId].fZMQID);
          pend.retries++;
          pend.lastRetry.Update();
        }
      }
    }else{
      remove = true; //Worker died - removing request
    }
    if(remove) {
      std::cout << "Removing msg: " << pend.msg << std::endl;
      it = fPendingRequests.erase(it);
      fHPCWorkers[pend.workerId].fDiscardedMsg++;
      ok = false;
    }
    else {
      it++;
    }
  }
  return ok;
}

void GeantEventDispatcher::FinishWorkers() {
  for(auto& it : fHPCWorkers){
    auto& w = it.second;
    SendFinishMsg(w);
  }
  while(fPendingRequests.size() > 0){
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
  req["wrk_id"] = job.fWorkerID;
  req["job_id"] = job.fUID;
  if(fPendingJobs.count(job.fUID)) {
    if(retToPool) {
      fJobPool->ReturnToPool(fPendingJobs[job.fUID]);
    }
    fPendingJobs.erase(job.fUID);
  }
  SendReq(req.dump(),fHPCWorkers[job.fWorkerID]);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
