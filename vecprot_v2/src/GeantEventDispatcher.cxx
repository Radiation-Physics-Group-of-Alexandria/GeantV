#include "GeantEventDispatcher.h"

#include <iostream>
#include <fstream>

#include "globals.h"
#include "Geant/Error.h"

#include "GeantJobPool.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GeantEventDispatcher::GeantEventDispatcher(GeantConfig *config)
    : fZMQContext(1), fSocket(fZMQContext, ZMQ_ROUTER), fConfig(config), fJobPool(config->jobPool), fWorkerCounter(0)
{
}

const json kErrorReply = "{\"type\":\"error\"}"_json;

json GeantEventDispatcher::HandleNewWorkerMsg(const std::string &address)
{
  fWorkerCounter++;
  GeantHPCWorker worker;
  worker.fID    = fWorkerCounter;
  worker.fZMQID = address;
  worker.fLastContact.Update();
  if (fHPCWorkerZMQIDtoUID[address] > 0 && fHPCWorkers.count(fHPCWorkerZMQIDtoUID[address]) > 0) {
    CleanWorkerState(fHPCWorkers[fHPCWorkerZMQIDtoUID[address]]);
    fHPCWorkers.erase(fHPCWorkerZMQIDtoUID[address]);
    fHPCHosts[address].fReconnects++;
  }
  if (fHPCHosts.count(address) == 0) {
    std::cout << "New host added from: " << address << '\n';
    fHPCHosts[address] = Host(address);
  }
  fHPCWorkers[worker.fID]             = worker;
  fHPCWorkerZMQIDtoUID[worker.fZMQID] = worker.fID;
  json reply;
  reply["type"]   = "new_wrk_rep";
  reply["wrk_id"] = worker.fID;
  std::cout << "New worker: " << worker.fID << " at: " << address << '\n';
  return reply;
}

json GeantEventDispatcher::HandleNewJobReq(json &req)
{
  int n      = req["events"];
  int wrk_id = req["wrk_id"];
  if (fHPCWorkers.count(wrk_id) == 0) {
    return kErrorReply;
  }
  auto &worker = fHPCWorkers[wrk_id];
  worker.fLastContact.Update();
  GeantHPCJob job = fJobPool->GetJob(n, worker);
  if (job.fEvents == 0) {
    bool foundDupeJob = FindDuplicateJob(worker, job, n);
    if (!foundDupeJob) return "{\"type\":\"wait\"}"_json;
  }
  job.fWorkerID = worker.fID;
  job.fDispatchTime.Update();
  fPendingJobs[job.fUID] = job;
  fHPCHosts[worker.fZMQID].fEventsGiven += job.fEvents;
  json reply;
  reply["type"]   = "job_rep";
  reply["wrk_id"] = wrk_id;
  reply["job_id"] = job.fUID;
  if (job.fType == JobType::HepMC) {
    reply["job_type"] = "hepmc";
    reply["files"]    = job.fHepMCJobs;
  } else {
    reply["job_type"] = "generate";
  }
  reply["event"] = job.fEvents;

  std::cout << "Giving " << job.fEvents << " to worker: " << job.fWorkerID << '\n';

  return reply;
}

bool GeantEventDispatcher::FindDuplicateJob(GeantHPCWorker &worker, GeantHPCJob &job, int n)
{
  bool foundDupeJob = false;
  for (auto pendIt = fPendingJobs.rbegin(); pendIt != fPendingJobs.rend(); pendIt++) {
    auto &pendJob = pendIt->second;
    if (n < pendJob.fEvents) continue;
    if (IsWorkerDoingJob(worker, pendJob)) continue;
    if (*pendJob.fDuplicateTries >= fConfig->fJobMaxDuplicateTries) continue;
    job          = fJobPool->GetDublicateJob(pendJob);
    foundDupeJob = true;
    *job.fDuplicateTries += 1;
    std::cout << "Duplicating job: " << pendJob.fUID << '\n';
    std::cout << "Dupe job id: ";
    for (auto dum : *job.fDublicateUIDs) {
      std::cout << dum << ' ';
    }
    std::cout << '\n';
    break;
  }
  return foundDupeJob;
}

bool GeantEventDispatcher::IsWorkerDoingJob(GeantHPCWorker &worker, GeantHPCJob &job)
{
  if (job.fWorkerID == worker.fID) return true;
  for (int jobUid : *job.fDublicateUIDs) {
    if (fPendingJobs.count(jobUid) > 0) {
      if (fPendingJobs[jobUid].fWorkerID == worker.fID) return true;
    }
  }
  return false;
}

json GeantEventDispatcher::HandleJobConfirm(json &req)
{
  int job_id = req["job_id"];
  int wrk_id = req["wrk_id"];
  if (fHPCWorkers.count(wrk_id) == 0) {
    return kErrorReply;
  }
  auto &worker = fHPCWorkers[wrk_id];

  worker.fLastContact.Update();

  if (fPendingJobs.count(job_id) == 0) return kErrorReply;
  auto &job = fPendingJobs[job_id];

  worker.fExpectedTimeForEvent =
      (worker.fTransportedEvents) * (worker.fExpectedTimeForEvent) + job.fDispatchTime.Since();
  worker.fTransportedEvents += job.fEvents;
  worker.fExpectedTimeForEvent /= worker.fTransportedEvents;

  fHPCHosts[worker.fZMQID].fEventsConfirmed += job.fEvents;
  fHPCHosts[worker.fZMQID].fAverageTimeForEvent = worker.fExpectedTimeForEvent;

  job.fDublicateUIDs->erase(job_id);
  for (auto dupeId : *job.fDublicateUIDs) {
    SendJobCancelMsg(fPendingJobs[dupeId], false);
  }
  fPendingJobs.erase(job_id);
  json rep;
  rep["type"] = "job_done_rep";
  std::cout << "Confirming job: " << job_id << std::endl;
  return rep;
}

json GeantEventDispatcher::HandleHeartbeat(json &req)
{
  int wrk_id = req["wrk_id"];
  if (fHPCWorkers.count(wrk_id) == 0) {
    return kErrorReply;
  }
  auto &worker = fHPCWorkers[wrk_id];
  worker.fLastContact.Update();
  return "{\"type\": \"hb_rep\"}"_json;
}

void GeantEventDispatcher::Initialize()
{
  BindSocket();
}

void GeantEventDispatcher::RunReqReplyLoop()
{
  ZmqTimer totalTime;
  totalTime.Update();
  while (!fJobPool->IsEmpty() || !fPendingJobs.empty()) {
    PollForMsg();
    ResendMsg();
    CleanDeadWorkers();
    CleanDeadJobs();
    UpdateWorkerStats();
  }
  auto timeSpent = std::chrono::duration_cast<std::chrono::seconds>(totalTime.Since()).count();
  FinishWorkers();
  std::cout << "Summary: " << '\n';
  for (auto &hostIt : fHPCHosts) {
    auto &host = hostIt.second;
    std::cout << "ZMQAddress: " << host.fZMQAddress << '\n';
    std::cout << "Reconnects: " << host.fReconnects << '\n';
    std::cout << "Events given: " << host.fEventsGiven << '\n';
    std::cout << "Events confirmed: " << host.fEventsConfirmed << '\n';
    std::cout << "Average time for confirming event: " << host.fAverageTimeForEvent.count() << " ms" << '\n';
    std::cout << "Last min load avg: " << host.fLastMinLoadAverage << '\n';
    std::cout << '\n';
  }
  std::cout << "Total time: " << timeSpent << "s" << '\n';
}

void GeantEventDispatcher::CleanDeadWorkers()
{
  for (auto w_it = fHPCWorkers.begin(); w_it != fHPCWorkers.end();) {
    auto w = &w_it->second;
    if (w->fLastContact.Since() > fConfig->fDeadWorkerDetectTime) {
      std::cout << "Worker is dead: " << w->fID << '\n';
      CleanWorkerState(*w);
      w_it = fHPCWorkers.erase(w_it);
    } else {
      ++w_it;
    }
  }
}

void GeantEventDispatcher::BindSocket()
{
  fZMQSocketPollItem = {fSocket, 0, ZMQ_POLLIN, 0};
  fSocket.bind("tcp://*:" + std::to_string(fConfig->fMasterPort));
  std::cout << "HPC: Event dispatcher bounded to tcp://*:" + std::to_string(fConfig->fMasterPort) << std::endl;
}

void GeantEventDispatcher::SendMessage(const std::string &msg, const std::string &type, size_t uid,
                                       const std::string &address, const std::string &localAddress)
{

  zmq::message_t m_addr(address.size());
  memcpy(m_addr.data(), address.c_str(), address.size());
  zmq::message_t m_loc_addr(localAddress.size());
  memcpy(m_loc_addr.data(), localAddress.data(), localAddress.size());
  zmq::message_t m_type(type.size());
  memcpy(m_type.data(), type.c_str(), type.size());
  zmq::message_t m_uid(sizeof(size_t));
  memcpy(m_uid.data(), &uid, sizeof(size_t));
  zmq::message_t message(msg.size());
  memcpy(message.data(), msg.c_str(), msg.size());

  fSocket.send(m_addr, ZMQ_SNDMORE);
  fSocket.send(m_loc_addr, ZMQ_SNDMORE);
  fSocket.send(m_type, ZMQ_SNDMORE);
  fSocket.send(m_uid, ZMQ_SNDMORE);
  fSocket.send(message);

  if (fConfig->fPrintMessages) {
    std::cout << "send: " << address;
    std::cout << " " << type << '\n';
    std::cout << msg << '\n';
  }
}

void GeantEventDispatcher::SendRep(const std::string &msg, size_t uid, const std::string &address,
                                   const std::string &localAddress)
{
  SendMessage(msg, "REP", uid, address, localAddress);
}

std::string GeantEventDispatcher::RecvReqWorker(const std::string &msg, const std::string &address)
{
  json req = json::parse(msg);
  json rep;
  rep["type"] = "error";

  auto req_type = req["type"];
  if (req_type == "new_wrk_req") {
    rep = HandleNewWorkerMsg(address);
  } else if (req_type == "job_req") {
    rep = HandleNewJobReq(req);
  } else if (req_type == "job_done_req") {
    rep = HandleJobConfirm(req);
  } else if (req_type == "hb") {
    rep = HandleHeartbeat(req);
  }
  return rep.dump();
}

void GeantEventDispatcher::SendReqWorker(const std::string &msg, GeantHPCWorker &worker)
{
  size_t msgUid            = CreateMessageUID(msg);
  fPendingRequests[msgUid] = MasterPendingMessage(msg, worker.fZMQID, true);
  worker.fPendingRequests.insert(msgUid);
  SendMessage(msg, "REQ", msgUid, worker.fZMQID, "WRK");
}
void GeantEventDispatcher::SendReqProcMgr(const std::string &msg, const std::string &address)
{
  size_t msgUid            = CreateMessageUID(msg);
  fPendingRequests[msgUid] = MasterPendingMessage(msg, address, false);
  SendMessage(msg, "REQ", msgUid, address, "MGR");
}

void GeantEventDispatcher::RecvRepWorker(const std::string &msg, GeantHPCWorker &worker)
{
  json rep = json::parse(msg);
  if (rep["type"] == "stat_load_rep") {
    RecvGetLoadMsg(rep, worker);
  }
}

void GeantEventDispatcher::PollForMsg()
{
  zmq::poll(&fZMQSocketPollItem, 1, std::chrono::milliseconds(10));
  if (fZMQSocketPollItem.revents & ZMQ_POLLIN) {
    zmq::message_t addr;
    zmq::message_t locAddr;
    zmq::message_t type;
    zmq::message_t muid;
    zmq::message_t message;
    fSocket.recv(&addr);
    fSocket.recv(&locAddr);
    fSocket.recv(&type);
    fSocket.recv(&muid);
    fSocket.recv(&message);
    std::string messageAddr((char *)addr.data(), (char *)addr.data() + addr.size());
    std::string messageLocAddr((char *)locAddr.data(), (char *)locAddr.data() + locAddr.size());
    std::string messageType((char *)type.data(), (char *)type.data() + type.size());
    size_t messageUid;
    memcpy(&messageUid, muid.data(), muid.size());
    std::string messageContent((char *)message.data(), (char *)message.data() + message.size());
    if (fConfig->fPrintMessages) {
      std::cout << "recv: " << messageAddr << " loc: " << messageLocAddr;
      std::cout << " " << messageType << '\n';
      std::cout << messageContent << '\n';
    }
    if (messageType == "REQ") {
      std::string replyContent = RecvReqWorker(messageContent, messageAddr);
      SendRep(replyContent, messageUid, messageAddr, "WRK");
    } else { // messageType == REP
      if (fPendingRequests.count(messageUid) > 0) {
        fPendingRequests.erase(messageUid);
        if (messageLocAddr == "WRK") {
          if (fHPCWorkerZMQIDtoUID.count(messageAddr) > 0 && fHPCWorkers.count(fHPCWorkerZMQIDtoUID[messageAddr]) > 0) {
            auto &worker = fHPCWorkers[fHPCWorkerZMQIDtoUID[messageAddr]];
            worker.fPendingRequests.erase(messageUid);
            RecvRepWorker(messageContent, worker);
          }
        } else { // message loc addr == "MGR"
          RecvRepProcMgr(messageContent, messageAddr);
        }
      }
    }
  }
}

bool GeantEventDispatcher::ResendMsg()
{
  bool ok = true;
  for (auto it = fPendingRequests.begin(); it != fPendingRequests.end();) {
    auto &pend = it->second;

    bool remove = false;
    if (pend.lastRetry.Since() > fConfig->fMessageResendTime) {
      if (pend.retries >= fConfig->fMessageResendRetries) {
        remove = true;
      } else {
        if (fConfig->fPrintMessages) {
          std::cout << "Resending msg: " << pend.msg << std::endl;
        }
        SendMessage(pend.msg, "REQ", it->first, pend.address, pend.toWorker ? "WRK" : "MGR");
        pend.retries++;
        pend.lastRetry.Update();
      }
    }

    if (pend.toWorker) {
      if (fHPCWorkers.count(fHPCWorkerZMQIDtoUID[pend.address]) == 0) {
        remove = true; // Worker died - removing request
      }
    }
    if (remove) {
      if (fConfig->fPrintMessages) {
        std::cout << "Removing msg: " << pend.msg << std::endl;
      }
      it = fPendingRequests.erase(it);
      if (pend.toWorker) {
        if (fHPCWorkers.count(fHPCWorkerZMQIDtoUID[pend.address]) != 0) {
          fHPCWorkers[fHPCWorkerZMQIDtoUID[pend.address]].fDiscardedMsg++;
        }
      }
      ok = false;
    } else {
      it++;
    }
  }
  return ok;
}

void GeantEventDispatcher::FinishWorkers()
{
  for (auto &it : fHPCWorkers) {
    auto &w = it.second;
    SendFinishMsg(w);
  }
  while (fPendingRequests.size() > 0) {
    PollForMsg();
    ResendMsg();
  }
}

void GeantEventDispatcher::SendFinishMsg(GeantHPCWorker &worker)
{
  SendReqWorker("{\"type\":\"finish_req\"}", worker);
}

void GeantEventDispatcher::SendJobCancelMsg(GeantHPCJob &job, bool retToPool)
{
  json req;
  req["type"]   = "job_cancel_req";
  req["wrk_id"] = job.fWorkerID;
  req["job_id"] = job.fUID;
  if (fPendingJobs.count(job.fUID)) {
    if (retToPool) {
      fJobPool->ReturnToPool(fPendingJobs[job.fUID]);
    }
    fPendingJobs.erase(job.fUID);
  }
  SendReqWorker(req.dump(), fHPCWorkers[job.fWorkerID]);
}

void GeantEventDispatcher::CleanDeadJobs()
{
  std::vector<int> jobsForDeleting;
  for (auto &jobIt : fPendingJobs) {
    auto &job = jobIt.second;
    if (job.fDispatchTime.Since() >
        (fConfig->fEventDeadlineMinTime + 2 * fHPCWorkers[job.fWorkerID].fExpectedTimeForEvent) *
            (job.fEvents + job.fHepMCJobs.size())) {
      jobsForDeleting.push_back(job.fUID);
    }
  }
  for (int delJobID : jobsForDeleting) {
    std::cout << "Deleting dead job: " << delJobID << '\n';
    if (fPendingJobs.count(delJobID) > 0) {
      fPendingJobs[delJobID].fDublicateUIDs->erase(delJobID);
      bool returnToPool = fPendingJobs[delJobID].fDublicateUIDs->size() == 0 ? true : false;
      SendJobCancelMsg(fPendingJobs[delJobID], returnToPool);
    }
  }
}

void GeantEventDispatcher::SendGetLoadMsg(GeantHPCWorker &worker)
{
  json req;
  req["type"]   = "stat_load_req";
  req["wrk_id"] = worker.fID;
  worker.fLastStatReq.Update();
  SendReqWorker(req.dump(), worker);
}

void GeantEventDispatcher::RecvGetLoadMsg(const json &msg, GeantHPCWorker &worker)
{
  if (msg["wrk_id"] != worker.fID) return;
  fHPCHosts[worker.fZMQID].fLastMinLoadAverage = msg["load"];
}

void GeantEventDispatcher::UpdateWorkerStats()
{
  for (auto &workerIt : fHPCWorkers) {
    auto &worker = workerIt.second;
    if (worker.fLastStatReq.Since() > std::chrono::seconds(20)) {
      SendGetLoadMsg(worker);
    }
  }
}

void GeantEventDispatcher::SendAbortMsg(const std::string &address)
{
  json req;
  req["type"] = "abort_req";
  SendReqProcMgr(req.dump(), address);
}

void GeantEventDispatcher::SendRestartMsg(const std::string &address)
{
  json req;
  req["type"] = "restart_req";
  SendReqProcMgr(req.dump(), address);
}

void GeantEventDispatcher::SendGetStatusMsg(const std::string &address)
{
  json req;
  req["type"] = "status_req";
  SendReqProcMgr(req.dump(), address);
}

void GeantEventDispatcher::RecvGetStatusMsg(const json &msg, const std::string &address)
{
  std::cout << "Worker " << address << " is running: " << msg["running"].get<bool>() << '\n';
}

void GeantEventDispatcher::CleanWorkerState(GeantHPCWorker &worker)
{
  std::cout << "Cleaning worker: " << worker.fID << " from " << worker.fZMQID << '\n';
  for (auto job_it = fPendingJobs.begin(); job_it != fPendingJobs.end();) {
    if (job_it->second.fWorkerID == worker.fID) {
      job_it->second.fDublicateUIDs->erase(job_it->second.fUID);
      if (job_it->second.fDublicateUIDs->size() == 0) {
        fJobPool->ReturnToPool(job_it->second);
      }
      job_it = fPendingJobs.erase(job_it);
    } else {
      ++job_it;
    }
  }
  if (fHPCWorkerZMQIDtoUID[worker.fZMQID] == worker.fID) fHPCWorkerZMQIDtoUID[worker.fZMQID] = -1;
}

void GeantEventDispatcher::RecvRepProcMgr(const std::string &msg, const std::string &address)
{
  json rep = json::parse(msg);
  if (rep["type"] == "status_rep") {
    RecvGetStatusMsg(rep, address);
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
