#include "GeantEventReceiver.h"

#include "Geant/Error.h"

#include "GeantEvent.h"
#include "GeantEventServer.h"
#include "GeantEventDispatcher.h"
#include "GeantRunManager.h"
#include "PrimaryGenerator.h"
#include "MCTruthMgr.h"

#ifdef USE_HPC
#include "zmq.hpp"
#include <json.hpp>
using nlohmann::json;
#endif

#include <iostream>
#include <cstdint>
#include <cstring>
#include <HepMCGeneratorMultFiles.h>
#include <GeantJobPool.h>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

int MaxReconnectTry = 3;
//______________________________________________________________________________
GeantEventReceiver::GeantEventReceiver(const std::string &serverHostName, const std::string &workerHostName, GeantRunManager *runmgr,
                                       GeantConfig *conf)
    : zmqContext(1), fSocket(zmqContext, ZMQ_DEALER), fServHname(serverHostName), config(conf),
      runManager(runmgr), isTransportCompleted(false), connected(false),
      connectTries(0),fWorkerHname(workerHostName),fServPort(conf->fMasterPort),
      eventDiff(0), fWorkerPid(getpid())
{
  fFetchAhead = conf->fNbuff + 2;
}

void GeantEventReceiver::BindSocket() {
  auto zmqId = fWorkerHname + std::to_string(fWorkerPid);
  fSocket.setsockopt(ZMQ_IDENTITY,zmqId.c_str(),zmqId.size());
  auto masterZMQAddress = std::string("tcp://") + fServHname + ":" + std::to_string(fServPort);
  fSocket.connect(masterZMQAddress);
  zmqPollItem = {fSocket, 0, ZMQ_POLLIN, 0};
  std::cout << "HPC: Event receiver connected to: " << masterZMQAddress << std::endl;
}

//______________________________________________________________________________
void GeantEventReceiver::Initialize()
{

  BindSocket();
  fReceivedEvents = 0;
}

void GeantEventReceiver::Run()
{
  std::thread commThread([this]{
    while(runManager->GetEventServer() == nullptr){
      sleep(1);
    }
    lastContact.Update();
    RunCommunicationThread();
  });
  runManager->RunSimulation();
  commThread.join();
}

void GeantEventReceiver::RunCommunicationThread() {
  while(!isTransportCompleted){
    if(connected){
      int diff = eventDiff.load();
      if(diff < fFetchAhead){
        if(lastJobAsk.Since() > std::chrono::milliseconds(5*1000)){
          SendJobRequest(fFetchAhead - diff);
          lastJobAsk.Update();
        }
      }

      {
        std::lock_guard<std::mutex> lock_guard(jobsMutex);
        for(auto it = jobs.begin(); it != jobs.end();){
          auto& job = it->second;
          if (job.left == 0){
            SendJobConfirm(job.id);
            it = jobs.erase(it);
          }else{
            ++it;
          }
        }
      }
    }else {
      if(lastMasterAsk.Since() > std::chrono::seconds(5)){
        if(connectTries >= 3) {
          isTransportCompleted = true;
          return;
        }
        SendMasterIntro();
        lastMasterAsk.Update();
      }
    }
    if (lastContact.Since() > std::chrono::seconds(20)){
      SendHB();
    }

    PollForMsg();
    bool ok = ResendMsg();
    if(!ok){
      DisconnectFromMaster();
    }

  }
}

void GeantEventReceiver::EventAdded() {
  eventDiff.fetch_add(1);
}

void GeantEventReceiver::EventTransported(int evt) {
  eventDiff.fetch_add(-1);
  std::lock_guard<std::mutex> lock_guard(jobsMutex);
  for(auto& job_it: jobs){
    auto& job = job_it.second;
    if(job.startID <= evt && evt < job.endID){
      --job.left;
      break;
    }
  }
}

void GeantEventReceiver::PollForMsg() {
  zmq::poll(&zmqPollItem,1,std::chrono::milliseconds(10));
  if (zmqPollItem.revents & ZMQ_POLLIN) {
    zmq::message_t type;
    zmq::message_t muid;
    zmq::message_t message;
    fSocket.recv(&type);
    fSocket.recv(&muid);
    fSocket.recv(&message);
    std::string messageType((char*)type.data(),(char*)type.data() + type.size());
    size_t messageUid;
    memcpy(&messageUid,muid.data(),muid.size());
    std::string messageContent((char*)message.data(),(char*)message.data() + message.size());
    std::cout << "recv: " << messageType << '\n';
    std::cout << messageContent << '\n';
    if(messageType == "REQ"){
      std::string replyContent = RecvReq(messageContent);
      SendRep(replyContent,messageUid);
    }else { //messageType == REP
      if(pendingRequests.count(messageUid) > 0){
        pendingRequests.erase(messageUid);
        RecvRep(messageContent);
      }
    }
  }
}

void GeantEventReceiver::SendMessage(const std::string &msg, const std::string &type, size_t uid) {
  zmq::message_t m_type(type.size());
  memcpy(m_type.data(),type.c_str(),type.size());
  zmq::message_t m_uid(sizeof(size_t));
  memcpy(m_uid.data(),&uid,sizeof(size_t));
  zmq::message_t message(msg.size());
  memcpy(message.data(),msg.c_str(),msg.size());

  fSocket.send(m_type,ZMQ_SNDMORE);
  fSocket.send(m_uid,ZMQ_SNDMORE);
  fSocket.send(message);

  lastContact.Update();

  std::cout << "send: " << type << '\n';
  std::cout << msg << '\n';
}

void GeantEventReceiver::SendReq(const std::string &msg) {
  if(pendingRequests.size() >= 3)
    return;
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  size_t msgUid = std::hash<std::string>{}(msg+std::to_string(time));
  pendingRequests[msgUid] = msg;
  SendMessage(msg,"REQ",msgUid);
}

void GeantEventReceiver::SendRep(const std::string &msg,size_t uid) {
  SendMessage(msg,"REP",uid);
}

bool GeantEventReceiver::ResendMsg() {
  bool ok = true;
  for(auto it = pendingRequests.begin(); it != pendingRequests.end();){
    auto& pend = it->second;
    bool remove = false;
    if(pend.lastRetry.Since() > std::chrono::seconds(5)){
      if(pend.retries >= 1){
        remove = true;
      }else {
        std::cout << "Resending msg: " << pend.msg << std::endl;
        SendMessage(pend.msg, "REQ", it->first);
        pend.retries++;
        pend.lastRetry.Update();
      }
    }
    if(remove) {
      std::cout << "Removing msg: " << pend.msg << std::endl;
      it = pendingRequests.erase(it);
      discardedMsg++;
      ok = false;
    }
    else {
      it++;
    }
  }
  return ok;
}

string GeantEventReceiver::RecvReq(std::string &msg) {
  json req = json::parse(msg);
  if (req["type"].get<std::string>() == "finish_req") {
    return HandleFinishMsg(req).dump();
  } else if (req["type"].get<std::string>() == "job_cancel_req"){
    return HandleJobCancelMsg(req).dump();
  }
  return "{}";
}

void GeantEventReceiver::RecvRep(std::string &msg) {
  json rep = json::parse(msg);
  if (rep["type"].get<std::string>() == "job_rep") {
    RecvJobRequest(rep);
  } else if (rep["type"].get<std::string>() == "new_wrk_rep"){
    RecvMasterIntro(rep);
  } else if (rep["type"].get<std::string>() == "error") {
    DisconnectFromMaster();
  }
}

void GeantEventReceiver::SendMasterIntro() {
  json req;
  req["type"] = "new_wrk_req";
  req["zmq_id"] =  fWorkerHname + std::to_string(fWorkerPid);
  connectTries++;
  SendReq(req.dump());
}

void GeantEventReceiver::RecvMasterIntro(const json &msg) {
    wrk_id = msg["wrk_id"].get<int>();
    connected = true;
    connectTries = 0;
    discardedMsg = 0;
}

void GeantEventReceiver::SendJobRequest(int num) {
  json req;
  req["type"] = "job_req";
  req["wrk_id"] = wrk_id;
  req["events"] = num;
  SendReq(req.dump());
}

void GeantEventReceiver::RecvJobRequest(const json &msg) {
  HPCJob new_job;
  new_job.id = msg["job_id"].get<int>();
  if (msg["job_type"].get<std::string>() == "hepmc") {
    std::vector<GeantHepMCJob> recv_jobs = msg["files"].get<std::vector<GeantHepMCJob>>();
    fReceivedEvents = 0;
    new_job.left = recv_jobs.size();
    new_job.startID = runManager->GetEventServer()->GetNload();
    new_job.endID = new_job.startID + recv_jobs.size();
    {
      std::lock_guard<std::mutex> lock_guard(jobsMutex);
      jobs[new_job.id] = new_job;
    }
    for (auto &job : recv_jobs) {
      auto multFileGenerator = (HepMCGeneratorMultFiles *) runManager->GetPrimaryGenerator();
      multFileGenerator->SetEventSource(job.filename, job.offset);
      for (int j = 0; j < job.amount; ++j) {
        runManager->GetEventServer()->AddEvent();
        ++fReceivedEvents;
      }
    }
  }
  if (msg["job_type"].get<std::string>() == "generate") {
    int events = msg["event"].get<int>();
    fReceivedEvents = 0;
    new_job.left = events;
    new_job.startID = runManager->GetEventServer()->GetNload();
    new_job.endID = new_job.startID + events;
    {
      std::lock_guard<std::mutex> lock_guard(jobsMutex);
      jobs[new_job.id] = new_job;
    }
    for (int i = 0; i < events; ++i) {
      runManager->GetEventServer()->AddEvent();
      ++fReceivedEvents;
    }
  }
  runManager->GetEventServer()->ActivateEvents();
  std::cout << "wrk recv events: " << fReceivedEvents << std::endl;
}

void GeantEventReceiver::SendHB(){
  if(connected){
    json hb;
    hb["type"] = "hb";
    hb["wrk_id"] = wrk_id;
    SendReq(hb.dump());
  }
}

void GeantEventReceiver::SendJobConfirm(int id){
  json req;
  req["type"] = "job_done_req";
  req["wrk_id"] = wrk_id;
  req["job_id"] = id;
  SendReq(req.dump());
}

void GeantEventReceiver::DisconnectFromMaster() {
  std::cout << "Disconnecting from master" << std::endl;
  connected = false;
  pendingRequests.clear();
  wrk_id = -1;
  jobs.clear();
}

json GeantEventReceiver::HandleFinishMsg(json &msg) {
  json rep;
  rep["type"] = "finish_rep";
  DisconnectFromMaster();
  isTransportCompleted = true;
  return rep;
}

json GeantEventReceiver::HandleJobCancelMsg(json &msg) {
  std::lock_guard<std::mutex> lock_guard(jobsMutex);
  int jobId = msg["job_id"].get<int>();
  if(jobs.count(jobId) > 0)
    jobs.erase(jobId);

  json rep;
  rep["type"] = "job_cancel_rep";
  rep["job_id"] = jobId;
  return rep;
}

}
}
