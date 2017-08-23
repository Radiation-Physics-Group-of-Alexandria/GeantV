#include "GeantEventReceiver.h"

#include <iostream>

#include "Geant/Error.h"
#include "HepMCGeneratorMultFiles.h"
#include "GeantJobPool.h"


namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantEventReceiver::GeantEventReceiver(const std::string &serverHostName, const std::string &workerHostName, GeantRunManager *runmgr,
                                       GeantConfig *conf)
    : fZMQContext(1), fSocket(fZMQContext, ZMQ_DEALER), fServHname(serverHostName), fServPort(conf->fMasterPort),
      fWorkerHname(workerHostName), fWorkerPid(getpid()),
      fGeantConfig(conf), fRunManager(runmgr),
      fEventDiff(0), fFetchAhead(conf->fNbuff + 2),fIsTransportCompleted(false),
      fReceivedEvents(0), fConnected(false), fWorkerID(-1), fConnectRetries(0),
      fDiscardedMsgs(0)
{
}

void GeantEventReceiver::BindSocket() {
  auto zmqId = fWorkerHname + std::to_string(fWorkerPid);
  fSocket.setsockopt(ZMQ_IDENTITY,zmqId.c_str(),zmqId.size());
  auto masterZMQAddress = std::string("tcp://") + fServHname + ":" + std::to_string(fServPort);
  fSocket.connect(masterZMQAddress);
  fZMQSocketPollItem = {fSocket, 0, ZMQ_POLLIN, 0};
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
    while(fRunManager->GetEventServer() == nullptr){
      sleep(1);
    }
    fLastContact.Update();
    RunCommunicationThread();
  });
  fRunManager->RunSimulation();
  commThread.join();
}

void GeantEventReceiver::RunCommunicationThread() {
  while(!fIsTransportCompleted){
    if(fConnected){
      int diff = fEventDiff.load();
      if(diff < fFetchAhead){
        if(fLastJobAsk.Since() > std::chrono::milliseconds(5*1000)){
          SendJobRequest(fFetchAhead - diff);
          fLastJobAsk.Update();
        }
      }

      {
        std::lock_guard<std::mutex> lock_guard(fJobsMutex);
        for(auto it = fJobs.begin(); it != fJobs.end();){
          auto& job = it->second;
          if (job.left == 0){
            SendJobConfirm(job.id);
            it = fJobs.erase(it);
          }else{
            ++it;
          }
        }
      }
    }else {
      if(fLastMasterAsk.Since() > std::chrono::seconds(5)){
        if(fConnectRetries >= 3) {
          fIsTransportCompleted = true;
          return;
        }
        SendMasterIntro();
        fLastMasterAsk.Update();
      }
    }
    if (fLastContact.Since() > std::chrono::seconds(20)){
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
  fEventDiff.fetch_add(1);
}

void GeantEventReceiver::EventTransported(int evt) {
  fEventDiff.fetch_add(-1);
  std::lock_guard<std::mutex> lock_guard(fJobsMutex);
  for(auto& job_it: fJobs){
    auto& job = job_it.second;
    if(job.startID <= evt && evt < job.endID){
      --job.left;
      break;
    }
  }
}

void GeantEventReceiver::PollForMsg() {
  zmq::poll(&fZMQSocketPollItem,1,std::chrono::milliseconds(10));
  if (fZMQSocketPollItem.revents & ZMQ_POLLIN) {
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
      if(fPendingRequests.count(messageUid) > 0){
        fPendingRequests.erase(messageUid);
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

  fLastContact.Update();

  std::cout << "send: " << type << '\n';
  std::cout << msg << '\n';
}

void GeantEventReceiver::SendReq(const std::string &msg) {
  if(fPendingRequests.size() >= 3)
    return;
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  size_t msgUid = std::hash<std::string>{}(msg+std::to_string(time));
  fPendingRequests[msgUid] = msg;
  SendMessage(msg,"REQ",msgUid);
}

void GeantEventReceiver::SendRep(const std::string &msg,size_t uid) {
  SendMessage(msg,"REP",uid);
}

bool GeantEventReceiver::ResendMsg() {
  bool ok = true;
  for(auto it = fPendingRequests.begin(); it != fPendingRequests.end();){
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
      it = fPendingRequests.erase(it);
      fDiscardedMsgs++;
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
  auto msgType = req["type"].get<std::string>();
  if (msgType == "finish_req") {
    return HandleFinishMsg(req).dump();
  } else if (msgType == "job_cancel_req"){
    return HandleJobCancelMsg(req).dump();
  } else if (msgType == "stat_load_req"){
    return HandleGetLoadMsg(req).dump();
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
  fConnectRetries++;
  SendReq(req.dump());
}

void GeantEventReceiver::RecvMasterIntro(const json &msg) {
    fWorkerID = msg["wrk_id"].get<int>();
    fConnected = true;
    fConnectRetries = 0;
    fDiscardedMsgs = 0;
}

void GeantEventReceiver::SendJobRequest(int num) {
  json req;
  req["type"] = "job_req";
  req["wrk_id"] = fWorkerID;
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
    new_job.startID = fRunManager->GetEventServer()->GetNload();
    new_job.endID = new_job.startID + recv_jobs.size();
    {
      std::lock_guard<std::mutex> lock_guard(fJobsMutex);
      fJobs[new_job.id] = new_job;
    }
    for (auto &job : recv_jobs) {
      auto multFileGenerator = (HepMCGeneratorMultFiles *) fRunManager->GetPrimaryGenerator();
      multFileGenerator->SetEventSource(job.fFilename, job.fOffset);
      for (int j = 0; j < job.fAmount; ++j) {
        fRunManager->GetEventServer()->AddEvent();
        ++fReceivedEvents;
      }
    }
  }
  if (msg["job_type"].get<std::string>() == "generate") {
    int events = msg["event"].get<int>();
    fReceivedEvents = 0;
    new_job.left = events;
    new_job.startID = fRunManager->GetEventServer()->GetNload();
    new_job.endID = new_job.startID + events;
    {
      std::lock_guard<std::mutex> lock_guard(fJobsMutex);
      fJobs[new_job.id] = new_job;
    }
    for (int i = 0; i < events; ++i) {
      fRunManager->GetEventServer()->AddEvent();
      ++fReceivedEvents;
    }
  }
  fRunManager->GetEventServer()->ActivateEvents();
  std::cout << "wrk recv events: " << fReceivedEvents << std::endl;
}

void GeantEventReceiver::SendHB(){
  if(fConnected){
    json hb;
    hb["type"] = "hb";
    hb["wrk_id"] = fWorkerID;
    SendReq(hb.dump());
  }
}

void GeantEventReceiver::SendJobConfirm(int id){
  json req;
  req["type"] = "job_done_req";
  req["wrk_id"] = fWorkerID;
  req["job_id"] = id;
  SendReq(req.dump());
}

void GeantEventReceiver::DisconnectFromMaster() {
  std::cout << "Disconnecting from master" << std::endl;
  fConnected = false;
  fPendingRequests.clear();
  fWorkerID = -1;
  fJobs.clear();
}

json GeantEventReceiver::HandleFinishMsg(json &msg) {
  json rep;
  rep["type"] = "finish_rep";
  DisconnectFromMaster();
  fIsTransportCompleted = true;
  return rep;
}

json GeantEventReceiver::HandleJobCancelMsg(json &msg) {
  std::lock_guard<std::mutex> lock_guard(fJobsMutex);
  int jobId = msg["job_id"].get<int>();
  if(fJobs.count(jobId) > 0)
    fJobs.erase(jobId);

  json rep;
  rep["type"] = "job_cancel_rep";
  rep["job_id"] = jobId;
  return rep;
}

json GeantEventReceiver::HandleGetLoadMsg(json &msg) {
  json rep;
  rep["type"] = "stat_load_rep";
  double loadAvg1min = 0.0;
  getloadavg(&loadAvg1min,1);
  rep["load"] = loadAvg1min;
  return rep;
}

}
}
