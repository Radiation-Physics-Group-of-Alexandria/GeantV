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
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REQ), fServHname(serverHostName), config(conf),
      runManager(runmgr), isTransportCompleted(false), connected(false),  zmqSocketIn(zmqContext,ZMQ_REP),
      connectTries(0), fWorkerPort(conf->fWorkerPort), fWorkerHname(workerHostName),fMasterPort(conf->fMasterPort),
      eventDiff(0)
{
  fFetchAhead = conf->fNbuff + 2;
}

void GeantEventReceiver::BindSocket() {
  zmqSocket.setsockopt(ZMQ_RCVTIMEO,5*1000); //5 sec
  auto masterZMQAddress = std::string("tcp://") + fServHname + ":" + std::to_string(fMasterPort);
  zmqSocket.connect(masterZMQAddress);
  std::cout << "HPC: Event receiver connected to: " << masterZMQAddress << std::endl;
}

//______________________________________________________________________________
void GeantEventReceiver::Initialize()
{

  BindSocket();
  zmqSocketIn.bind("tcp://*:"+std::to_string(fWorkerPort));
  std::cout << "Worker bounded to " << "tcp://*:"+std::to_string(fWorkerPort) << '\n';

  fReceivedEvents = 0;
}

void GeantEventReceiver::Run()
{
  std::thread commThread([this]{
    while(runManager->GetEventServer() == nullptr){
      sleep(1);
    }
    lastContact = std::chrono::system_clock::now();
    RunCommunicationThread();
  });
  runManager->RunSimulation();
  commThread.join();
}

bool GeantEventReceiver::SendMsg(const std::string& req, std::string& rep) {
  zmq::message_t request(req.size());
  memcpy(request.data(),req.c_str(),req.size());
  std::cout << "wrk: msg req " << req << std::endl;

  zmqSocket.send(request);
  zmq::message_t reply;
  bool ok = zmqSocket.recv(&reply);
  if (!ok) return false;

  std::string rep_msg((char* )reply.data(),(char* )reply.data() + reply.size());
  rep = rep_msg;
  std::cout << "wrk: msg rep" << rep << std::endl;

  lastContact = std::chrono::system_clock::now();
  return true;
}

bool GeantEventReceiver::ConnectToMaster(){
  if (connectTries >= MaxReconnectTry){
    return false;
  }
  std::cout << "wrk: connect to master" << std::endl;
  connectTries++;
  json req;
  req["type"] = "new_wrk_req";
  auto workerZMQAddress = std::string("tcp://") + fWorkerHname + ":" + std::to_string(fWorkerPort);
  req["adr"] = workerZMQAddress;
  std::string msg = req.dump();
  std::string rep_msg;
  bool ok = SendMsg(msg,rep_msg);
  if (!ok) return false;
  auto rep = json::parse(rep_msg);
  if (rep["type"].get<std::string>() == "new_wrk_rep"){
    wrk_id = rep["wrk_id"].get<int>();
    connected = true;
    connectTries = 0;
    return true;
  } else {
    return false;
  }
}

int GeantEventReceiver::RequestJob(int num){
  std::cout << "wrk: req job" << std::endl;

  json req;
  req["type"] = "job_req";
  req["wrk_id"] = wrk_id;
  req["events"] = num;
  std::string msg = req.dump();
  std::string rep_msg;
  bool ok = SendMsg(msg,rep_msg);
  if (!ok) return 0;

  auto rep = json::parse(rep_msg);
  if (rep["type"].get<std::string>() == "job_rep"){
    //wrk_id = rep["wrk_id"].get<int>();
    //connected = true;
    HPCJob new_job;
    new_job.id = rep["job_id"].get<int>();
    if(rep["job_type"].get<std::string>() == "hepmc") {
      std::vector<GeantHepMCJob> recv_jobs = rep["files"].get<std::vector<GeantHepMCJob>>();
      fReceivedEvents = 0;
      new_job.left = recv_jobs.size();
      new_job.startID = runManager->GetEventServer()->GetNload();
      new_job.endID = new_job.startID + recv_jobs.size();
      {
        std::lock_guard<std::mutex> lock_guard(job_mutex);
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
    if(rep["job_type"].get<std::string>() == "generate") {
      int events = rep["event"].get<int>();
      fReceivedEvents = 0;
      new_job.left = events;
      new_job.startID = runManager->GetEventServer()->GetNload();
      new_job.endID = new_job.startID + events;
      {
        std::lock_guard<std::mutex> lock_guard(job_mutex);
        jobs[new_job.id] = new_job;
      }
      for (int i = 0; i < events; ++i) {
        runManager->GetEventServer()->AddEvent();
        ++fReceivedEvents;
      }
    }
    runManager->GetEventServer()->ActivateEvents();
    std::cout << "wrk recv events: " << fReceivedEvents << std::endl;
    return 1;
  } else if(rep["type"].get<std::string>() == "wait"){
    return -1;
  } else{
    return 0;
  }
}

bool GeantEventReceiver::ConfirmJob(int id){
  std::cout << "wrk: confirm job" << std::endl;
  json req;
  req["type"] = "job_done_req";
  req["wrk_id"] = wrk_id;
  req["job_id"] = id;
  std::string msg = req.dump();
  std::string rep_msg;
  bool ok = SendMsg(msg,rep_msg);
  if (!ok) return false;
  return true;
}

void GeantEventReceiver::Reconnect() {
  connected = false;
  wrk_id = -1;
  zmqSocket.close();
  zmqSocket = zmq::socket_t(zmqContext,ZMQ_REQ);
  BindSocket();
}

void GeantEventReceiver::SendHB(){
  std::lock_guard<std::mutex> lockGuard(zmqMutex);
    if(connected){
      json hb;
      hb["type"] = "hb";
      hb["wrk_id"] = wrk_id;
      std::string req = hb.dump();
      std::string rep;
      bool ok = SendMsg(req,rep);
      if (!ok) Reconnect();
    }
}

int GeantEventReceiver::AskForNewEvent(int num)
{
  std::lock_guard<std::mutex> lockGuard(zmqMutex);
  std::cout << "wrk ask new event"<<std::endl;
  reconn:
  bool ok = true;
  if(!connected){
    ok = ConnectToMaster();
  }
  if (!ok){
    std::cout << "wrk connect error" << std::endl;
    if(connectTries >= MaxReconnectTry){
      isTransportCompleted = true;
      return 0;
    }
    Reconnect();
    goto reconn;
  }

  int code = RequestJob(num);
  if (code == 0){
    std::cout << "wrk request job error" << std::endl;
    Reconnect();
    goto reconn;
  }
  return fReceivedEvents;
}
void GeantEventReceiver::RunCommunicationThread() {
  while(!isTransportCompleted){
    sleep(3);
    int diff = eventDiff.load();
    if(diff < fFetchAhead){
      AskForNewEvent(fFetchAhead - diff);
    }
    {
      std::lock_guard<std::mutex> lock_guard(job_mutex);
      for(auto it = jobs.begin(); it != jobs.end();){
        auto& job = it->second;
        if (job.left == 0){
          ConfirmJob(job.id);
          it = jobs.erase(it);
        }else{
          ++it;
        }

      }
    }
    auto now = std::chrono::system_clock::now();
    auto since = std::chrono::duration_cast<std::chrono::seconds>(now - lastContact);
    if (since > std::chrono::seconds(20)){
      SendHB();
    }

    ReceiveFromMaster();

  }
}

void GeantEventReceiver::EventAdded() {
  eventDiff.fetch_add(1);
}

void GeantEventReceiver::EventTransported(int evt) {
  eventDiff.fetch_add(-1);
  std::lock_guard<std::mutex> lock_guard(job_mutex);
  for(auto& job_it: jobs){
    auto& job = job_it.second;
    if(job.startID <= evt && evt < job.endID){
      --job.left;
      break;
    }
  }
}

void GeantEventReceiver::ReceiveFromMaster() {
 zmq_pollitem_t items [] = { { zmqSocketIn, 0, ZMQ_POLLIN, 0 } };
 zmq::poll(items,1,std::chrono::milliseconds(30));
    if (items[0].revents & ZMQ_POLLIN) {

      zmq::message_t request;
      zmqSocketIn.recv(&request);
      std::string requestMsg((char* )request.data(),(char* )request.data()+request.size());
      std::cout << "worker from mast req: " << requestMsg << std::endl;
      json req = json::parse(requestMsg);
      auto req_type = req["type"].get<std::string>();

      json rep;
      if (req_type=="finish"){
        rep = ReceiveFinish(req);
      } else if(req_type =="cancel"){
        rep = ReceiveCancel(req);
      }

      std::string responce = rep.dump();
      std::cout << "worker to master rep: " << responce << std::endl;

      zmq::message_t reply(responce.size());
      memcpy(reply.data(), responce.c_str(), responce.size());
      zmqSocketIn.send(reply);
  }
}

json GeantEventReceiver::ReceiveFinish(json& msg) {
  isTransportCompleted = true;

  json j;
  j["type"] = "finish_rep";
  return j;
}

json GeantEventReceiver::ReceiveCancel(json& msg) {
  std::lock_guard<std::mutex> lock_guard(job_mutex);

  json j;
  j["type"] = "cancel_rep";

  int job_id = msg["job_id"].get<int>();
  int msg_wrk_id = msg["wrk_id"].get<int>();

  if(msg_wrk_id != wrk_id) return j;
  if(jobs.count(job_id) == 0) return j;

  jobs.erase(job_id);

  return j;
}

}
}
