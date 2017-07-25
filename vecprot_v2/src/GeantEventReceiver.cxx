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
GeantEventReceiver::GeantEventReceiver(std::string serverHostName, GeantConfig *conf, GeantRunManager *runmgr)
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REQ), fServHname(serverHostName), config(conf),
      runManager(runmgr), isTransportCompleted(false), connected(false), currentJob(-1), zmqSocketIn(zmqContext,ZMQ_REP),
      connectTries(0)
{
  //TODO:
  fWorkerHname = "*:6666";
}

void GeantEventReceiver::BindSocket() {
  zmqSocket.setsockopt(ZMQ_RCVTIMEO,5*1000); //5 sec
  zmqSocket.connect(std::string("tcp://") + fServHname + std::string(":5678"));
  std::cout << "HPC: Event receiver connected to: " << (std::string("tcp://") + fServHname + std::string(":5678"))
            << std::endl;
}

//______________________________________________________________________________
void GeantEventReceiver::Initialize()
{

  BindSocket();
  //zmqSocketIn.bind("tcp://"+fWorkerHname);

  fReceivedEvents = 0;


}

void GeantEventReceiver::Run()
{
  std::thread hbThread([this]{
    while(!isTransportCompleted) {
      sleep(10);
      SendHB();
    }
  });
  runManager->RunSimulation();
  hbThread.join();
}

bool GeantEventReceiver::SendMsg(const std::string& req, std::string& rep) {
  zmq::message_t request(1024);
  memcpy(request.data(),req.c_str(),req.size()+1);
  std::cout << "wrk: msg req " << req << std::endl;

  zmqSocket.send(request);
  zmq::message_t reply;
  bool ok = zmqSocket.recv(&reply);
  if (!ok) return false;

  char rep_msg[1024];
  memcpy(rep_msg, reply.data(), std::min((size_t)1024,reply.size()));
  rep = rep_msg;
  std::cout << "wrk: msg rep" << rep << std::endl;
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
  req["adr"] = "tcp://"+fReceiveAddr;
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
    currentJob = rep["job_id"].get<int>();
    if(rep["job_type"].get<std::string>() == "hepmc") {
      std::vector<GeantHepMCJob> jobs = rep["files"].get<std::vector<GeantHepMCJob>>();
      fReceivedEvents = 0;
      for (auto &job : jobs) {
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
      for (int i = 0; i < events; ++i) {
        runManager->GetEventServer()->AddEvent();
        ++fReceivedEvents;
      }
    }
    runManager->GetEventServer()->ActivateEvents();
    std::cout << "wrk recv events: " << fReceivedEvents << std::endl;
    return 1;
  } else if(rep["type"].get<std::string>() == "wait"){
    sleep(10);
    return -1;
  } else{
    return 0;
  }
}

bool GeantEventReceiver::ConfirmJob(){
  if(currentJob < 0) return true;
  std::cout << "wrk: confirm job" << std::endl;
  json req;
  req["type"] = "job_done_req";
  req["wrk_id"] = wrk_id;
  req["job_id"] = currentJob;
  std::string msg = req.dump();
  std::string rep_msg;
  bool ok = SendMsg(msg,rep_msg);
  if (!ok) return false;
  currentJob = -1;
  return true;
}

void GeantEventReceiver::Reconnect() {
  connected = false;
  currentJob = -1;
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

  ok = ConfirmJob();
  if (!ok){
    std::cout << "wrk confirm job error" << std::endl;
    Reconnect();
    goto reconn;
  }

  bool reaskJob = false;
  int code;
  do {
    code = RequestJob(num);
    if (code == -1) reaskJob = true;
    else reaskJob = false;
  }while(reaskJob);
  if (code == 0){
    std::cout << "wrk request job error" << std::endl;
    Reconnect();
    goto reconn;
  }
  return fReceivedEvents;
}
}
}
