//===--- GeantEventReceiver.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEventReceiver.h
 * @brief Implementation of node event service (based on ZMQ) in Geant-V prototype
 * @author Oksana Shadura, Vitalii Drohan
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_RECEIVER
#define GEANT_RECEIVER
#ifdef USE_HPC
#include <zmq.hpp>
#include <json.hpp>
using nlohmann::json;
#endif

#include <string>

#include "GeantRunManager.h"
#include "GeantConfig.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantEventDispatcher;

struct HPCJob{
  int id;
  int left;
  int startID,endID;
};

class GeantEventReceiver {
public:
  GeantEventReceiver(const std::string &serverHostName, const std::string &workerHostName, GeantRunManager *runmgr,
                       GeantConfig *conf);

  void Initialize();
  void Run();
  int AskForNewEvent(int num);
  void SendHB();
  void RunCommunicationThread();
  void EventAdded();
  void EventTransported(int evt);
  void ReceiveFromMaster();

  bool GetIsTransportCompleted(){ return isTransportCompleted;}

private:
  std::atomic_int eventDiff;
  int fFetchAhead;

  std::mutex zmqMutex;
  zmq::context_t zmqContext;
  zmq::socket_t zmqSocket;
  zmq::socket_t zmqSocketIn;

  std::string fServHname;
  int fMasterPort;
  std::string fWorkerHname;
  int fWorkerPort;
  std::string fReceiveAddr;
  GeantConfig *config;
  GeantRunManager *runManager;

  bool isTransportCompleted;
  int fReceivedEvents = 0;

  std::map<int, HPCJob> jobs;
  std::mutex job_mutex;
  std::chrono::time_point<std::chrono::system_clock> lastContact;

  bool connected;
  int connectTries;
  bool ConnectToMaster();
  bool ConfirmJob(int id);

  int RequestJob(int num);
  void Reconnect();
  void BindSocket();
  bool SendMsg(const std::string& req, std::string& rep);

  json ReceiveFinish(json& msg);
  json ReceiveCancel(json& msg);
  int wrk_id;
};
}
}

#endif
