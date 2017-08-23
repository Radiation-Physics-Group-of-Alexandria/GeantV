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
#include "zmq_util.h"

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
  void RunCommunicationThread();

  void EventAdded();
  void EventTransported(int evt);

  bool GetIsTransportCompleted(){ return isTransportCompleted;}

private:
  zmq::context_t zmqContext;

  std::string fServHname;
  int fServPort;

  std::string fWorkerHname;
  int fWorkerPid;

  GeantConfig *config;
  GeantRunManager *runManager;

  std::atomic_int eventDiff;
  int fFetchAhead;
  bool isTransportCompleted;
  int fReceivedEvents = 0;

  std::map<int, HPCJob> jobs;
  std::mutex jobsMutex;
  ZmqTimer lastContact;
  ZmqTimer lastJobAsk;
  ZmqTimer lastMasterAsk;


  zmq::socket_t fSocket;
  bool connected;
  int connectTries;
  int discardedMsg;
  void SendJobConfirm(int id);
  void SendMasterIntro();
  void RecvMasterIntro(const json& msg);
  void SendJobRequest(int num);
  void RecvJobRequest(const json& msg);
  void SendHB();

  zmq_pollitem_t zmqPollItem;
  void BindSocket();
  void DisconnectFromMaster();
  void SendMessage(const std::string& msg, const std::string& type, size_t uid);
  void SendReq(const std::string& msg);
  void SendRep(const std::string& msg, size_t uid);
  string RecvReq(std::string &msg);
  void RecvRep(std::string& msg);
  void PollForMsg();
  bool ResendMsg();

  std::map<size_t, PendingMessage> pendingRequests;

  json HandleFinishMsg(json &msg);
  json HandleJobCancelMsg(json &msg);
  int wrk_id;
};
}
}

#endif
