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
#include <string>
#include <zmq.hpp>
#include <json.hpp>
#include "GeantRunManager.h"
#include "GeantConfig.h"
#include "zmq_util.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using nlohmann::json;
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

  bool GetIsTransportCompleted(){ return fIsTransportCompleted;}

private:
  zmq::context_t fZMQContext;
  zmq::socket_t fSocket;
  zmq_pollitem_t fZMQSocketPollItem;

  std::string fServHname;
  int fServPort;
  std::string fWorkerHname;
  int fWorkerPid;

  GeantConfig *fGeantConfig;
  GeantRunManager *fRunManager;

  std::atomic_int fEventDiff;
  int fFetchAhead;
  bool fIsTransportCompleted;
  int fReceivedEvents;

  std::map<int, HPCJob> fJobs;
  std::mutex fJobsMutex;
  ZmqTimer fLastContact;
  ZmqTimer fLastJobAsk;
  ZmqTimer fLastMasterAsk;

  bool fConnected;
  int fWorkerID;
  int fConnectRetries;
  std::map<size_t, PendingMessage> fPendingRequests;
  int fDiscardedMsgs;

  void SendJobConfirm(int id);
  void SendMasterIntro();
  void RecvMasterIntro(const json& msg);
  void SendJobRequest(int num);
  void RecvJobRequest(const json& msg);
  void SendHB();

  json HandleFinishMsg(json &msg);
  json HandleJobCancelMsg(json &msg);
  json HandleGetLoadMsg(json &msg);

  void BindSocket();
  void DisconnectFromMaster();
  void SendMessage(const std::string& msg, const std::string& type, size_t uid);
  void SendReq(const std::string& msg);
  void SendRep(const std::string& msg, size_t uid);
  std::string RecvReq(std::string &msg);
  void RecvRep(std::string& msg);
  void PollForMsg();
  bool ResendMsg();


};
}
}

#endif
