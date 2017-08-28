//===--- GeantEventDispatcher.h - Geant-V --------------------------*- C++
//-*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEventDispatcher.h
 * @brief Implementation of Event Dispatcher in Geant-V prototype
 * @author Oksana Shadura, Vitalii Drohan
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_EVENT_DISPATCHER_H
#define GEANT_EVENT_DISPATCHER_H

#include <map>
#include <deque>
#include <zmq.hpp>
#include <json.hpp>

#include "GeantEventServer.h"
#include "GeantConfig.h"
#include "GeantRunManager.h"
#include "GeantEventReceiver.h"

#include "GeantJobPool.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using nlohmann::json;

class GeantHPCJobPool;

struct Host {
  Host(const std::string &zmqAddress) : fZMQAddress(zmqAddress){};
  Host() : fZMQAddress(""){};
  std::string fZMQAddress;
  int fReconnects                                = 0;
  int fEventsGiven                               = 0;
  int fEventsConfirmed                           = 0;
  std::chrono::milliseconds fAverageTimeForEvent = std::chrono::milliseconds(0);
  double fLastMinLoadAverage                     = 0.0;
};

class GeantEventDispatcher {
public:
  GeantEventDispatcher(GeantConfig *config);

  void Initialize();
  void RunReqReplyLoop();

private:
  zmq::context_t fZMQContext;
  zmq::socket_t fSocket;
  zmq_pollitem_t fZMQSocketPollItem;

  GeantConfig *fConfig;
  GeantHPCJobPool *fJobPool;

  std::map<int, GeantHPCJob> fPendingJobs;
  int fWorkerCounter;
  std::map<int, GeantHPCWorker> fHPCWorkers;
  std::map<std::string, int> fHPCWorkerZMQIDtoUID;
  std::map<std::string, Host> fHPCHosts;

  std::map<size_t, MasterPendingMessage> fPendingRequests;

  json HandleNewWorkerMsg(const std::string &address);
  json HandleNewJobReq(json &req);
  bool FindDuplicateJob(GeantHPCWorker &worker, GeantHPCJob &job, int n);
  json HandleJobConfirm(json &req);
  json HandleHeartbeat(json &req);

  void SendFinishMsg(GeantHPCWorker &worker);
  void SendJobCancelMsg(GeantHPCJob &job, bool retToPool = true);
  void SendGetLoadMsg(GeantHPCWorker &worker);
  void RecvGetLoadMsg(const json &msg, GeantHPCWorker &worker);
  void SendGetStatusMsg(const std::string &address);
  void RecvGetStatusMsg(const json &msg, const std::string &address);
  void SendRestartMsg(const std::string &address);
  void SendAbortMsg(const std::string &address);

  void UpdateWorkerStats();
  void CleanDeadWorkers();
  void CleanDeadJobs();
  void FinishWorkers();

  void BindSocket();
  void SendMessage(const std::string &msg, const std::string &type, size_t uid, const std::string &address,
                   const std::string &localAddress);
  void SendReqWorker(const std::string &msg, GeantHPCWorker &worker);
  void SendReqProcMgr(const std::string &msg, const std::string &address);
  void SendRep(const std::string &msg, size_t uid, const std::string &address, const std::string &localAddress);
  std::string RecvReqWorker(const std::string &msg, const std::string &address);
  void RecvRepWorker(const std::string &msg, GeantHPCWorker &worker);
  void RecvRepProcMgr(const std::string &msg, const std::string &address);
  void PollForMsg();
  bool ResendMsg();
  void CleanWorkerState(GeantHPCWorker &worker);

  bool IsWorkerDoingJob(GeantHPCWorker &worker, GeantHPCJob &job);
};
}
}

#endif // GEANT_EVENT_DISPATCHER_H
