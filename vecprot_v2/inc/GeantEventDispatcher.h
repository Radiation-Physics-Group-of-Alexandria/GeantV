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

struct MCEventSource{
  std::string fFileName;
  int fOffset;
  int fEventAmount;
  int fDispatched;
};

struct Host{
  std::string hostname;
};
class GeantHPCJobPool;


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
  GeantHPCJobPool* fJobPool;

  std::map<int,GeantHPCJob> fPendingJobs;
  int fWorkerCounter;
  std::map<int,GeantHPCWorker> fHPCWorkers;
  std::map<std::string,int> fHPCWorkerZMQIDtoUID;
  std::vector<Host> fHosts;

  std::map<size_t, MasterPendingMessage> fPendingRequests;

  json HandleNewWorkerMsg(json &req);
  json HandleNewJobReq(json &req);
  json HandleJobConfirm(json &req);
  json HandleHeartbeat(json &req);

  void SendFinishMsg(GeantHPCWorker& worker);
  void SendJobCancelMsg(GeantHPCJob& job, bool retToPool = true);

  void CleanDeadWorkers();
  void FinishWorkers();

  void BindSocket();
  void SendMessage(const std::string& msg, const std::string& type, size_t uid,const std::string& address);
  void SendReq(const std::string &msg, GeantHPCWorker &worker);
  void SendRep(const std::string& msg, size_t uid,const std::string& address);
  std::string RecvReq(const std::string &msg);
  void RecvRep(const std::string &msg);
  void PollForMsg();
  bool ResendMsg();

  bool IsWorkerDoingJob(GeantHPCWorker& worker, GeantHPCJob& job);
};
}
}

#endif // GEANT_EVENT_DISPATCHER_H
