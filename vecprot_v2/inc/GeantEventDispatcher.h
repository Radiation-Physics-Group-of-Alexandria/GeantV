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

#include "HepMC/Reader.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Search/FindParticles.h"

#include "GeantEventServer.h"
#include "GeantConfig.h"
#include "GeantRunManager.h"
#include "GeantEventReceiver.h"

#ifdef USE_HPC
#include "zmq.hpp"
#include <json.hpp>
using nlohmann::json;

#include "GeantJobPool.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

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
  zmq::context_t zmqContext;
  zmq::socket_t fSocket;


  int dispatchedEvents;

  GeantConfig *fConfig;


  GeantHPCJobPool* jobPool;
  int workerCounter;
  std::map<int,GeantHPCJob> pendingJobs;
  std::map<int,GeantHPCWorker> workers;
  std::map<std::string,int> workerId;
  std::vector<Host> fHosts;

  json HandleNewWorkerMsg(json &req);
  json HandleNewJobReq(json &req);
  json HandleJobConfirm(json &req);
  json HandleHeartbeat(json &req);

  void SendFinishMsg(GeantHPCWorker& worker);
  void SendJobCancelMsg(GeantHPCJob& job, bool retToPool = true);

  void CleanDeadWorkers();
  void FinishWorkers();

  zmq_pollitem_t zmqPollItem;
  std::map<size_t, MasterPendingMessage> pendingRequests;
  void BindSocket();
  void SendMessage(const std::string& msg, const std::string& type, size_t uid,const std::string& address);
  void SendReq(const std::string &msg, GeantHPCWorker &worker);
  void SendRep(const std::string& msg, size_t uid,const std::string& address);
  string RecvReq(const std::string &msg);
  void RecvRep(const string &msg);
  void PollForMsg();
  bool ResendMsg();

  bool IsWorkerDoingJob(GeantHPCWorker& worker, GeantHPCJob& job);
};
}
}

#endif // GEANT_EVENT_DISPATCHER_H
