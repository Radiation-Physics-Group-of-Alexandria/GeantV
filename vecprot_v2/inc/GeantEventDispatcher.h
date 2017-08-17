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
  zmq::socket_t zmqSocket;

  int dispatchedEvents;

  GeantConfig *fConfig;


  GeantHPCJobPool* jobPool;
  int workerCounter;
  std::vector<GeantHPCJob> pendingJobs;
  std::vector<GeantHPCWorker> workers;
  std::vector<Host> fHosts;

  json NewWorker(json& req);
  json JobReq(json& req);
  json JobConfirm(json& req);
  json HeartBeat(json& req);
  void CleanDeadWorkers();
};
}
}

#endif // GEANT_EVENT_DISPATCHER_H
