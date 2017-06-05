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
#endif

#include <string>

#include "GeantRunManager.h"
#include "GeantConfig.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantEventDispatcher;

class GeantEventReceiver {
public:
  GeantEventReceiver(std::string serverHostName, GeantConfig *conf, GeantRunManager *runmgr);

  void Initialize();
  void Run();
  int AskForNewEvent(int num);

  bool GetIsTransportCompleted(){ return isTransportCompleted;}

private:
  zmq::context_t zmqContext;
  zmq::socket_t zmqSocket;

  std::atomic_int askForEventLock;

  std::string fServHname;
  GeantConfig *config;
  GeantRunManager *runManager;

  bool isTransportCompleted;
  int fReceivedEvents = 0;
};
}
}

#endif
