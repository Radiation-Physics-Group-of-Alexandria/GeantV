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
  void SendHB();

  bool GetIsTransportCompleted(){ return isTransportCompleted;}

private:
  std::mutex zmqMutex;
  zmq::context_t zmqContext;
  zmq::socket_t zmqSocket;
  zmq::socket_t zmqSocketIn;

  std::string fServHname;
  std::string fWorkerHname;
  std::string fReceiveAddr;
  GeantConfig *config;
  GeantRunManager *runManager;

  bool isTransportCompleted;
  int fReceivedEvents = 0;

  bool connected;
  int currentJob;
  int connectTries;
  bool ConnectToMaster();
  bool ConfirmJob();

  int RequestJob(int num);
  void Reconnect();
  void BindSocket();
  bool SendMsg(const std::string& req, std::string& rep);
  int wrk_id;
};
}
}

#endif
