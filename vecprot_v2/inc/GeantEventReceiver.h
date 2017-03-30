//===--- GeantEventReceiver.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEventReceiver.h
 * @brief Implementation of node event service (based on ZMQ) in Geant-V prototype
 * @author Oksana Shadura
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_RECEIVER
#define GEANT_RECEIVER
#ifdef USE_HPC
#include <zmq.h>
// #include "hwloc.h"
#endif

#include <string>

#include "GeantEvent.h"
#include "GeantEventServer.h"
#include "GeantEventDispatcher.h"
#include "GeantRunManager.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantEventDispatcher;

class GeantEventReceiver{
public:
    GeantEventReceiver(){}
    GeantEventReceiver(GeantRunManager *runmgr);

    ~GeantEventReceiver() = default;
    GeantEventReceiver(const GeantEventReceiver& other) = default;
    GeantEventReceiver(GeantEventReceiver&& other) = default;
    GeantEventReceiver& operator=(const GeantEventReceiver& other) = default;
    GeantEventReceiver& operator=(GeantEventReceiver&& other) = default;

    GEANT_FORCE_INLINE
  	bool IsStoppedTranfer() const { return fStopTransfer; }

  	void SetStoppedTranfer(bool switcher)  { fStopTransfer = switcher; }

  /**
   * @brief [brief description]
   * @details [long description]
   */  
    void Initialize();

  /**
   * @brief [brief description]
   * @details [long description]
   */
    void StartSubscriberRouterDealer(const std::string &hostname);

  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param event [description]
   */
    void CheckEventTransported(int event);

  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param number [description]
   */
  void StartPullPollTest(int number);

  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param hostname [description]
   */
  void StartPublishSubscribeTest(const std::string &hostname);

  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param hostname [description]
   */
  void StartSubscriber(const std::string &hostname);

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void StopSubscribe();

  void SetCounter(int counter) { fEventReceiveCounter = counter; }

  GeantEventServer *GetSlaveEventServer(){ return fEventServer; }

private:
  GeantEventDispatcher *fMasterDispatcher = nullptr;  /** Pointer to Master Event Dispatcher */ 
  bool fStopTransfer;					                        /** Bool switch to stop ZMQ transmistion */
  std::atomic_int fEventReceiveCounter;                 /** Current pointer to last event in a file */
  GeantEventServer *fEventServer = nullptr;           /** Should be a map of EventServers - rank - std::string - ptr */
  GeantRunManager *fRunMgrSlave = nullptr;            /** Run manager starting  MPI nodes */
};

}}

#endif
