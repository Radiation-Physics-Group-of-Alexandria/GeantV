//===--- GeantEventDispatcher.h - Geant-V --------------------------*- C++
//-*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEventDispatcher.h
 * @brief Implementation of Event Dispatcher in Geant-V prototype
 * @author Oksana Shadura
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
#include "GeantVZMQMessenger.h"
#include "GeantEventReceiver.h"


namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantHPCMGraph;
class GeantVWhiteBoard;
class GeantVZMQMessanger;
class EventServer;
class GeantRunManager;
class GeantEventReceiver;

class GeantEventDispatcher {

public:
  GeantEventDispatcher() {}
  /**
   * @brief [brief description]
   * @details [long description]
   *
   * @param runmgr master node dispatcher
   */
  GeantEventDispatcher(GeantRunManager *runmgr, std::string topiczmq);
  ~GeantEventDispatcher() = default;
  GeantEventDispatcher(const GeantEventDispatcher &other) = default;
  GeantEventDispatcher(GeantEventDispatcher &&other) = default;
  GeantEventDispatcher &operator=(const GeantEventDispatcher &other) = default;
  GeantEventDispatcher &operator=(GeantEventDispatcher &&other) = default;

  /**
   * @brief [brief description]
   * @details [long description]
   *
   * @param argc [description]
   * @param argv [description]
   */
  void InitializeDiscovery(int argc, char *argv[]);

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void EventFileCheck();

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void StartPublisher(const std::string &hostname);

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void ReadMessage();

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void EventMessage();

  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
  */
  void StartPublisherSubscriberTest();

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void StartPullPollTest();

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void Finalize();

  /**
   * @brief [brief description]
   * @details [long description]
   */
  void StopPublish();

  void DispLock();

  void DispUnlock();

  int GetRank() { return rank; }

  int GetNHosts() { return fNHosts; }

  GeantEventServer *GetMasterEventServer(){ return fEventServer; }

  GeantEventReceiver *GetEventReceiver(){ return fReceiver; }

  int GetNloadsDispatcher(){ return fNloadsDispatcher; } 

  void SetNloadsDispatcher( int n) { fNloadsDispatcher = n; }

private:
  GeantRunManager *fRunMgr = nullptr;       /** Run manager starting Event Dispatcher */
  GeantEventServer *fEventServer = nullptr; /** Should be a map of EventServers - rank - std::string - ptr */
  //////////////////////////////////////////////
  /// In future will be extended as a Whiteboard class
  //////////////////////////////////////////////
  std::atomic_int fNloadsDispatcher;        /** Number of events that had been passed */ 
  GeantEventReceiver *fReceiver = nullptr;   /** HPC node receivers - should be map of rank - string - ptr */
  std::atomic_flag fDispatcherLock;
  ///////////////////////////////////////////////////////////////
  /// ZMQ things
  ///////////////////////////////////////////////////////////////
  zmq::context_t *fContext;   /* ZMQ context */
  zmq::socket_t *fSocket;     /* ZMQ socket */
  volatile bool fShutdown;    /* Shutdown switch */
  std::string fEndpoint;      /* Endpoint */
  std::string fZMQTopic;      /* ZMQ topic */
  //////////////////////////////////////////////////////////////
  /// MPI
  //////////////////////////////////////////////////////////////
  std::map<int, std::string> fNodeList; /** Our node table:: 0 - master, others: ranks */
  int rank, fNHosts;                    /** MPI ranks and number of hosts in a game */
};

}}

#endif // GEANT_EVENT_DISPATCHER_H
