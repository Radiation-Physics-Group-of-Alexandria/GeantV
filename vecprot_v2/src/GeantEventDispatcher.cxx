#include "GeantEventDispatcher.h"
#include "GeantEventServer.h"
#include "GeantEventReceiver.h"

#include "globals.h"
#include "Geant/Error.h"

#include "GeantTrack.h"
#include "GeantEvent.h"
#include "GeantRunManager.h"
#include "LocalityManager.h"
#include "PrimaryGenerator.h"
#include "GeantTaskData.h"
#include "GeantBasket.h"
#include "MCTruthMgr.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/SimpleNavigator.h"
#include "volumes/PlacedVolume.h"
#include "management/GeoManager.h"
#else
#ifdef USE_ROOT
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#endif
#endif

#ifdef USE_HPC
#include "zmq.hpp"
#include "mpi.h"
#endif

#include <thread>
#include <functional>
#include <iostream>
#include <cstdint>
#include <cstring>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GeantEventDispatcher::GeantEventDispatcher(GeantConfig *config)
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REP), dispatchedEvents(0), fConfig(config)

{
}

void GeantEventDispatcher::Initialize()
{
  zmqSocket.bind("tcp://*:5678");
  std::cout << "HPC: Event dispatcher bounded to \"tcp://*:5678\"" << std::endl;
}

void GeantEventDispatcher::RunReqReplyLoop()
{
  while (dispatchedEvents != fConfig->fNtotal) {
    zmq::message_t request(2);
    zmqSocket.recv(&request); // Event requested

    zmq::message_t reply(2);
    memcpy(reply.data(), "OK", 2);
    zmqSocket.send(reply);
    ++dispatchedEvents;
    std::cout << "HPC: Event dispatcher event sent out to worker" << std::endl;
  }

  for (int i = 0; i < fConfig->fNClients; ++i) {
    zmq::message_t request;
    zmqSocket.recv(&request); // Event requested, but no left

    zmq::message_t reply(2);
    memcpy(reply.data(), "NO", 2);
    zmqSocket.send(reply);
    std::cout << "HPC: Event dispatcher no event msg sent to worker" << std::endl;
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
