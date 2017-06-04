#include "GeantEventReceiver.h"

#include "Geant/Error.h"

#include "GeantEvent.h"
#include "GeantEventServer.h"
#include "GeantEventDispatcher.h"
#include "GeantRunManager.h"
#include "PrimaryGenerator.h"
#include "MCTruthMgr.h"

#ifdef USE_HPC
#include "zmq.hpp"
#include "mpi.h"
#endif

#include <iostream>
#include <cstdint>
#include <cstring>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantEventReceiver::GeantEventReceiver(std::string serverHostName, GeantConfig *conf, GeantRunManager *runmgr)
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REQ), askForEventLock(0), fServHname(serverHostName), config(conf),
      runManager(runmgr), isEventReceived(false), isTransportCompleted(false)
{
}

//______________________________________________________________________________
void GeantEventReceiver::Initialize()
{
  zmqSocket.connect(std::string("tcp://") + fServHname + std::string(":5678"));
  std::cout << "HPC: Event receiver connected to: " << (std::string("tcp://") + fServHname + std::string(":5678"))
            << std::endl;
}

void GeantEventReceiver::Run()
{
  runManager->RunSimulation();
}

bool GeantEventReceiver::AskForNewEvent()
{
  if (isTransportCompleted) return false;

  isEventReceived = false;
  if (askForEventLock.fetch_add(1) == 0) { // only one thread is using ZMQ for communications
    zmq::message_t request(2);
    std::cout << "HPC: worker sent event request" << std::endl;
    zmqSocket.send(request);

    zmq::message_t reply;
    zmqSocket.recv(&reply);
    char msg[3];
    msg[2] = '\0';
    memcpy(msg, reply.data(), 2);
    std::cout << "HPC: worker received response from server: " << msg << std::endl;
    if (std::string(msg) == "OK") {
      runManager->GetEventServer()->AddEvent();
      runManager->GetEventServer()->ActivateEvents();
      isEventReceived = true;
    } else {
      isEventReceived      = false;
      isTransportCompleted = true;
    }
    askForEventLock.store(0);
  } else {
    while (askForEventLock.load() != 0) { // other threads wait
    }
  }

  return isEventReceived;
}
}
}
