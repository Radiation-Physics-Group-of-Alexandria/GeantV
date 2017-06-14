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
#include <HepMCGeneratorMultFiles.h>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantEventReceiver::GeantEventReceiver(std::string serverHostName, GeantConfig *conf, GeantRunManager *runmgr)
    : zmqContext(1), zmqSocket(zmqContext, ZMQ_REQ), fServHname(serverHostName), config(conf),
      runManager(runmgr), isTransportCompleted(false)
{
}

//______________________________________________________________________________
void GeantEventReceiver::Initialize()
{
  std::ifstream inputFile(config->fEventListFilename) ;
  std::string line;
  int nEvents = 0;
  while(inputFile >> line){ //line format: pathToFile:eventOffset:eventAmount
    auto col1 = line.find(':',0);
    auto col2 = line.find(':',col1+1);
    nEvents += std::stoi(line.substr(col2+1));
  }
  config->fNtotal = nEvents;


  zmqSocket.connect(std::string("tcp://") + fServHname + std::string(":5678"));
  std::cout << "HPC: Event receiver connected to: " << (std::string("tcp://") + fServHname + std::string(":5678"))
            << std::endl;
  fReceivedEvents = 0;
}

void GeantEventReceiver::Run()
{
  runManager->RunSimulation();
}

int GeantEventReceiver::AskForNewEvent(int num)
{
  if (isTransportCompleted) return 0;
      fReceivedEvents = 0;

      zmq::message_t request(10);
      std::string strReq = std::to_string(num);
      memcpy(request.data(),strReq.c_str(),strReq.size());
      std::cout << "HPC: worker sent event request" << std::endl;
      zmqSocket.send(request);

      zmq::message_t reply;
      zmqSocket.recv(&reply);
      char msg[4096+1+10+1+10+1];
      memcpy(msg, reply.data(), 4096+1+10+1+10);
      std::string strMsg = msg;

      std::cout << "HPC: worker received response from server: " << msg << std::endl;
      if (strcmp(msg,"NO") == 0) {
        isTransportCompleted = true;
      } else {
        auto multFileGenerator = (HepMCGeneratorMultFiles*)runManager->GetPrimaryGenerator();
        auto col1 = strMsg.find(':',0);
        auto col2 = strMsg.find(':',col1+1);

        std::string eventFileName = strMsg.substr(0, col1);
        int offset = std::stoi(strMsg.substr(col1+1, col2-col1-1));
        int amount = std::stoi(strMsg.substr(col2+1));
        multFileGenerator->SetEventSource(eventFileName, offset);
        for (int j = 0; j < amount; ++j) {
          runManager->GetEventServer()->AddEvent();
          ++fReceivedEvents;
        }
        runManager->GetEventServer()->ActivateEvents();
      }
  return fReceivedEvents;
}
}
}
