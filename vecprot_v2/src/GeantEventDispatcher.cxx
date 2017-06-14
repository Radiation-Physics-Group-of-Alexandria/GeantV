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
#include <fstream>
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

  std::ifstream inputFile(fConfig->fEventListFilename) ;
  std::string line;
  MCEventSource source;
  while(inputFile >> line){ //line format: pathToFile:eventOffset:eventAmount
    auto col1 = line.find(':',0);
    auto col2 = line.find(':',col1+1);

    source.fFileName = line.substr(0, col1);
    source.fOffset = std::stoi(line.substr(col1+1, col2-col1-1));
    source.fEventAmount = std::stoi(line.substr(col2+1));
    source.fDispatched = 0;

    fEventSources.push_back(source);
  }

  currentSource = 0;

}

void GeantEventDispatcher::RunReqReplyLoop()
{
  while(currentSource < fEventSources.size()){
    auto& currSource = fEventSources[currentSource];

    zmq::message_t request(2);
    zmqSocket.recv(&request); // Event requested
    char requestMsg[10] = "0";
    memcpy(requestMsg,request.data(),10);
    int requestedAmount = 0;
    sscanf(requestMsg,"%d",&requestedAmount);

    int eventGivenOut = std::min(currSource.fEventAmount - currSource.fDispatched, requestedAmount);


    zmq::message_t reply(4096+1+10+1+10); //Max path name + : + offset + : + amount
    std::string replyString = currSource.fFileName + ":" + std::to_string(currSource.fOffset + currSource.fDispatched) +
        ":" + std::to_string(eventGivenOut);
    memcpy(reply.data(), replyString.c_str(), replyString.size()+1); //copy null byte
    zmqSocket.send(reply);
    ++dispatchedEvents;

    currSource.fDispatched += eventGivenOut;
    if(currSource.fDispatched == currSource.fEventAmount) currentSource++;
  }

  for (int i = 0; i < fConfig->fNClients; ++i) {
    zmq::message_t request;
    zmqSocket.recv(&request); // Event requested, but no left

    zmq::message_t reply(4096+1+10+1+10);
    memcpy(reply.data(), "NO", 3);
    zmqSocket.send(reply);
    std::cout << "HPC: Event dispatcher no event msg sent to worker" << std::endl;
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
