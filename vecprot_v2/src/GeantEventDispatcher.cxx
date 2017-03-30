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

//______________________________________________________________________________
GeantEventDispatcher::GeantEventDispatcher(GeantRunManager *runmgr, std::string str)
    : fRunMgr(runmgr), fZMQTopic(str), fContext(), fNodeList(), fNloadsDispatcher(0),
     fReceiver(), fNHosts(), rank(){
      fDispatcherLock.clear();
     }

//____________________________________________________________________________
  void GeantEventDispatcher::DispLock() {
    while (fDispatcherLock.test_and_set(std::memory_order_acquire))
      std::cout << "DispServer is locked.."<<std::endl;
  }

//____________________________________________________________________________
  void GeantEventDispatcher::DispUnlock() {
    fDispatcherLock.clear(std::memory_order_release); 
  }

//______________________________________________________________________________
void GeantEventDispatcher::EventFileCheck(){}

//______________________________________________________________________________
void GeantEventDispatcher::InitializeDiscovery(int argc, char *argv[]){
  // Rank in MPI set of hosts
  int fNSubsribers, fNPublishers, proc;

  std::cout << "Starting node discovery.." << std::endl;

  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    assert("MPI_Init failed");
    std::abort();
  }
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  MPI_Comm_size(MPI_COMM_WORLD, &fNHosts);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Number of hosts in the game: " << fNHosts << std::endl;
    assert(fNHosts >= 2);
  }

  char procname[MPI_MAX_PROCESSOR_NAME];
  if (MPI_SUCCESS != MPI_Get_processor_name(procname, &proc)) {
    strcpy(procname, "unknown");
  }
  for (int i = 0; i < fNHosts; ++i) {
    fNodeList.emplace(std::make_pair(i, std::string(procname)));
  }
  std::cout << "Master processor: " << procname << std::endl;
  // Creating map of nodes
  for (const auto &p : fNodeList) {
    std::cout << "Node list:" << std::endl;
    std::cout << p.first << " => " << p.second << "\n";
  }

  if (fNHosts == 1) {
    // Until we dont have other node - we are using separate C++11 thread for Publisher
    std::cout << "MPI Rank = 0, but we don't have slaves => single node transport (ZMQ)" << std::endl;
    fReceiver = new GeantEventReceiver(fRunMgr);
    fReceiver->Initialize();
    fEventServer = new GeantEventServer((fRunMgr->GetConfig()->fNtotal) / (fRunMgr->GetConfig()->fNClients), fRunMgr);
    std::thread publisher([&] { this->StartPublisher(fNodeList[0]); });
    std::thread subscriber([&] { fReceiver->StartSubscriberRouterDealer(fNodeList[0]); });
    std::this_thread::sleep_for(std::chrono::seconds(10));
    fReceiver->SetStoppedTranfer(true);
    publisher.join();
    subscriber.join();
    for (int i = 0; i < (fRunMgr->GetConfig()->fNtotal) / (fRunMgr->GetConfig()->fNClients); ++i)
      fEventServer->AddEvent();
  }
  // Scheduling jobs
  if (rank == 0) {
    std::cout << "MPI Rank = 0" << std::endl;
    // Starting publish information from master
    std::thread publisher([&] { this->StartPublisher(fNodeList[rank]); });
    publisher.join();
    ////////////////////////////////////////////////////////////////////////
    /// Starting Event Server for master node (can't separate Propagator and WorkloadManager)
    ////////////////////////////////////////////////////////////////////////
    fEventServer = new GeantEventServer(((fRunMgr->GetConfig()->fNtotal) / (fRunMgr->GetConfig()->fNClients)), fRunMgr);
    // To be moved to check feeder part in wmg -> EventServer::CheckNewEvent()
    for (int i = 0; i < (fRunMgr->GetConfig()->fNtotal) / (fRunMgr->GetConfig()->fNClients); ++i){
      fEventServer->AddEvent();
      fNloadsDispatcher.fetch_add(1);
      std::cout << "MASTER: Number of event that was already added " << this->GetNloadsDispatcher() << std::endl;
    }
  } else {
    std::cout << "MPI Rank =! 0" << std::endl;
    // Starting subscribers on nodes (clients)
    // Change idea of fRunMgr->GetConfig()->fNtotal -> we need to tune a percentage of events spread on nodes
    fReceiver = new GeantEventReceiver(fRunMgr);
    fReceiver->Initialize();
    std::thread subscriber([&] {fReceiver->StartSubscriberRouterDealer(fNodeList[rank]);});
    subscriber.join();
    if(fNloadsDispatcher.load() == fRunMgr->GetConfig()->fNtotal)
      fReceiver->SetStoppedTranfer(true);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

//______________________________________________________________________________
void GeantEventDispatcher::ReadMessage(){
}

//______________________________________________________________________________
void GeantEventDispatcher::EventMessage(){
}

//______________________________________________________________________________
void GeantEventDispatcher::StartPublisher(const std::string &hostname){
  // The context is the container for all sockets in a single process,
  // and acts as the transport for inproc sockets, which are the fastest way to connect threads in one process.
  // TBD: will be moved to separate header GeantVZMQMessanger.h

  std::string miy = "*";
  std::string protocol = "tcp://";
  std::string port     = ":5559";
  std::string address  = protocol + miy + port;
  std::cout << "Starting publisher on EventDispatcher " << address << std::endl;

  int rc;
  void *ctx = zmq_ctx_new();
  assert(ctx != NULL);
  // A socket of type ZMQ_ROUTER is an advanced socket type
  // used for extending request/reply sockets.
  // When receiving messages a ZMQ_ROUTER socket shall prepend a message part containing the identity of
  // the originating peer to the message before passing it to the application.
  void *socket = zmq_socket(ctx, ZMQ_ROUTER);
  assert(socket != NULL);
  rc = zmq_bind(socket, address.c_str());
  // zmq_bind() function shall return zero if successful
  assert(rc == 0);

  zmq_msg_t id;
  zmq_msg_init(&id);
  zmq_msg_recv(&id, socket, 0);
  std::string id_copy;
  id_copy.append((char *)zmq_msg_data(&id), zmq_msg_size(&id));
  zmq_msg_close(&id);

  //////////////////////////////
  zmq_msg_t d;
  zmq_msg_init(&d);
  zmq_msg_recv(&d, socket, 0);
  zmq_msg_close(&d);

  zmq_msg_t tmp;
  zmq_msg_init(&tmp);
  zmq_msg_recv(&tmp, socket, 0);
  zmq_msg_close(&tmp);

  for (int i = 0; i < (fRunMgr->GetConfig()->fNtotal)/(fRunMgr->GetConfig()->fNClients); ++i) {
    //////////////////////////////
    zmq_msg_t id;
    zmq_msg_init_size(&id, id_copy.size());
    memcpy(zmq_msg_data(&id), id_copy.c_str(), id_copy.size());
    zmq_msg_send(&id, socket, ZMQ_SNDMORE);
    zmq_msg_close(&id);
    /////////////////////////////
    //ØMQ "multi-part message": allows to concatenate multiple messages into a single message
    zmq_msg_t d;
    zmq_msg_init(&d);
    zmq_msg_send(&d, socket, ZMQ_SNDMORE);
    zmq_msg_close(&d);
    /////////////////////////////
    fNloadsDispatcher.fetch_add(1);
    zmq_msg_t dataevent;
    zmq_msg_init_size(&dataevent, sizeof(fNloadsDispatcher));
    std::memcpy(&fNloadsDispatcher,zmq_msg_data(&dataevent), sizeof(int));
    zmq_msg_send(&dataevent, socket, 0);
    zmq_msg_close(&dataevent);
  }
  zmq_close(socket);
  zmq_ctx_term(ctx);
}

//______________________________________________________________________________
void GeantEventDispatcher::StartPublisherSubscriberTest(){
  const std::string url = "tcp://127.0.0.1:5559";
  int rc;
  void *ctx = zmq_ctx_new();
  assert(ctx != NULL);
  void *socket = zmq_socket(ctx, ZMQ_PUB);
  assert(socket != NULL);
  rc = zmq_bind(socket, url.c_str());
  assert(rc == 0);
  for(int i = 0; i < (fRunMgr->GetConfig()->fNtotal)/(fRunMgr->GetConfig()->fNClients); ++i){
    fNloadsDispatcher.fetch_add(1);
    zmq_msg_t dataevent;
    zmq_msg_init_size(&dataevent, sizeof(fNloadsDispatcher));
    std::memcpy(&fNloadsDispatcher, zmq_msg_data(&dataevent), sizeof(int));
    zmq_msg_send(&dataevent, socket, 0);
    zmq_msg_close(&dataevent);
  }
  zmq_close(socket);
  zmq_ctx_term(ctx);
}

//______________________________________________________________________________
void GeantEventDispatcher::StartPullPollTest(){
	const std::string url = "tcp://127.0.0.1:5559";
	int rc;
	void *ctx = zmq_ctx_new();
	assert(ctx != NULL);
	void *socket = zmq_socket(ctx, ZMQ_PUSH);
	assert(socket != NULL);
	rc = zmq_bind(socket, url.c_str());
	assert(rc == 0);
	//it's important to wait for all connections established
	sleep(1000);
	for(int i = 1; i < (fRunMgr->GetConfig()->fNtotal)/(fRunMgr->GetConfig()->fNClients); ++i){
    fNloadsDispatcher.fetch_add(1);
		zmq_msg_t dataevent;
    zmq_msg_init_size(&dataevent, sizeof(fNloadsDispatcher));
    std::memcpy(&fNloadsDispatcher, zmq_msg_data(&dataevent), sizeof(int));
    zmq_msg_send(&dataevent, socket, 0);
    zmq_msg_close(&dataevent);
	}
	zmq_close(socket);
	zmq_ctx_term(ctx);
}

//______________________________________________________________________________
void GeantEventDispatcher::StopPublish(){}

//______________________________________________________________________________
void GeantEventDispatcher::Finalize(){
  MPI_Finalize();
}

} // GEANT_IMPL_NAMESPACE
} // Geant
