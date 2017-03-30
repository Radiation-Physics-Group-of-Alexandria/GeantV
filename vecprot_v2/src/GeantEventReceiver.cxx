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

namespace Geant{
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantEventReceiver::GeantEventReceiver(GeantRunManager *runmgr)
  :fRunMgrSlave(runmgr), fStopTransfer(false){
  }

//______________________________________________________________________________
void GeantEventReceiver::Initialize(){
  fEventServer = new GeantEventServer((fRunMgrSlave->GetConfig()->fNtotal) / (fRunMgrSlave->GetConfig()->fNClients), fRunMgrSlave);
}

//______________________________________________________________________________
void GeantEventReceiver::StartSubscriber(const std::string &hostname){}

//______________________________________________________________________________
void GeantEventReceiver::StartSubscriberRouterDealer(const std::string &hostname){
  std::string host = hostname;
  std::cout << "Starting subscriber for " << hostname << std::endl;
  //TBD: will be moved to separate header GeantVZMQMessanger.h
  //std::string miy = "*";
  std::string protocol = "tcp://";
  std::string port = ":5559";
  std::string address = protocol + hostname + port;

  int rc;
  void *ctx = zmq_ctx_new();
  assert(ctx != NULL);
  void *socket = zmq_socket(ctx, ZMQ_DEALER);
  assert(socket != NULL);
  rc = zmq_connect(socket, address.c_str());
  // zmq_bind() function shall return zero if successful
  assert(rc == 0);

  zmq_msg_t d;
  zmq_msg_init(&d);
  zmq_msg_send(&d, socket, 0);
  zmq_msg_close(&d);

  zmq_msg_t start;
  zmq_msg_init(&start);
  zmq_msg_send(&start, socket, 0);
  zmq_msg_close(&start);

  zmq_pollitem_t item[] = { socket, 0, ZMQ_POLLIN, 0 };
  
  while (!fStopTransfer){
    int len = zmq_poll(item, 1, 100);
    if (len == 0){
        continue;
      }
    zmq_msg_t delimiter;
    zmq_msg_init(&delimiter);
    zmq_msg_recv(&delimiter, socket, 0);
    zmq_msg_close(&delimiter);
    int dummyeventcounter;
    zmq_msg_t dataevent;
    zmq_msg_init(&dataevent);
    zmq_msg_recv(&dataevent, socket, 0);
    std::cout << "Event was accepted :D" << std::endl;
    std::memcpy(&dummyeventcounter, zmq_msg_data(&dataevent), sizeof(dummyeventcounter));
    fRunMgrSlave->GetEventDispatcher()->SetNloadsDispatcher(fRunMgrSlave->GetEventDispatcher()->GetNloadsDispatcher() + 1);
    std::cout << "SLAVE: Number of event that was already added " << fRunMgrSlave->GetEventDispatcher()->GetNloadsDispatcher()<< std::endl;
    // To be moved to check feeder part in wmg -> EventServer::CheckNewEvent()
    fEventServer->AddEvent();
    zmq_msg_close(&dataevent);
  } 
  zmq_close(socket);
  zmq_ctx_term(ctx);
}

//______________________________________________________________________________
void GeantEventReceiver::StartPullPollTest(int number){
  bool stop = false;
  const std::string url = "tcp://127.0.0.1:5555";
  std::mutex print_lock;

  int rc;
  void *ctx = zmq_ctx_new();
  assert(ctx != NULL);

  void *socket = zmq_socket(ctx, ZMQ_PULL);
  assert(socket != NULL);
  rc = zmq_connect(socket, url.c_str());
  assert(rc == 0);
  //zmq_setsockopt(socket, ZMQ_SUBSCRIBE, "", 0);

  zmq_pollitem_t items[] = { socket, 0, ZMQ_POLLIN, 0 };
  while (!stop)
  {
    int len = zmq_poll(items, 1, 100);
    if (len == 0){
        continue;
    }
    if (items[0].revents & ZMQ_POLLIN)
    {
      std::string result;
      zmq_msg_t msg;
      zmq_msg_init(&msg);
      zmq_msg_recv(&msg, socket, 0);
      result.append(static_cast<char*>(zmq_msg_data(&msg)), zmq_msg_size(&msg));
      print_lock.lock();
      std::cout << "Worker - " << number << " : " << result << std::endl;
      print_lock.unlock();
      zmq_msg_close(&msg);
    }
  }
  zmq_close(socket);
  zmq_ctx_term(ctx);
}

//______________________________________________________________________________
void GeantEventReceiver::StartPublishSubscribeTest(const std::string &hostname){
  zmq::context_t context(1);
  zmq::socket_t subscriber(context, ZMQ_SUB);
  std::string protocol = "tcp://";
  std::string port = ":5559";
  std::string address = protocol + hostname + port;
  subscriber.connect(address);
  subscriber.setsockopt(ZMQ_SUBSCRIBE,"",0);
  zmq::message_t message;
  subscriber.recv(&message);
}

//______________________________________________________________________________
void GeantEventReceiver::StopSubscribe(){
  //zmq_close(socket);
  //zmq_ctx_term(ctx);
}

}}