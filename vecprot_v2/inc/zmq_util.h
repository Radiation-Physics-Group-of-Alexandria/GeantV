//===--- zmq_util.h - Geant-V --------------------------*- C++
//-*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file zmq_util.h
 * @brief Some usefull functions related to HPC communications.
 * @author Vitalii Drohan
 */
//===----------------------------------------------------------------------===//
#ifndef GEANTV_ZMQ_UTIL_H
#define GEANTV_ZMQ_UTIL_H

#include <string>
#include <chrono>

class ZmqTimer{
public:
  ZmqTimer() : lastTime(std::chrono::system_clock::now()) {};
  void Update(){ lastTime = std::chrono::system_clock::now(); }
  std::chrono::milliseconds Since(){
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - lastTime);
  }
private:
  std::chrono::time_point<std::chrono::system_clock> lastTime;
};

struct PendingMessage{
  std::string msg;
  ZmqTimer lastRetry;
  int retries;
  PendingMessage(const std::string& message) : msg{message}, lastRetry{},
                                               retries{0}{}
  PendingMessage() : msg{""}, lastRetry{}, retries{0}{}
};

struct MasterPendingMessage : public PendingMessage{
  MasterPendingMessage(const std::string& message, const std::string& address, bool toWorker) : PendingMessage(message),
                                                                            address{address}, toWorker{toWorker} {}
  MasterPendingMessage() : PendingMessage(), address{""}, toWorker{true} {}
  std::string address;
  bool toWorker;
};

inline size_t CreateMessageUID(const std::string& message){
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  return std::hash<std::string>{}(message+std::to_string(time));
}


#endif //GEANTV_ZMQ_UTIL_H
