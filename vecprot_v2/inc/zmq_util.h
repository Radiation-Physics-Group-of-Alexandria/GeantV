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
  MasterPendingMessage(const std::string& message, int id) : PendingMessage(message), workerId{id} {}
  MasterPendingMessage() : PendingMessage(), workerId{0} {}
  int workerId;
};


#endif //GEANTV_ZMQ_UTIL_H
