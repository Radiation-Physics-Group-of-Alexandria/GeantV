//===--- HPCProcessManager.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file HPCProcessManager.h
 * @brief Process manager can restart managed application upon request.
 * @author Vitalii Drohan
 */
//===----------------------------------------------------------------------===//
#ifndef GEANTV_HPCPROCESSMANAGER_H
#define GEANTV_HPCPROCESSMANAGER_H

#include <zmq.hpp>
#include <json.hpp>
#include <atomic>

using nlohmann::json;

class HPCProcessManagerProxy {
public:
  HPCProcessManagerProxy(const std::string &masterHostname, int masterPort);
  void RunProxy();
  void Bind();
  void Close();

private:
  zmq::context_t context;
  zmq::socket_t frontend;
  zmq::socket_t backend;
  std::atomic<bool> closed;
  std::string masterHostname;
  int masterPort;
};

class HPCProcessManager {
public:
  HPCProcessManager(const std::string &executable, const std::string &params, const std::string &logfile);

  void RunApp();
  void PollMsg();

  bool Finished() { return isFinished; }

private:
  void KillWorker();

  zmq::context_t context;
  zmq::socket_t socket;
  zmq::pollitem_t pollitem;

  std::string executable;
  std::string params;
  std::string logfile;
  pid_t geantPid     = -1;
  bool isRunning     = false;
  bool isFinished    = false;
  bool restartOnFail = true;

  void ForceFinish();

  std::string RecvReq(const std::string &msg);
  json HandleStatusMsg();
  json HandleAbortMsg();
  json HandleRestartMsg();
};

#endif // GEANTV_HPCPROCESSMANAGER_H
