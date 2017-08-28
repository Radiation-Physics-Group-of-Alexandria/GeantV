#include "HPCProcessManager.h"
#include <string>
#include <iostream>
#include <atomic>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <sstream>
#include <thread>
#include <getopt.h>
#include <err.h>

std::string executable;
std::string params;
std::string logfile = "";
std::string masterHostname;
int masterPort;

static struct option options[] = {
    {"executable", required_argument, 0, 'e'},  {"params", required_argument, 0, 'p'},
    {"log file", required_argument, 0, 'l'},    {"master host name", required_argument, 0, 'H'},
    {"master port", required_argument, 0, 'P'}, {0, 0, 0, 0}};
void help()
{
  printf("\nUsage: process manager [OPTIONS] \n\n");

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

int main(int argc, char **argv)
{
  if (argc == 1) {
    help();
    return 0;
  }
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "e:p:l:H:P:", options, &optidx);

    if (c == -1) break;

    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'e':
      executable = optarg;
      break;
    case 'p':
      params = optarg;
      break;
    case 'l':
      logfile = optarg;
      break;
    case 'H':
      masterHostname = optarg;
      break;
    case 'P':
      masterPort = (int)strtol(optarg, NULL, 10);
      break;
    default:
      errx(1, "unknown option %c", c);
    }
  }
  HPCProcessManager procMgr(executable, params, logfile);
  HPCProcessManagerProxy procMgrProxy(masterHostname, masterPort);

  procMgrProxy.Bind();
  std::thread proxyThread([&]() { procMgrProxy.RunProxy(); });

  procMgr.RunApp();

  procMgrProxy.Close();
  proxyThread.join();
  return 0;
}

HPCProcessManagerProxy::HPCProcessManagerProxy(const std::string &masterHostname, int masterPort)
    : context(1), frontend(context, ZMQ_DEALER), backend(context, ZMQ_ROUTER), closed(false),
      masterHostname(masterHostname), masterPort(masterPort)
{
}

void HPCProcessManagerProxy::RunProxy()
{
  zmq::pollitem_t pollItems[] = {{frontend, 0, ZMQ_POLLIN, 0}, {backend, 0, ZMQ_POLLIN, 0}};
  while (!closed) {
    zmq::poll(pollItems, 2, std::chrono::milliseconds(10));
    if (pollItems[0].revents & ZMQ_POLLIN) {
      while (1) {
        zmq::message_t inMsg;
        frontend.recv(&inMsg);
        if (inMsg.more()) {
          backend.send(inMsg, ZMQ_SNDMORE);
        } else {
          backend.send(inMsg);
          break;
        }
      }
    }
    if (pollItems[1].revents & ZMQ_POLLIN) {
      while (1) {
        zmq::message_t inMsg;
        backend.recv(&inMsg);
        if (inMsg.more()) {
          frontend.send(inMsg, ZMQ_SNDMORE);
        } else {
          frontend.send(inMsg);
          break;
        }
      }
    }
  }

  frontend.setsockopt(ZMQ_LINGER, 0);
  backend.setsockopt(ZMQ_LINGER, 0);
}

void HPCProcessManagerProxy::Close()
{
  closed = true;
}

void HPCProcessManagerProxy::Bind()
{
  char hostname[256];
  gethostname(hostname, 256);
  auto zmqId = hostname + std::string("+") + std::to_string(getpid());
  frontend.setsockopt(ZMQ_IDENTITY, zmqId.c_str(), zmqId.size());
  frontend.connect("tcp://" + masterHostname + ":" + std::to_string(masterPort));
  backend.bind("ipc:///tmp/geantV_back.pipe");
}

HPCProcessManager::HPCProcessManager(const std::string &executable, const std::string &params,
                                     const std::string &logfile)
    : context(zmq::context_t(1)), socket(zmq::socket_t(context, ZMQ_DEALER)), executable(executable), params(params),
      logfile(logfile)
{
  socket.setsockopt(ZMQ_IDENTITY, "MGR", 3);
  auto masterZMQAddress = std::string("ipc://") + "/tmp/geantV_back.pipe";
  socket.connect(masterZMQAddress);
  pollitem = {socket, 0, ZMQ_POLLIN, 0};
}

void HPCProcessManager::PollMsg()
{
  zmq::poll(&pollitem, 1, std::chrono::milliseconds(10));
  if (pollitem.revents & ZMQ_POLLIN) {
    zmq::message_t type;
    zmq::message_t muid;
    zmq::message_t message;
    socket.recv(&type);
    socket.recv(&muid);
    socket.recv(&message);
    std::string messageType((char *)type.data(), (char *)type.data() + type.size());
    size_t messageUid;
    memcpy(&messageUid, muid.data(), muid.size());
    std::string messageContent((char *)message.data(), (char *)message.data() + message.size());
    if (messageType == "REQ") {
      std::string replyContent = RecvReq(messageContent);
      zmq::message_t repType(3);
      memcpy(repType.data(), "REP", 3);
      zmq::message_t repMuid(sizeof(messageUid));
      memcpy(repMuid.data(), &messageUid, sizeof(messageUid));
      zmq::message_t repMessage(replyContent.size());
      memcpy(repMessage.data(), replyContent.c_str(), replyContent.size());
      socket.send(repType, ZMQ_SNDMORE);
      socket.send(repMuid, ZMQ_SNDMORE);
      socket.send(repMessage);
    }
  }
}

std::vector<std::string> splitParamString(const std::string &prm)
{
  std::istringstream iss(prm);
  std::vector<std::string> tokens{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
  return tokens;
}

void HPCProcessManager::RunApp()
{
  const int kMaxRetries = 3;
  int retries = 0;
  while (restartOnFail && retries < kMaxRetries) {
    retries += 1;
    int status;
    geantPid = fork();
    if (geantPid == -1) {
      std::cerr << "Fork failed" << '\n';
      break;
    }

    if (geantPid == 0) { // child
      if (logfile != "") {
        int fd = open(logfile.c_str(), O_APPEND | O_CREAT | O_WRONLY, 0666);
        dup2(fd, 1);
        dup2(fd, 2);
      }

      auto splitParams  = splitParamString(params);
      const char **argv = new const char *[splitParams.size() + 2];
      argv[0]           = executable.c_str();
      for (int i = 1; i <= splitParams.size(); ++i) {
        argv[i] = splitParams[i - 1].c_str();
      }
      argv[splitParams.size() + 1] = NULL;
      execv(executable.c_str(), (char **)argv);
    }

    isRunning = true;
    std::cout << "Child pid: " << geantPid << '\n';
    pid_t c_pid = 0;
    while (c_pid == 0) {
      c_pid = waitpid(geantPid, &status, WNOHANG);
      PollMsg();
    }
    isRunning = false;
    if (WIFEXITED(status)) {
      break;
    }
  }
  std::cout << "FINISHED" << std::endl;
  isFinished = true;
}

void HPCProcessManager::KillWorker()
{
  if (geantPid == -1) return;
  kill(geantPid, SIGTERM);
  sleep(2);
  kill(geantPid, SIGKILL);
}

void HPCProcessManager::ForceFinish()
{
  restartOnFail = false;
  KillWorker();
}

std::string HPCProcessManager::RecvReq(const std::string &msg)
{
  json req     = json::parse(msg);
  auto msgType = req["type"].get<std::string>();
  if (msgType == "status_req") {
    return HandleStatusMsg().dump();
  } else if (msgType == "abort_req") {
    return HandleAbortMsg().dump();
  } else if (msgType == "restart_req") {
    return HandleRestartMsg().dump();
  }
  return "{}";
}

json HPCProcessManager::HandleStatusMsg()
{
  json rep;
  rep["type"]    = "status_rep";
  rep["running"] = isRunning;
  return rep;
}

json HPCProcessManager::HandleAbortMsg()
{
  ForceFinish();
  json rep;
  rep["type"] = "abort_rep";
  return rep;
}

json HPCProcessManager::HandleRestartMsg()
{
  KillWorker();
  json rep;
  rep["type"] = "restart_rep";
  return rep;
}
