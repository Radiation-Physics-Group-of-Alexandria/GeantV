#include "GeantDistributeManager.h"
#include "Geant/Error.h"
#include <iostream>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
GeantDistributeManger::GeantDistributeManger(GeantRunManager *runManager, GeantConfig *config)
    : fRunMgr(runManager), fConfig(config)
{
}

void GeantDistributeManger::InitializeDistributedApplication(int argc, char *argv[])
{
  char hostname[256];
  gethostname(hostname, sizeof(hostname));

  if (fConfig->fMasterNode) {
    fMaster = new GeantEventDispatcher(fConfig);
    fWorker = new GeantEventReceiver("localhost", fConfig, fRunMgr);
    std::cout << "On server node" << std::endl;

  } else {
    fWorker = new GeantEventReceiver(fConfig->fMasterHostname, fConfig, fRunMgr);
    std::cout << "On client node" << std::endl;
  }
  fRunMgr->SetEventReceiver(fWorker);
}

void GeantDistributeManger::RunDistributedSimulation()
{
  std::thread eventDispThread;
  if (fMaster != nullptr) {
    eventDispThread = std::thread([&] {
      fMaster->Initialize();
      fMaster->RunReqReplyLoop();
    });
  }
  if (fWorker != nullptr) {
    fWorker->Initialize();
    fWorker->Run();
  }
  if (fMaster != nullptr) {
    eventDispThread.join();
  }
}

GeantDistributeManger::~GeantDistributeManger()
{
}
}
}
