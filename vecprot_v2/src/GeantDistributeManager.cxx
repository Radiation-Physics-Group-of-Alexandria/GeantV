#include "GeantDistributeManager.h"
#include "Geant/Error.h"
#include <iostream>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
GeantDistributeManger::GeantDistributeManger(GeantRunManager *runManager, GeantConfig *config)
    : fRunMgr(runManager), fConfig(config)
{
}

void GeantDistributeManger::InitializeDistributedApplication()
{
  if (fConfig->fMasterNode) {
    fMaster = new GeantEventDispatcher(fConfig);
    std::cout << "On server node" << std::endl;

  } else {
    fWorker = new GeantEventReceiver(fRunMgr, fConfig);
    std::cout << "On client node" << std::endl;
  }
  fRunMgr->SetEventReceiver(fWorker);
}

void GeantDistributeManger::RunDistributedSimulation()
{
  std::thread eventDispThread;
  if (fMaster != nullptr) {
    if (fConfig->fWorkerBootstrap != "") {
      system(fConfig->fWorkerBootstrap.c_str());
    }
    fMaster->Initialize();
    fMaster->RunReqReplyLoop();
  }
  if (fWorker != nullptr) {
    fWorker->Initialize();
    fWorker->Run();
  }
}

GeantDistributeManger::~GeantDistributeManger()
{
}
}
}
