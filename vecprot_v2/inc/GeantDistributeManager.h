//===--- GeantDistributeManager.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantDistributeManager.h
 * @brief Implementation of distribute run manager
 * @author Vitalii Drohan
 */
//===----------------------------------------------------------------------===//

#ifndef GEANTV_GEANTDISTRIBUTEMANAGER_H
#define GEANTV_GEANTDISTRIBUTEMANAGER_H

#include "GeantRunManager.h"
//#include "GeantVZMQMessenger.h"
//#include "GeantEventReceiver.h"
#include "GeantConfig.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantEventReceiver;
class GeantEventDispatcher;

class GeantDistributeManger {
public:
  GeantDistributeManger(GeantRunManager *runManager, GeantConfig *config);
  void InitializeDistributedApplication(int argc, char *argv[]);
  void RunDistributedSimulation();
  ~GeantDistributeManger();

private:
  GeantRunManager *fRunMgr      = nullptr;
  GeantEventReceiver *fWorker   = nullptr;
  GeantEventDispatcher *fMaster = nullptr;

  char (*hostnames)[256];
  GeantConfig *fConfig;
};
}
}

#endif // GEANTV_GEANTDISTRIBUTEMANAGER_H
