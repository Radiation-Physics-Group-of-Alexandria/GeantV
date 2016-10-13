#include "EventLoopTask.h"

#include "GeantRunManager.h"
#include "GeantApplicationTBB.h"
#include "Geant/Error.h"

#include "tbb/task_scheduler_init.h"

using namespace Geant;

//______________________________________________________________________________
EventLoopTask::EventLoopTask(GeantRunManager *runmgr)
{
// Constructor
  fRunMgr = runmgr;
  fApplication = dynamic_cast<GeantApplicationTBB*>(runmgr->GetUserApplication());
  if (!fApplication) {
    Fatal("InitializerTask", "No TBB user application connected to the run manager");
  }
}

//______________________________________________________________________________
tbb::task* EventLoopTask::execute ()
{
  if (!fRunMgr->IsInitialized()) fRunMgr->Initialize();
  return fApplication->SpawnUserEventFeeder(fRunMgr->GetEventServer());
}
