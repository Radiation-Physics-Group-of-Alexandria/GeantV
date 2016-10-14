#include "EventLoopTask.h"

#include "GeantRunManager.h"
#include "GeantApplicationTBB.h"
#include "Geant/Error.h"

#include "tbb/task_scheduler_init.h"

using namespace Geant;

//______________________________________________________________________________
EventLoopTask::EventLoopTask(GeantRunManager *runmgr, int nevents)
              :fRunMgr(runmgr), fNevents(nevents)
{
// Initialize an event loop task to transport nevents
  fApplication = dynamic_cast<GeantApplicationTBB*>(runmgr->GetUserApplication());
  if (!fApplication) {
    Fatal("InitializerTask", "No TBB user application connected to the run manager");
  }
}

//______________________________________________________________________________
tbb::task* EventLoopTask::execute ()
{
  // Steer the event loop 
  if (!fRunMgr->IsInitialized()) fRunMgr->Initialize();
  // The loop below will trigger calls to the primary generator. The user
  // code can assemble the events there using the NextEvent interface offered
  // by the primary generator.
  for (auto i=0; i<fNevents; ++i)
    fRunMgr->GetEventServer()->AddEvent();
  
  // Now execute the full GeantV loop.  
  fRunMgr->RunSimulation();
  
  // Spawn the user end of run task to complete the loop
  return fApplication->SpawnUserEndRunTask();
}
