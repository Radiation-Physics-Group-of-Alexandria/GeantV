//===--- GeantVApplication.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantApplicationTBB.h
 * @brief Base class for TBB based user application
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_APPLICATION_TBB
#define GEANT_APPLICATION_TBB

#include "GeantVApplication.h"
#include "tbb/task.h"
#include "EventLoopTask.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantEventServer;

/** @brief GeantVApplication class */
class GeantApplicationTBB : GeantVApplication {
public:

  /** @brief GeantApplicationTBB constructor */	
  GeantApplicationTBB(GeantRunManager *runmgr) : GeantVApplication(runmgr) {}
  
  /** @brief GeantApplicationTBB destructor */
  virtual ~GeantApplicationTBB() {}

  /** @brief Function of initialization */
  virtual bool Initialize() = 0;

  /**
   * @brief User callback function for scoring
   * 
   * @param tid  Thread ID
   * @param npart Number of tracks
   * @param tracks Set of tracks
   */
  virtual void StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) = 0;

  /**
   * @brief Function of digitization
   * 
   * @param event Event for digitization
   */
  virtual void Digitize(GeantEvent *event) = 0;

  /** @brief User FinishRun function */
  virtual void FinishRun() = 0;
  
  /** @brief User actions in terms of TBB tasks */
  virtual tbb::task *SpawnUserEndRunTask() = 0;

  /** @brief Spawn the event loop task */
  tbb::task *SpawnEventLoopTask(int nevents) {
    tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
    EventLoopTask & evloopTask = *new(cont.allocate_child()) EventLoopTask( fRunMgr, nevents );
    return & evloopTask;
  }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
