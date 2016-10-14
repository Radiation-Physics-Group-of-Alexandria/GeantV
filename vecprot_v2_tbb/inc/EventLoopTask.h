#ifndef GEANT_TBB_EVENT_LOOP_TASK
#define GEANT_TBB_EVENT_LOOP_TASK

/**
 * @file EventLoopTask.h
 * @brief Task spawned by the user application in its own event loop.
 * @author Andrei Gheata
 */

#include "tbb/task.h"
#include "GeantFwd.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class GeantRunManager;
  class GeantApplicationTBB;
}}


class EventLoopTask : public tbb::task
{
private:
  Geant::GeantRunManager     *fRunMgr      = nullptr; /** Run manager */
  Geant::GeantApplicationTBB *fApplication = nullptr; /** TBB user application */
  int fNevents = 0;                                   /** Number of events to be transported */

public:
  EventLoopTask(Geant::GeantRunManager *runmgr, int nevents);
  ~EventLoopTask() {}

  tbb::task* execute();
};

#endif //GEANT_TBB_EVENT_LOOP_TASK
