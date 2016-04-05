#ifndef GEANT_TBB_FEEDER_TASK
#define GEANT_TBB_FEEDER_TASK

#include "ThreadData.h"
#include "TransportTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "GeantEvent.h"
#include "GeantTaskData.h"
#include "tbb/task.h"

class FeederTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;

public:
  FeederTask (Geant::GeantTaskData *td);
  ~FeederTask ();

  tbb::task* execute ();

};

#endif //GEANT_TBB_FEEDER_TASK