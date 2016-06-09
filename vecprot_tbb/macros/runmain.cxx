#include <iostream>
#include "GeantMainPropagator.h"

// Number of threads should be specified one less then desired for consistency
// There are no more separate thread manager which was added before
// Simply nthreads+1 are occupied
void run(int nthreads = 15, bool graphics = true, const char *geomfile = "../geometry/cms.root") {
  GeantMainPropagator *mainprop = GeantMainPropagator::Instance();

  // For your info. Originally min_feeder was
  // std::max<int>(50, 2*nthreads)

  // n_threads, events_total, events_buffered, tracks_average, max_per_basket
  // min_feeder, n_events_to_prioritize, threshold_to_start_DispTask
  mainprop->SetParams(nthreads, 1500, 100, 500., 10, 50, 5, 100);

  mainprop->Start(geomfile, graphics);
}

int main() {
  run();
  std::cout << "done" << std::endl;
}
