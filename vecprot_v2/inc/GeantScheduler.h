//===--- GeantScheduler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantScheduler.h
 * @brief Implementation of the GeantV scheduler running in a single thread. Collects tracks
 * from all threads via an input queue and fills baskets corresponding to each
 * volume, which are then injected in the main work queue.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SCHEDULER
#define GEANT_SCHEDULER

#include <atomic>
#include <vector>
#include <stddef.h>
class concurrent_queue;
class GeantBasket;
class GeantBasketMgr;

#include "Geant/Typedefs.h"
#include "GeantFwd.h"

/**
 * @brief Class GeantScheduler
 * @details Dispatcher running in a single thread. Collects tracks
 * from all threads via an input queue and fills baskets corresponding to each
 * volume, which are then injected in the main work queue.
 *
 */
class GeantScheduler {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

protected:
  int fNvolumes;                                 /** Number of active volumes in the geometry */
  int fNpriority;                                /** Number of priority baskets held */
  GeantBasketMgr **fBasketMgr;                   /** Array of basket managers */
  GeantBasketMgr *fGarbageCollector;             /** Garbage collector manager */
  int *fNstvol;                                  /**[fNvolumes] Number of steps per volume */
  int *fIstvol;                                  /**[fNvolumes] Sorted index of number of steps per volume */
  int *fNvect;                                   /**[256] Number of tracks basketized in vectors of given size */
  std::atomic_long fNsteps;                      /** Total number of tracks steps */
  std::atomic_int fCrtMgr;                       /** Current basket manager being garbage collected */
  std::atomic_bool fCollecting;                  /** Flag marking colecting tracks for priority events */
  std::atomic_flag fLearning = ATOMIC_FLAG_INIT; /** Flag marking the learning phase */
  std::atomic_flag fGBCLock = ATOMIC_FLAG_INIT;  /** Flag marking that garbage collector is busy */
  std::vector<Volume_t const *> fVolumes;        /** List of logical volumes */

private:
  /**
   * @brief GeantScheduler copy constructor
   * @details Not allowed
   */
  GeantScheduler(const GeantScheduler &);

  /**
   * @brief Function operator=
   * @details Not allowed
   */
  GeantScheduler &operator=(const GeantScheduler &);

public:
  /** GeantScheduler default constructor */
  GeantScheduler();

  /** GeantScheduler destructor */
  virtual ~GeantScheduler();

  /** @brief Activate basket managers based on distribution of steps */
  void ActivateBasketManagers();

  /**
   * @brief Schedule a new track
   *
   * @param track Track to be scheduled
   */
  int AddTrack(GeantTrack &track, GeantTaskData *td);

  /**
   * @brief Re-schedule all tracks from an output basket
   * @details Transported baskets contain tracks exiting the current volume,
   * tracks killed and new tracks generated by physics along the step. These tracks
   * have to be re-scheduled for transport or bookkeeped for removal.
   *
   * @param output Transported basket
   * @param ntot Total number of tracks
   * @param nnew Number of new tracks
   * @param nkilled Number of killed tracks
   * @param td Thread data
   */
  int AddTracks(GeantTrack_v &output, int &ntot, int &nnew, int &nkilled, GeantTaskData *td);

  /** @brief Function to adjust the basket size automatically */
  void AdjustBasketSize();

  /** @brief Function to create initially baskets */
  void CreateBaskets();

  /**
   * @brief Getter for the array of basket managers
   * return Array of basket managers
   * */
  GeantBasketMgr **GetBasketManagers() const { return fBasketMgr; }

  /**
   * @brief Getter for the garbage collector
   * @return Garbage collector manager
   */
  GeantBasketMgr *GetGarbageCollector() const { return fGarbageCollector; }

  /**
   * @brief Getter for total number of steps
   *
   * @return Number of steps
   */
  long GetNsteps() const { return fNsteps.load(); }

  /**
   * @brief Getter for collecting flag
   * @return Value of fCollecting flag
   */
  bool IsCollecting() const { return fCollecting.load(); }

  /**
   * @brief Setter for collecting flag
   */
  void SetCollecting(bool flag) { fCollecting.store(flag); }

  /**
   * @brief Getter for learning flag
   * @return Value of fLearning flag
   */
  bool IsLearning() {
    bool learning = fLearning.test_and_set(std::memory_order_acquire);
    if (!learning)
      fLearning.clear(std::memory_order_release);
    return learning;
  }

  /**
   * @brief Setter for the learning flag
   */
  void SetLearning(bool flag) {
    if (flag)
      fLearning.test_and_set(std::memory_order_acquire);
    else
      fLearning.clear(std::memory_order_release);
  }

  /**
   * @brief Getter for array fNvect
   * @return Pointer to fNvect array
   */
  int *GetNvect() { return fNvect; }

  /**
   * @brief Function to return N priority baskets per volume
   * @return Number of priority baskets held
   */
  int GetNpriority() const { return fNpriority; }

  /**
   * @brief Function to return N volumes
   * @return Number of active volumes in the geometry
   */
  int GetNvolumes() const { return fNvolumes; }

  std::vector<Volume_t const *> &GetVolumes() { return fVolumes; }

  /** @brief Garbage collection function */
  int GarbageCollect(GeantTaskData *td, bool force = false);

  /** @brief Function to print size */
  void PrintSize() const;

  /** @brief Get the number of tracks that can be reused with the same thread */
  int ReusableTracks(GeantTrack_v &tracks) const;

  /** @brief Copy reusable tracks to the input of the basket */
  int CopyReusableTracks(GeantTrack_v &tracks, GeantTrack_v &input, int nmax) const;

  /** @brief Function to returns size */
  size_t Sizeof() const;
 
   void Sort( int n,int *a,int *out) {
      for (int i=0;i<n;i++) { out[i] = a[i]; }
      std::sort(out,out + n, std::greater<int>());

   } 
  
};
#endif
