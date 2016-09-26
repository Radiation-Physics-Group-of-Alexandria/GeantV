//===--- GeantPropagator.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantPropagator.h
 * @brief Implementation of propogator in Geant-V prototype.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PROPAGATOR
#define GEANT_PROPAGATOR

#ifndef GEANT_TRACK
#include "GeantTrackVec.h"
#endif

#include <vector>
#include <atomic>
#include <mutex>

#include "Geant/Typedefs.h"
#ifdef USE_ROOT
class TTree;
class TFile;
class TStopwatch;
#else
#include "base/Stopwatch.h"
#endif
#include "base/BitSet.h"
using veccore::BitSet;
class PhysicsProcess;
class GeantEvent;
class GeantBasket;
class GeantBasketMgr;
class WorkloadManager;
class GeantVApplication;
class GeantVTaskMgr;
class PrimaryGenerator;
class MCTruthMgr;
class TaskBroker;
#if USE_VECPHYS ==1
class VecPhysOrchestrator;
class VecPhysCompton;
#endif

#include "GeantFwd.h"

class GeantPropagator {
public:
  /**
   * @brief Monitoring type
   */
  enum EGeantMonitoringType {
    kMonQueue = 0,
    kMonMemory,
    kMonBasketsPerVol,
    kMonVectors,
    kMonConcurrency,
    kMonTracksPerEvent,
    kMonTracks
  };
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;
  // data members to be made private
  int fNthreads;                                   /** Number of worker threads */
  int fNevents;                                    /** Number of buffered events */
  int fNtotal;                                     /** Total number of events to be transported */
  std::atomic<long> fNtransported;                 /** Number of transported tracks */
  std::atomic<long> fNprimaries;                   /** Number of primary tracks */
  std::atomic<long> fNsteps;                       /** Total number of steps */
  std::atomic<long> fNsnext;                       /** Total number of calls to getting distance to next boundary */
  std::atomic<long> fNphys;                        /** Total number of steps to physics processes */
  std::atomic<long> fNmag;                         /** Total number of partial steps in magnetic field */
  std::atomic<long> fNsmall;                       /** Total number of small steps taken */
  std::atomic<long> fNcross;                       /** Total number of boundaries crossed */
  std::atomic_flag fFeederLock = ATOMIC_FLAG_INIT; /** Atomic flag to protect the particle feeder */
  std::atomic_int fPriorityEvents;                 /** Number of prioritized events */
  BitSet *fDoneEvents;                             /** Array of bits marking done events */
  int fNprocesses;                                 /** Number of active physics processes */
  int fNstart;                                     /** Cumulated initial number of tracks */
  int fMaxTracks;                                  /** Maximum number of tracks per event */
  int fMaxThreads;                                 /** Maximum number of threads */
  int fNminThreshold;                              /** Threshold for starting transporting a basket */
  int fDebugEvt;                                   /** Event to debug */
  int fDebugTrk;                                   /** Track to debug */
  int fDebugStp;                                   /** Step to start debugging */
  int fDebugRep;                                   /** Number of steps to debug */
  int fMaxSteps;                                   /** Maximum number of steps per track */
  int fNperBasket;                                 /** Number of tracks per basket */
  int fMaxPerBasket;                               /** Maximum number of tracks per basket */
  int fMaxPerEvent;                                /** Maximum number of tracks per event */
  int fMaxDepth;                                   /** Maximum geometry depth */
  int fLearnSteps;                                 /** Number of steps needed for the learning phase */
  int fLastEvent;                                  /** Last transported event */
  float fPriorityThr;                              /** Threshold for prioritizing events */
  int fNstepsKillThr;                              /** Threshold in number of steps to kill a track */
  int fNminReuse; /** Minimum number of transported tracks to be reused without re-basketizing */
  int fNTransportTask;
  double fMaxRes;    /** Maximum resident memory allowed [MBytes] */
  double fMaxVirt;   /** Maximum virtual memory allowed [MBytes] */
  double fNaverage;  /** Average number of tracks per event */
  double fVertex[3]; /** Vertex position */
  double fEmin;      /** Min energy threshold */
  double fEmax;      /** Max energy threshold */
  double fBmag;      /** Magnetic field */
  double fEpsilonRK; /** Relative error in RK integration */

  bool fUsePhysics;            /** Enable/disable physics */
  bool fUseRungeKutta;         /** Enable/disable Runge-Kutta integration in field */
  bool fUseDebug;              /** Use debug mode */
  bool fUseGraphics;           /** Graphics mode */
  bool fUseStdScoring;         /** Use standard scoring */
  bool fTransportOngoing;      /** Flag for ongoing transport */
  bool fSingleTrack;           /** Use single track transport mode */
  bool fFillTree;            /** Enable I/O */
  int fTreeSizeWriteThreshold; /** Maximum size of the tree (before automatic writing) **/
  bool fConcurrentWrite;     /** switch between single and mutlithreaded writing */
  bool fUseMonitoring;         /** Monitoring different features */
  bool fUseAppMonitoring;      /** Monitoring the application */
  std::mutex fTracksLock;          /** Mutex for adding tracks */

  WorkloadManager *fWMgr;             /** Workload manager */
  GeantVApplication *fApplication;    /** User application */
  GeantVApplication *fStdApplication; /** Standard application */
  GeantVTaskMgr     *fTaskMgr;        /** GeantV task manager */

  #ifdef USE_ROOT
  TStopwatch *fTimer; /** Timer */
  #else
  vecgeom::Stopwatch *fTimer; /** Timer */
  #endif

  PhysicsProcess *fProcess;              /** For now the only generic process pointing to the tabulated physics */
#if USE_VECPHYS ==1
    PhysicsProcess *fVectorPhysicsProcess;        /** interface to vector physics final state sampling */
    VecPhysOrchestrator *fVecPhysOrchestrator;     /**Orchestrator of the vectorized physics*/
#endif
    
  //   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
  GeantTrack_v *fStoredTracks;         /** Stored array of tracks (history?) */
  PrimaryGenerator *fPrimaryGenerator; /** Primary generator */
  MCTruthMgr *fTruthMgr;               /** MCTruth manager */

  // Data per event
  int *fNtracks;               /** ![fNevents] Number of tracks {array of [fNevents]} */
  GeantEvent **fEvents;        /** ![fNevents]    Array of events */
  GeantTaskData **fThreadData; /** ![fNthreads] Data private to threads */

  static GeantPropagator *fgInstance;

  /** @brief Initialization function */
  void Initialize();

  /** @brief Initialization function */
  void InitializeAfterGeom();

  /** @brief Initialize classes for RK Integration */
  void PrepareRkIntegration();

  /** @brief Function for loading geometry */
  bool LoadGeometry(const char *filename = "geometry.root");

  /** @brief Function for loading VecGeom geometry */
  bool LoadVecGeomGeometry();

  /** @brief Function to initialize VecGeom navigators */
  void InitNavigators();

public:
  /** @brief GeantPropagator constructor */
  GeantPropagator();

  /** @brief GeantPropagator destructor */
  virtual ~GeantPropagator();

  /**
   * @brief Function that returns the number of prioritized events (C++11)
   * @return Number of prioritized events
   */
  int GetNpriority() const { return fPriorityEvents.load(); }

  /**
   * @brief Function that returns the number of transported tracks (C++11)
   * @return Number of transported tracks
   */
  long GetNtransported() const { return fNtransported.load(); }

  /**
   * @brief Function that returns a temporary track object per thread
   * @details Temporary track for the current caller thread
   *
   * @param tid Track ID
   */
  GeantTrack &GetTempTrack(int tid = -1);

  /**
   * @brief Function to add a track to the scheduler
   *
   * @param track Track that should be added
   */
  int AddTrack(GeantTrack &track);

  /**
   * @brief Function to dispatch a track
   *
   * @param track Track that should be dispatched
   */
  int DispatchTrack(GeantTrack &track, GeantTaskData *td);

  /**
   * @brief  Function for marking a track as stopped
   *
   * @param tracks Track array container
   * @param itr Track id
   */
  void StopTrack(const GeantTrack_v &tracks, int itr);

  /**
   * @brief Feeder for importing tracks
   *
   * @param td Thread data object
   * @param init Flag specifying if this is the first call
   * @return Number of injected baskets
   */
  int Feeder(GeantTaskData *td);

  /** @brief Check if transport is feeding with new tracks. */
  inline bool IsFeeding() {
    bool feeding = fFeederLock.test_and_set(std::memory_order_acquire);
    if (feeding)
      return true;
    fFeederLock.clear(std::memory_order_release);
    return false;
  }

  /**
   * @brief Function for importing tracks
   *
   * @param nevents Number of events
   * @param average Average number of tracks
   * @param startevent Start event
   * @param startslot Start slot
   */
  int ImportTracks(int nevents, int startevent, int startslot, GeantTaskData *td);

  /**
   * @brief Instance function returning the singleton pointer
   *
   * @param ntotal Total number of tracks
   * @param nbuffered Number of buffered tracks
   */
  VECCORE_ATT_HOST_DEVICE
  static GeantPropagator *Instance(int ntotal = 0, int nbuffered = 0, int nthreads = 0);

  /**
   * @brief Propose the physics step for an array of tracks
   *
   * @param ntracks Number of ttracks
   * @param tracks Vector of tracks
   * @param td Thread data
   */
  void ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Apply multiple scattering process
   *
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks
   * @param td Thread data
   */
  void ApplyMsc(int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  //   PhysicsProcess  *Process(int iproc) const {return fProcesses[iproc];}

  /**
   * @brief Getter for the process
   * @return  Generic process pointing to the tabulated physics
   */
  PhysicsProcess *Process() const { return fProcess; }

  /**
   * @brief Entry point to start simulation with GeantV
   *
   * @param geomfile Geometry file
   * @param nthreads Number of threads
   * @param graphics Graphics (by default False)
   * @param single Transport single tracks rather than vectors (by default False)
   * @param Execution using TBB tasks instead of static threads (by default False)
   */
  void PropagatorGeom(const char *geomfile = "geometry.root", int nthreads = 4, bool graphics = false,
                      bool single = false);

  /** @brief Function returning the number of monitored features */
  int GetMonFeatures() const;

  /** @brief Check if a monitoring feature is enabled */
  bool IsMonitored(GeantPropagator::EGeantMonitoringType feature) const;

  /** @brief Enable monitoring a feature */
  void SetMonitored(EGeantMonitoringType feature, bool flag = true);

  /** @brief Setter for the global transport threshold */
  void SetNminThreshold(int thr);

  /** @brief  Getter for task broker */
  TaskBroker *GetTaskBroker();

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief Function checking if transport is completed */
  bool TransportCompleted() const { return ((int)fDoneEvents->FirstNullBit() >= fNtotal); }

  /** @brief Try to acquire the lock */
  bool TryLock() { return (fFeederLock.test_and_set(std::memory_order_acquire)); }

  /** @brief Release the lock */
  void ReleaseLock() { fFeederLock.clear(std::memory_order_release); }

private:
  /** @brief Copy constructor not implemented */
  GeantPropagator(const GeantPropagator &);

  /** @brief Assignment operator not implemented */
  GeantPropagator &operator=(const GeantPropagator &);

};
#endif
