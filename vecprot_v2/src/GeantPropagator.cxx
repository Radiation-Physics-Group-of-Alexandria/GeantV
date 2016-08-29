// A simple propagator taking as input a set of particles located in a given
// volume AND the global matrix of the volume.
// The ProposeStep() method choses between a "scattering" process with no eloss
// and a "ionization" process and generates a random "physical" step. In this simple
// model all particlea undertake the same list of processes
// The ScatteringProcess() method emulates scattering and changes the particle
// direction randomly in a forward cone with opening angle proportional with 1/p
// The IonizationProcess() method simulates an energy deposition with an amount
// epsil*int(1+K*rnd) (epsil, 2*epsil, ..., K*epsil)
// In future we can add absorption and decay to check the stack model...
// The PropagateInField(step) method propagates in a uniform magnetic field with
// an amount equal to step
// The Transport() method computes the safety and snext and compares with the
// physical step. If (safety>pstep), PropagateInField(pstep) is called and the
// physics process is simulated. Otherwise, PropagateInField(safety) is called
// and safety subtracted from the physical step, then the procedure is repeated
// until C*snext/4 < 1E-6 (tangent of angle with sagita, C=1/R is the curvature)
//
#include "GeantPropagator.h"

#ifdef USE_ROOT
#include "TTimer.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include <fenv.h>
#endif


#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NewSimpleNavigator.h"
#include "navigation/SimpleABBoxNavigator.h"
#include "navigation/SimpleABBoxLevelLocator.h"
#include "navigation/HybridNavigator2.h"
#ifdef USE_ROOT
#include "management/RootGeoManager.h"
#endif
#include "volumes/PlacedVolume.h"
#else
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#endif

#include "Geant/Error.h"
#include "GeantTrackVec.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "StdApplication.h"
#include "GeantFactoryStore.h"
#include "GeantEvent.h"
#include "GeantScheduler.h"
#include "PrimaryGenerator.h"
#include "MCTruthMgr.h"

#ifdef USE_CALLGRIND_CONTROL
#include <valgrind/callgrind.h>
#endif

// #define RUNGE_KUTTA  1

// The classes for integrating in a non-uniform magnetic field
// #ifdef   RUNGE_KUTTA
#include "TUniformMagField.h"
#include "FieldEquationFactory.h"
#include "StepperFactory.h"
#include "GUIntegrationDriver.h"

#include "GUFieldPropagator.h"
#include "GUFieldPropagatorPool.h"

#include "UserDetectorConstruction.h"
#include "FieldPropagatorFactory.h"
// #endif

using namespace Geant;
using namespace vecgeom;

GeantPropagator *gPropagator = 0;
GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
    : fNthreads(1), fNevents(100), fNtotal(1000), fNtransported(0), fNprimaries(0), fNsteps(0),
      fNsnext(0), fNphys(0), fNmag(0), fNsmall(0), fFeederLock(ATOMIC_FLAG_INIT), fPriorityEvents(0), fDoneEvents(0),
      fNprocesses(3), fNstart(0), fMaxTracks(0), fMaxThreads(100), fNminThreshold(10), fDebugEvt(-1), fDebugTrk(-1),
      fDebugStp(-1), fDebugRep(-1), fMaxSteps(10000), fNperBasket(16), fMaxPerBasket(256), fMaxPerEvent(0),
      fMaxDepth(0), fLearnSteps(0), fLastEvent(0), fPriorityThr(0), fMaxRes(0), fMaxVirt(0), fNaverage(0), fVertex(),
      fEmin(1.E-4), // 100 KeV
      fEmax(10),    // 10 Gev
      fEpsilonRK(0.0003),
      // fBmag(0.),
      fUsePhysics(true), fUseRungeKutta(false), fInitialisedRKIntegration(false), fUseDebug(false), fUseGraphics(false), fUseStdScoring(false),
      fTransportOngoing(false), fSingleTrack(false), fFillTree(false), fTreeSizeWriteThreshold(100000), fConcurrentWrite(true), fUseMonitoring(false), fUseAppMonitoring(false), fTracksLock(),  
      fWMgr(0), fApplication(0), fStdApplication(0), fTimer(0), fProcess(0), fVectorPhysicsProcess(0), fStoredTracks(0),
      fPrimaryGenerator(0), fNtracks(0), fEvents(0), fThreadData(0) {
  // Constructor
  fVertex[0] = fVertex[1] = fVertex[2] = 0.;
  fBfieldArr[0] = fBfieldArr[1] = fBfieldArr[2] = 0.;

  fgInstance = this;
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator() {
  // Destructor
  int i;
  delete fProcess;
  BitSet::ReleaseInstance(fDoneEvents);
#if USE_VECPHYS == 1
  delete fVectorPhysicsProcess;
#endif

  if (fEvents) {
    for (i = 0; i < fNevents; i++)
      delete fEvents[i];
    delete[] fEvents;
  }

  if (fThreadData) {
    for (i = 0; i < fNthreads; i++)
      delete fThreadData[i];
    delete[] fThreadData;
  }
  delete fTimer;
  delete fWMgr;
  delete fApplication;
}

//______________________________________________________________________________
int GeantPropagator::AddTrack(GeantTrack &track) {
  // Add a new track in the system. returns track number within the event.
  int slot = track.fEvslot;
  track.fParticle = fEvents[slot]->AddTrack();
  
  // call MCTruth manager if it has been instantiated
  if(fTruthMgr) fTruthMgr->AddTrack(track);
  
  //   fNtracks[slot]++;
  fNtransported++;
  return track.fParticle;
}

//______________________________________________________________________________
int GeantPropagator::DispatchTrack(GeantTrack &track, GeantTaskData *td) {
  // Dispatch a registered track produced by the generator.
  return fWMgr->GetScheduler()->AddTrack(track, td);
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(const GeantTrack_v &tracks, int itr) {
  // Mark track as stopped for tracking.
  //   Printf("Stopping track %d", track->particle);

  // stoping track in MCTruthManager
  if(fTruthMgr)
    {
      if(tracks.fStatusV[itr] == kKilled) fTruthMgr->EndTrack(tracks, itr);
    }
  
  if (fEvents[tracks.fEvslotV[itr]]->StopTrack())
    fPriorityEvents++;
}

//______________________________________________________________________________
GeantTrack &GeantPropagator::GetTempTrack(int tid) {
  // Returns a temporary track support for the physics processes, unique per
  // thread which can be used to add tracks produced by physics processes.
  if (tid < 0)
    tid = WorkloadManager::Instance()->ThreadId();
  if (tid > fNthreads)
    Geant::Fatal("GeantPropagator::GetTempTrack", "Thread id %d is too large (max %d)", tid, fNthreads);
  GeantTrack &track = fThreadData[tid]->fTrack;
  track.Clear();
  return track;
}

//______________________________________________________________________________
int GeantPropagator::Feeder(GeantTaskData *td) {
  // Feeder called by any thread to inject the next event(s)
  // Only one thread at a time
  if (fFeederLock.test_and_set(std::memory_order_acquire))
    return -1;
  int nbaskets = 0;
  if (!fLastEvent) {
    nbaskets = ImportTracks(fNevents, 0, 0, td);
    fLastEvent = fNevents;
    fFeederLock.clear(std::memory_order_release);
    return nbaskets;
  }
  // Check and mark finished events
  for (int islot = 0; islot < fNevents; islot++) {
    GeantEvent *evt = fEvents[islot];
    if (fDoneEvents->TestBitNumber(evt->GetEvent()))
      continue;
    if (evt->Transported()) {
      fPriorityEvents--;
      evt->Print();
      
      // closing event in MCTruthManager
      if(fTruthMgr) fTruthMgr->CloseEvent(evt->GetEvent());
      
      // Digitizer (todo)
      int ntracks = fNtracks[islot];
   #ifdef USE_ROOT
      Printf("= digitizing event %d with %d tracks pri=%d", evt->GetEvent(), ntracks, fPriorityEvents.load());
   #else
      printf("= digitizing event %d with %d tracks pri=%d", evt->GetEvent(), ntracks, fPriorityEvents.load());
   #endif
      //            propagator->fApplication->Digitize(evt->GetEvent());
      fDoneEvents->SetBitNumber(evt->GetEvent());
      if (fLastEvent < fNtotal) {
      #ifdef USE_ROOT
        Printf("=> Importing event %d", fLastEvent);
      #else
        printf("=> Importing event %d", fLastEvent);
      #endif
        nbaskets += ImportTracks(1, fLastEvent, islot, td);
        fLastEvent++;
      }
    }
  }

  fFeederLock.clear(std::memory_order_release);
  return nbaskets;
}

//______________________________________________________________________________
int GeantPropagator::ImportTracks(int nevents, int startevent, int startslot, GeantTaskData *thread_data) {
// Import tracks from "somewhere". Here we just generate nevents.
#ifdef USE_VECGEOM_NAVIGATOR
  using vecgeom::SimpleNavigator;
  using vecgeom::Vector3D;
  using vecgeom::Precision;
  using vecgeom::GeoManager;
#endif

  Volume_t *vol = 0;
  int ntracks = 0;
  int ntotal = 0;
  int ndispatched = 0;
  GeantTaskData *td = thread_data;
  if (td == 0) {
    int tid = WorkloadManager::Instance()->ThreadId();
    td = fThreadData[tid];
    td->fTid = tid;
  }
  VolumePath_t *startpath = VolumePath_t::MakeInstance(fMaxDepth);
  GeantBasketMgr *basket_mgr = 0;

#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::SimpleNavigator nav;
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif

  static bool init = true;
  if (init)
    init = false;
  int event = startevent;
  GeantEventInfo eventinfo;
  for (int slot = startslot; slot < startslot + nevents; slot++) {
    eventinfo = fPrimaryGenerator->NextEvent();
// Set initial track states
#ifdef USE_VECGEOM_NAVIGATOR
    startpath->Clear();
    nav.LocatePoint(GeoManager::Instance().GetWorld(),
                    Vector3D<Precision>(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert), *startpath, true);
    vol = const_cast<Volume_t *>(startpath->Top()->GetLogicalVolume());
    basket_mgr = static_cast<GeantBasketMgr *>(vol->GetBasketManagerPtr());
#else
    TGeoNode *node = nav->FindNode(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert);
    vol = node->GetVolume();
    basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
    startpath->InitFromNavigator(nav);
#endif
    basket_mgr->SetThreshold(fNperBasket);
    td->fVolume = vol;
    ntracks = eventinfo.ntracks;
    ntotal += ntracks;
    fNprimaries += ntracks;
    if (!fEvents[slot])
      fEvents[slot] = new GeantEvent();
    fEvents[slot]->SetSlot(slot);
    fEvents[slot]->SetEvent(event);
    fEvents[slot]->Reset();
    // Set priority threshold to non-default value
    if (fPriorityThr > 0)
      fEvents[slot]->SetPriorityThr(fPriorityThr);

    // start new event in MCTruthMgr
    if(fTruthMgr) fTruthMgr->OpenEvent(fEvents[slot]->GetEvent());

    for (int i = 0; i < ntracks; i++) {
      GeantTrack &track = td->GetTrack();
      track.SetPath(startpath);
      track.SetNextPath(startpath);
      track.SetEvent(event);
      track.SetEvslot(slot);
      fPrimaryGenerator->GetTrack(i, track);
      if (!track.IsNormalized())
        track.Print("Not normalized");
      track.fBoundary = false;
      track.fStatus = kAlive;
      AddTrack(track);
      ndispatched += DispatchTrack(track, thread_data);
    }
    event++;
  }

  VolumePath_t::ReleaseInstance(startpath);
  Geant::Print("ImportTracks", "Imported %d tracks from events %d to %d. Dispatched %d baskets.", ntotal, startevent,
               startevent + nevents - 1, ndispatched);
  return ndispatched;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeantPropagator *GeantPropagator::Instance(int ntotal, int nbuffered, int nthreads) {
  // Single instance of the propagator
  if (fgInstance)
    return fgInstance;
  if (ntotal <= 0 || nbuffered <= 0) {
#ifdef USE_ROOT
    Printf("GeantPropagator::Instance: Number of transported/buffered events should be positive");
#else
    printf("GeantPropagator::Instance: Number of transported/buffered events should be positive");
#endif
    return 0;
  }
  fgInstance = new GeantPropagator();
  fgInstance->fNtotal = ntotal;
  fgInstance->fNevents = nbuffered;
  fgInstance->fNthreads = nthreads;
  if (nbuffered > ntotal) {
#ifdef USE_ROOT
    Printf("GeantPropagator::Instance: Number of buffered events changed to %d", ntotal);
#else
    printf("GeantPropagator::Instance: Number of buffered events changed to %d", ntotal);
#endif
    fgInstance->fNevents = ntotal;
  }
  // Initialize workload manager
  fgInstance->fWMgr = WorkloadManager::Instance(nthreads);
  // Instantiate factory store
  GeantFactoryStore::Instance(nbuffered);
  return fgInstance;
}

//______________________________________________________________________________
void GeantPropagator::Initialize() {
  // Initialization
  fMaxPerEvent = 5 * fNaverage;
  fMaxTracks = fMaxPerEvent * fNevents;

  // Initialize arrays here.
  gPropagator = GeantPropagator::Instance();
  fDoneEvents = BitSet::MakeInstance(fNtotal);
  if (!fProcess) {
    Geant::Fatal("GeantPropagator::Initialize", "The physics process has to be initialized before this");
    return;
  }
  // Initialize the process(es)
  fProcess->Initialize();
#if USE_VECPHYS == 1
  fVectorPhysicsProcess->Initialize();
#endif

  if (!fNtracks) {
    fNtracks = new int[fNevents];
    memset(fNtracks, 0, fNevents * sizeof(int));
  }
}

//______________________________________________________________________________
void GeantPropagator::InitializeAfterGeom() {
  // Initialization, part two.

  // Add some empty baskets in the queue
  fWMgr->CreateBaskets(); // geometry must be created by now

  if( !fUserDetectorCtion ){
    Printf("- GeantPropagator::InitializeAfterGeom - %s.\n",
           " no User Detector Construction found." );
    Printf("     Default object created with field= %f %f %f\n",
            fBfieldArr[0], fBfieldArr[1], fBfieldArr[2]
          );
    fUserDetectorCtion= new UserDetectorConstruction();
    fUserDetectorCtion->UseConstantMagField( fBfieldArr );

    // CMSDetectorConstruction* CMSdetector= new CMSDetectorConstruction();
    // CMSdetector->SetFileForField("CMSmagneticField.txt");
    // fUserDetectorCtion= CMSdetector;
  }

  if( !fInitialisedRKIntegration ) {
    GeantPropagator::PrepareRkIntegration();
  }

  if (!fThreadData) {
    fThreadData = new GeantTaskData *[fNthreads];
    for (int i = 0; i < fNthreads; i++) {
      fThreadData[i] = new GeantTaskData(fNthreads, fMaxDepth, fMaxPerBasket);
      fThreadData[i]->fTid = i;
      if( fUseRungeKutta ){
        // delete fThreadData[i]->fFieldPropagator;  // -> No, do not own it
        fThreadData[i]->fFieldPropagator = ObtainThreadRkPropagator(i);
      }
      fThreadData[i]->fFieldObj = ObtainField(i);
      fThreadData[i]->fBfieldIsConst = fUserDetectorCtion->IsFieldUniform();
    }
  }
  // Initialize application
  if (fUseStdScoring) {
    fStdApplication = new StdApplication();
    fStdApplication->Initialize();
  }
  fApplication->Initialize();
}

void GeantPropagator::PrepareRkIntegration() {

  // using GUFieldPropagatorPool = ::GUFieldPropagatorPool;
  Printf(" GV-Propagator calling Detector's CreateFieldAndSolver");
  Printf("    fUserDetectorCtion= %p ", fUserDetectorCtion );
  
  // bool createdField= 
  fUserDetectorCtion->CreateFieldAndSolver(fUseRungeKutta);
  
  Printf("GV-Propagator CreateFieldAndSolver called.\n");

  // Create clones for other threads
  if( fUseRungeKutta )
    GUFieldPropagatorPool::Instance()->Initialize(fNthreads);

  fInitialisedRKIntegration= true;
}

#ifdef USE_ROOT
//______________________________________________________________________________
void GeantPropagator::InitNavigators() {
#if USE_VECGEOM_NAVIGATOR == 1
  for (auto &lvol : GeoManager::Instance().GetLogicalVolumesMap()) {
    if (lvol.second->GetDaughtersp()->size() < 4) {
      lvol.second->SetNavigator(NewSimpleNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 5) {
      lvol.second->SetNavigator(SimpleABBoxNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 10) {
      lvol.second->SetNavigator(HybridNavigator<>::Instance());
      HybridManager2::Instance().InitStructure((lvol.second));
    }
    lvol.second->SetLevelLocator(SimpleABBoxLevelLocator::GetInstance());
  }
#endif
}

GUFieldPropagator *
GeantPropagator::ObtainThreadRkPropagator(unsigned int ThreadId) {
  GUFieldPropagator *fieldPropagator= 0;

  Printf("GeantPropagator::ObtainThreadRkPropagator called with tid=%d\n",
          ThreadId );

  if( ! fInitialisedRKIntegration ) {
    Geant::Error("PrepareThreadRkPropagator",
            "Not initialised RkIntegration before calling GeantPropagator::PrepareThreadRkPropagator");
  }
  if( fInitialisedRKIntegration && fUseRungeKutta ){
    static GUFieldPropagatorPool* fpPool= GUFieldPropagatorPool::Instance();
    assert( fpPool );

    fieldPropagator = fpPool->GetPropagator(ThreadId);
    assert( fieldPropagator );
  }
  return fieldPropagator;
}

GUVField* GeantPropagator::ObtainField(unsigned int threadId) {
  GUVField* pField = nullptr;
  if( fInitialisedRKIntegration && fUseRungeKutta ){
    pField= GUFieldPropagatorPool::Instance()->GetField(threadId);
  } else {
    // Create a uniform field ?
    // pField = new TUniformMagField( vecgeom::Vector3D<float> zeroField );
  }
  return pField;
}

/**
 * function to setup the VecGeom geometry from a TGeo geometry ( if gGeoManager ) exists
 */
//______________________________________________________________________________
bool GeantPropagator::LoadVecGeomGeometry() {
#ifdef USE_VECGEOM_NAVIGATOR
  if (vecgeom::GeoManager::Instance().GetWorld() == NULL) {
    printf("Now loading VecGeom geometry\n");
    vecgeom::RootGeoManager::Instance().LoadRootGeometry();
    printf("Loading VecGeom geometry done\n");
    printf("Have depth %d\n", vecgeom::GeoManager::Instance().getMaxDepth());
    std::vector<vecgeom::LogicalVolume *> v1;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(v1);
    printf("Have logical volumes %ld\n", v1.size());
    std::vector<vecgeom::VPlacedVolume *> v2;
    vecgeom::GeoManager::Instance().getAllPlacedVolumes(v2);
    printf("Have placed volumes %ld\n", v2.size());
    //    vecgeom::RootGeoManager::Instance().world()->PrintContent();

    if (fWMgr->GetTaskBroker())
      Printf("Now upload VecGeom geometry to Coprocessor(s)\n");
    return fWMgr->LoadGeometry();
  }
  if (fWMgr && fWMgr->GetTaskBroker()) {
    printf("Now upload VecGeom geometry to Coprocessor(s)\n");
    return fWMgr->LoadGeometry();
  }
  InitNavigators();
  return true;
#else
  return false;
#endif
}
#endif

//______________________________________________________________________________
bool GeantPropagator::LoadGeometry(const char *filename) {
// Load the detector geometry from file, unless already loaded.

#ifdef USE_ROOT
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
#else
  TGeoManager *geom = (gGeoManager) ? gGeoManager : TGeoManager::Import(filename);
#endif
  if (geom) {
#ifdef USE_VECGEOM_NAVIGATOR
    LoadVecGeomGeometry();
    fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
#else
    fMaxDepth = TGeoManager::GetMaxLevels();
#endif
    return true;
  }
  Geant::Error("GeantPropagator::LoadGeometry", "Cannot load geometry from file %s", filename);
  return false;
#else

  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
     geom->LoadGeometryFromSharedLib(filename);
     fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
     return true;
  }
  return false;
#endif
}

//______________________________________________________________________________
void GeantPropagator::ApplyMsc(int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Apply multiple scattering for charged particles.
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
#else
    mat = td->fVolume->GetMaterial();
#endif
  fProcess->ApplyMsc(mat, ntracks, tracks, td);
}

//______________________________________________________________________________
void GeantPropagator::ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Generate all physics steps for the tracks in trackin.
  // Reset the current step length to 0
  for (int i = 0; i < ntracks; ++i) {
    tracks.fStepV[i] = 0.;
    tracks.fEdepV[i] = 0.;
  }
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
  ;
#else
    mat = td->fVolume->GetMaterial();
#endif
  fProcess->ComputeIntLen(mat, ntracks, tracks, 0, td);
}

//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(const char *geomfile, int nthreads, bool graphics, bool single) {
  // Propagate fNevents in the volume containing the vertex.
  // Simulate 2 physics processes up to exiting the current volume.
  static bool called = false;
  fUseGraphics = graphics;
  fNthreads = nthreads;
  fSingleTrack = single;
  std::cout << " GeantPropagator::PropagatorGeom called with app= " << fApplication << std::endl;
  if (!fApplication) {
    printf("No user application attached - aborting");
    return;
  }

  //  Read and create field early - to check, as it is failing ... 
  if( !fInitialisedRKIntegration ) {
    Printf("- GeantPropagator::PropagatorGeom> Reading and creating field early - to check"); 
    GeantPropagator::PrepareRkIntegration();
    fInitialisedRKIntegration= true;
  }
  
  Initialize();
  // Initialize geometry and current volume
  if (!LoadGeometry(geomfile))
    return;
  InitializeAfterGeom();
  // Initialize application
  fApplication->Initialize();
  if (called) {
    printf("Sorry, you can call this only once per session.");
    return;
  }
  called = true;

  fPrimaryGenerator->InitPrimaryGenerator();
  //   int itrack;

  if (fSingleTrack)
    Printf("==== Executing in single track loop mode using %d threads ====", fNthreads);
  else
    Printf("==== Executing in vectorized mode using %d threads ====", fNthreads);
  if (fUsePhysics)
    Printf("  Physics ON with %d processes", fNprocesses);
  else
    Printf("  Physics OFF");
  if (fUseRungeKutta)
    Printf("  Runge-Kutta integration ON with epsilon= %g", fEpsilonRK );
  else
    Printf("  Runge-Kutta integration OFF");

  // Import the input events. This will start also populating the main queue
  if (!fEvents) {
    fEvents = new GeantEvent *[fNevents];
    memset(fEvents, 0, fNevents * sizeof(GeantEvent *));
  }

  //  Feeder(fThreadData[0]);

  // Loop baskets and transport particles until there is nothing to transport anymore
  fTransportOngoing = true;
  WorkloadManager::Instance()->SetMaxThreads(nthreads);
#ifdef USE_ROOT
  if (fUseMonitoring) {
    TCanvas *cmon = new TCanvas("cscheduler", "Scheduler monitor", 900, 600);
    cmon->Update();
  }
  if (fUseAppMonitoring) {
    TCanvas *capp = new TCanvas("capp", "Application canvas", 700, 800);
    capp->Update();
  }

  fTimer = new TStopwatch();
#ifdef USE_CALLGRIND_CONTROL
  CALLGRIND_START_INSTRUMENTATION;
#endif
#else
  fTimer = new vecgeom::Stopwatch();
#endif

  fWMgr->StartThreads();
  fTimer->Start();
  // Wake up the main scheduler once to avoid blocking the system
  //  condition_locker &sched_locker = fWMgr->GetSchLocker();
  //  sched_locker.StartOne();
  fWMgr->WaitWorkers();
  fWMgr->JoinThreads();
  fTimer->Stop();
#ifdef USE_CALLGRIND_CONTROL
  CALLGRIND_STOP_INSTRUMENTATION;
  CALLGRIND_DUMP_STATS;
#endif
#ifdef USE_ROOT
  double rtime = fTimer->RealTime();
  double ctime = fTimer->CpuTime();
#else
  double rtime = fTimer->Elapsed();
  double ctime = fTimer->Elapsed();  // TO FIX 
#endif
  //   fTimer->Print();
  double speedup = ctime / rtime;
  double efficiency = speedup / nthreads;
//   fWMgr->Print();

#ifdef GEANTV_OUTPUT_RESULT_FILE
  const char *geomname = geomfile;
  if (strstr(geomfile, "http://root.cern.ch/files/"))
    geomname = geomfile + strlen("http://root.cern.ch/files/");
#endif
  //  int nsteps = fWMgr->GetScheduler()->GetNsteps();
  Printf("=== Transported: %ld primaries/%ld tracks,  total steps: %ld, snext calls: %ld, "
         "phys steps: %ld, mag. field steps: %ld, small steps: %ld bdr. crossings: %ld  RT=%gs, CP=%gs",
         fNprimaries.load(), fNtransported.load(), fNsteps.load(), fNsnext.load(), fNphys.load(), fNmag.load(),
         fNsmall.load(), fNcross.load(), rtime, ctime);
  Printf("   nthreads=%d speed-up=%f  efficiency=%f", nthreads, speedup, efficiency);
  //  Printf("Queue throughput: %g transactions/sec", double(fWMgr->FeederQueue()->n_ops()) / rtime);
  fApplication->FinishRun();
  if (fStdApplication)
    fStdApplication->FinishRun();
#ifdef USE_VECGEOM_NAVIGATOR
  printf("=== Navigation done using VecGeom ====");
#else
  printf("=== Navigation done using TGeo    ====");
#endif
  //  Printf("Navstate pool usage statistics:");
  //   fWMgr->NavStates()->statistics();
}

//______________________________________________________________________________
int GeantPropagator::GetMonFeatures() const {
  // Get the number of monitored features
  return fWMgr->GetMonFeatures();
}

//______________________________________________________________________________
bool GeantPropagator::IsMonitored(EGeantMonitoringType feature) const {
  // Check if a given feature is monitored
  return fWMgr->IsMonitored(feature);
}

//______________________________________________________________________________
void GeantPropagator::SetMonitored(EGeantMonitoringType feature, bool flag) {
  // Enable monitoring a feature
  fWMgr->SetMonitored(feature, flag);
}

//______________________________________________________________________________
void GeantPropagator::SetNminThreshold(int thr) {
  // Setter for the global transport threshold
  fWMgr->SetNminThreshold(thr);
}

//______________________________________________________________________________
void GeantPropagator::SetTaskBroker(TaskBroker *broker) {
  // Setter for task broker
  fWMgr->SetTaskBroker(broker);
}

//______________________________________________________________________________
TaskBroker *GeantPropagator::GetTaskBroker() {
  // Getter for task broker
  return fWMgr->GetTaskBroker();
}
