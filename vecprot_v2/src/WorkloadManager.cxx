#include "WorkloadManager.h"
#include "Geant/Error.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "GeantTrack.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantTaskData.h"
#include "PhysicsProcess.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"
#include "GeantVApplication.h"
#if USE_VECGEOM_NAVIGATOR
#include "base/TLS.h"
#include "management/GeoManager.h"
#include "materials/Medium.h"
#else
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#endif
#include "TaskBroker.h"

// added by WP for output handling
#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif
#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif
#include "TThread.h"
#include "GeantFactoryStore.h"
#include "TThreadMergingFile.h"

using namespace Geant;
using std::max;

WorkloadManager *WorkloadManager::fgInstance = 0;

//______________________________________________________________________________
WorkloadManager::WorkloadManager(int nthreads)
    : fNthreads(nthreads), fNbaskets(0), fBasketGeneration(0), fNbasketgen(0), fNidle(nthreads), fNminThreshold(10),
      fNqueued(0), fBtogo(0), fSchId(nthreads), fStarted(false), fStopped(false), fFeederQ(0), fTransportedQ(0),
      fDoneQ(0), fListThreads(), fFlushed(false), fFilling(false), fMonQueue(0), fMonMemory(0), fMonBasketsPerVol(0),
      fMonVectors(0), fMonConcurrency(0), fMonTracksPerEvent(0), fMonTracks(0), fMaxThreads(0), fScheduler(0), fBroker(0), fWaiting(0),
      fSchLocker(), fGbcLocker(), fLastEvent(0), fOutputIO(0) {
  // Private constructor.
  fFeederQ = new Geant::priority_queue<GeantBasket *>(1 << 16);
  fTransportedQ = new Geant::priority_queue<GeantBasket *>(1 << 16);
  fDoneQ = new Geant::priority_queue<GeantBasket *>(1 << 10);
  fgInstance = this;
  fScheduler = new GeantScheduler();
  fWaiting = new int[fNthreads + 1];
  memset(fWaiting, 0, (fNthreads + 1) * sizeof(int));

  fOutputIO = new dcqueue<TBufferFile*>;
  fMergingServer = new TThreadMergingServer(fOutputIO);
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager() {
  // Destructor.
  delete fFeederQ;
  delete fTransportedQ;
  delete fDoneQ;
  delete fScheduler;
  delete[] fWaiting;
  //   delete fNavStates;
  fgInstance = 0;
  delete fOutputIO;
}

//______________________________________________________________________________
int WorkloadManager::ThreadId() {
#if USE_VECGEOM_NAVIGATOR
  return BaseTLS::ThreadId();
#else  
  gGeoManager->SetMultiThread();
  return TGeoManager::ThreadId();
#endif  
}

//______________________________________________________________________________
void WorkloadManager::CreateBaskets() {
  // Create the array of baskets
  VolumePath_t *blueprint = 0;
  int maxdepth = GeantPropagator::Instance()->fMaxDepth;
  Geant::Info("CreateBaskets","Max depth: %d", maxdepth);
  blueprint = VolumePath_t::MakeInstance(maxdepth);
  //   fNavStates = new GeantObjectPool<VolumePath_t>(1000*fNthreads, blueprint);
  //   fNavStates = new rr_pool<VolumePath_t>(16*fNthreads, 1000, blueprint);
  fScheduler->CreateBaskets();
  VolumePath_t::ReleaseInstance(blueprint);

  if (fBroker)
    fBroker->CreateBaskets();
}

//______________________________________________________________________________
WorkloadManager *WorkloadManager::Instance(int nthreads) {
  // Return singleton instance.
  if (fgInstance)
    return fgInstance;
  if (!nthreads) {
    ::Error("WorkloadManager::Instance", "No instance yet so you should provide number of threads.");
    return 0;
  }
  return new WorkloadManager(nthreads);
}

//______________________________________________________________________________
void WorkloadManager::Print(const char *) const {
  //
}

void WorkloadManager::SetTaskBroker(TaskBroker *broker) {
  // Register the broker if it is valid.
  if (broker && broker->IsValid())
    fBroker = broker;
}

#if USE_VECGEOM_NAVIGATOR == 1
//______________________________________________________________________________
bool WorkloadManager::LoadGeometry(vecgeom::VPlacedVolume const *const volume) {
  /**
   * @brief Tell the task broker(s) to load the geometry.
   *
   * @param Volume to load
   */
  if (fBroker)
    return fBroker->UploadGeometry(volume);
  return true;
}
#endif

//______________________________________________________________________________
void WorkloadManager::StartThreads() {
  // Start the threads
  fStarted = true;
  if (!fListThreads.empty())
    return;
  int ith = 0;
  if (fBroker) {
     if (fBroker->GetNstream() > (unsigned int)fNthreads) {
      ::Fatal("StartThreads", "The task broker is using too many threads (%d out of %d)", fBroker->GetNstream(),
              fNthreads);
    }
    Geant::Info("StartThreads", "Running with a coprocessor broker (using %d threads).",fBroker->GetNstream()+1);
    fListThreads.emplace_back(WorkloadManager::TransportTracksCoprocessor, fBroker);
    ith += fBroker->GetNstream() + 1;
    if (ith == fNthreads && fBroker->IsSelective()) {
       Fatal("WorkloadManager::StartThreads","All %d threads are used by the coprocessor broker but it can only process a subset of particles.",fNthreads);
       return;
    }
  }
  // Start CPU transport threads
  for (; ith < fNthreads; ith++) {
    fListThreads.emplace_back(WorkloadManager::TransportTracks);
  }
  // Start output thread
  if (GeantPropagator::Instance()->fFillTree) {
    fListThreads.emplace_back(WorkloadManager::OutputThread);
  }
  // Start monitoring thread
  if (GeantPropagator::Instance()->fUseMonitoring) {
    fListThreads.emplace_back(WorkloadManager::MonitoringThread);
  }
  // Start garbage collector
  if (GeantPropagator::Instance()->fMaxRes > 0) {
    fListThreads.emplace_back(WorkloadManager::GarbageCollectorThread);
  }
}

//______________________________________________________________________________
void WorkloadManager::JoinThreads() {
  //
  int tojoin = fNthreads;
  if (fBroker)
    tojoin -= fBroker->GetNstream();
  for (int ith = 0; ith < tojoin; ith++)
    fFeederQ->push(0);

  for (auto &t : fListThreads) {
    t.join();
  }
}

//______________________________________________________________________________
void WorkloadManager::WaitWorkers() {
  // Waiting point for the main thread until work gets done.
  fFilling = false;
  int ntowait = fNthreads;
  if (fBroker)
    ntowait -= fBroker->GetNstream();
  GeantBasket *signal;
  while (ntowait) {
    fDoneQ->wait_and_pop(signal);
    ntowait--;
    Geant::Print("", "=== %d workers finished", fNthreads - ntowait);
  }
  //   fBasketGeneration++;
}

//______________________________________________________________________________
void *WorkloadManager::MainScheduler(void *) {
  // Garbage collector thread, called by a single thread.
  Geant::Print("","=== Scheduler: stopping threads and exiting ===");
  return 0;
}

//______________________________________________________________________________
static inline void MaybeCleanupBaskets(GeantTaskData *td) {
  if (td->NeedsToClean())
    td->CleanBaskets(0);
  else {
    // Check if there are at least 10 free baskets for this thread
    if (td->fPool.size() < 10) {
       td->fBmgr->CreateEmptyBaskets(10, td);
    }
  }
}

//______________________________________________________________________________
static inline void MaybeCleanupBaskets(GeantTaskData *td, GeantBasket *basket) {
  if (td->NeedsToClean())
    td->CleanBaskets(0);
  else {
    // Check if there are at least 10 free baskets for this thread
    if (td->fPool.size() < 10) {
      if (basket->GetBasketMgr())
        basket->GetBasketMgr()->CreateEmptyBaskets(2, td);
      else
        td->fBmgr->CreateEmptyBaskets(10, td);
    }
  }
}

//______________________________________________________________________________
/** @brief Call Feeder (if needed) and check exit condition. */
WorkloadManager::FeederResult WorkloadManager::CheckFeederAndExit(GeantBasketMgr &prioritizer,
                                                  GeantPropagator &propagator,
                                                  GeantTaskData &td) {

  if (!prioritizer.HasTracks() && (propagator.GetNpriority() || GetNworking() == 1)) {
    bool didFeeder = propagator.Feeder(&td);
    // Check exit condition
    if (propagator.TransportCompleted()) {
      int nworkers = propagator.fNthreads;
      for (int i = 0; i < nworkers; i++)
        FeederQueue()->push(0);
      TransportedQueue()->push(0);
      Stop();
      //         sched_locker.StartOne(); // signal the scheduler who has to exit
      //         gbc_locker.StartOne();
      return FeederResult::kStopProcessing;
    }
    if (didFeeder) return FeederResult::kFeederWork;
  }
  return FeederResult::kNone;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracks() {
  // Thread propagating all tracks from a basket.
  //      char slist[256];
  //      TString sslist;
  //   const int max_idle = 1;
  //   int indmin, indmax;
  static std::atomic<int> counter(0);
  int ntotnext, ncross, nbaskets;
  int ntotransport;
  int nextra_at_rest = 0;
  int generation = 0;
  //   int ninjected = 0;
  int nnew = 0;
  int ntot = 0;
  int nkilled = 0;
  int nphys = 0;
  int nout = 0;
  int ngcoll = 0;
  GeantBasket *basket = 0;
  int tid = Instance()->ThreadId();
  Geant::Print("","=== Worker thread %d created ===", tid);
  GeantPropagator *propagator = GeantPropagator::Instance();
  Geant::GeantTaskData *td = propagator->fThreadData[tid];
  td->fTid = tid;
  int nworkers = propagator->fNthreads;
  WorkloadManager *wm = WorkloadManager::Instance();
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  int *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = new GeantBasketMgr(sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(propagator->fNperBasket);
  prioritizer->SetFeederQueue(feederQ);

  // IO handling

  bool concurrentWrite = GeantPropagator::Instance()->fConcurrentWrite && GeantPropagator::Instance()->fFillTree;
  int treeSizeWriteThreshold = GeantPropagator::Instance()->fTreeSizeWriteThreshold;
    
  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  GeantFactory<MyHit> *myhitFactory = factoryStore->GetFactory<MyHit>(16);

  TThread t;
  TThreadMergingFile* file=0;
  TTree *tree=0;
  GeantBlock<MyHit>* data=0;
  
  if (concurrentWrite)
    {
      file = new TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");     
      tree = new TTree("Tree","Simulation output");
      
      tree->Branch("hitblockoutput", "GeantBlock<MyHit>", &data);
      
      // set factory to use thread-local queues
      myhitFactory->queue_per_thread = true;
    }
  
  // Start the feeder
  propagator->Feeder(td);
  Material_t *mat = 0;
  int *waiting = wm->GetWaiting();
//  condition_locker &sched_locker = wm->GetSchLocker();
//  condition_locker &gbc_locker = wm->GetGbcLocker();
// The last worker wakes the scheduler
//  if (tid == nworkers-1) sched_locker.StartOne();
//   int nprocesses = propagator->fNprocesses;
//   int ninput, noutput;
//   bool useDebug = propagator->fUseDebug;
//   Geant::Print("","(%d) WORKER started", tid);
// Create navigator if none serving this thread.
#ifdef USE_VECGEOM_NAVIGATOR
// Suppose I do not need a new navigator for VecGeom... otherwise how it ever worked... ?
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
  //   int iev[500], itrack[500];
  // TGeoBranchArray *crt[500], *nxt[500];
  while (1) {
    // Call the feeder if in priority mode
    auto feedres = wm->CheckFeederAndExit(*prioritizer, *propagator, *td);
    if (feedres == FeederResult::kFeederWork) {
       ngcoll = 0;
    } else if (feedres == FeederResult::kStopProcessing) {
       break;
    }

    waiting[tid] = 1;
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;
    // If prioritizer has work, just do it
    if (prioritizer->HasTracks()) {
      basket = prioritizer->GetBasketForTransport(td);
      ngcoll = 0;
    } else {
      if ((nbaskets < 1) && (!propagator->IsFeeding())) {
        sch->GarbageCollect(td);
        ngcoll++;
      }
      // Too many garbage collections - enter priority mode
      if ((ngcoll > 5) && (wm->GetNworking() <= 1)) {
        ngcoll = 0;
        for (int slot = 0; slot < propagator->fNevents; slot++)
          if (propagator->fEvents[slot]->Prioritize())
            propagator->fPriorityEvents++;
        while ((!sch->GarbageCollect(td, true)) && (feederQ->size_async() == 0))
          ;
      }
      wm->FeederQueue()->wait_and_pop(basket);
    }
    // Check exit condition: null basket in the queue
    if (!basket)
      break;
    waiting[tid] = 0;
    MaybeCleanupBaskets(td,basket);
    ++counter;
    // ntotransport = basket->GetNinput(); // all tracks to be transported
    // ninput = ntotransport;
    int ninput = basket->GetNinput(); // all tracks to be transported
    ntotransport= ninput;
    
    GeantTrack_v &input = basket->GetInputTracks();

// #define ENERGY_BALANCE 1

#ifdef ENERGY_BALANCE    
    double sumEin= 0.0;
    for (int ix = 0; ix < input.GetNtracks(); ix++) {
       // TrackStatus_t status= input.fStatusV[ix];
       // bool       isOutside= input.fPathV[ix]->IsOutside();

       int pdg = input.fPDGV[ix];
#ifdef USE_VECGEOM_NAVIGATOR
       const Particle_t *const &partPDG = &Particle_t::GetParticle(pdg);
#else
       TParticlePDG *partPDG = TDatabasePDG::Instance()->GetParticle(pdg);
#endif
       double     restMass = partPDG->Mass();
       
       double   eFullTrack = input.fEV[ix];
       double       eTrack = eFullTrack ;
       // double         kinE = eTrack - restMass;     // Only kinetic energy
       
       // Quick and dirty rule: produced matter is 'free', antimatter requires 2x rest-mass
       //  Must exclude antiparticles which were input.
       if( (pdg ==  11) || (pdg ==  2212) || (pdg ==  2112) ) {  // e-, proton & neutron
          eTrack -= restMass;                 // Only kinetic energy
       }
       if( (pdg == -11) || (pdg == -2212) || (pdg == -2112) ) {  // e-, proton & neutron
          // 'Anti'-particles generated in the interaction
          eTrack += restMass;                 // + 'Creation' energy
       }
       // eTrack += input.fEdepV[ix];  // Does track enter with Edep > 0 ?
       sumEin += eTrack;
    }
    double sumEkilled= 0.0, sumEexiting= 0.0, sumEphysics= 0.0, sumEdep= 0.0, sumEout= 0.0; // sumEremains= 0.0;
    double eBalance =  0.0;
    double eConv = 1.0e+3; // Print in MeV
#endif
    GeantTrack_v &output = *td->fTransported;

    if (!ntotransport)
      goto finish; // input list empty
    //      Geant::Print("","======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport,
    //      counter.load());
    //      basket->Print();
    //      Geant::Print("","==========================================");
    //      propagator->fTracksPerBasket[tid] = ntotransport;
    
    td->fVolume = 0;
    mat = 0;
    if (!basket->IsMixed()) {
      td->fVolume = basket->GetVolume();
#ifdef USE_VECGEOM_NAVIGATOR
      mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
#else
      mat = td->fVolume->GetMaterial();
#endif
      if (ntotransport < 257)
        nvect[ntotransport] += ntotransport;
    } else {
      nvect[1] += ntotransport;
    }

    // Select the discrete physics process for all particles in the basket
    if (propagator->fUsePhysics) {
      propagator->ProposeStep(ntotransport, input, td);
      // Apply msc for charged tracks
      propagator->ApplyMsc(ntotransport, input, td);
    }

    ncross = 0;
    generation = 0;

    while (ntotransport) {
      // Interrupt condition here. Work stealing could be also implemented here...
      generation++;
      //  Geant::Print("","====== WorkloadManager: Input tracks - top ");
      //     input.PrintTracks();

      // Propagate all remaining tracks
      if (basket->IsMixed())
        ncross += input.PropagateTracksScalar(td);
      else
        ncross += input.PropagateTracks(td);
      ntotransport = input.GetNtracks();
    }
    // All tracks are now in the output track vector. Possible statuses:
    // kCrossing - particles crossing boundaries
    // kPhysics - particles reaching the point where the discrete physics process
    //            will happen.
    // kExitingSetup - particles exiting the geometry
    // kKilled - particles that could not advance in geometry after several tries

    // Geant::Print("","====== WorkloadManager: tracks - type of step selected ");
    // output.PrintTracks();

    // Update track time for all output tracks (this should vectorize)
    nout = output.GetNtracks();
    for (auto itr = 0; itr < nout; ++itr) 
      output.fTimeV[itr] += output.TimeStep(itr, output.fStepV[itr]);
    
    // Post-step actions by continuous processes for all particles. There are no
    // new generated particles at this point.
    if (propagator->fUsePhysics) {
      nphys = 0;
      nextra_at_rest = 0;
      // count phyics steps here
      for (auto itr = 0; itr < nout; ++itr)
        if (output.fStatusV[itr] == kPhysics)
          nphys++;

      propagator->Process()->Eloss(mat, output.GetNtracks(), output, nextra_at_rest, td);
      //         if (nextra_at_rest) Geant::Print("","Extra particles: %d", nextra_at_rest);
      // Now we may also have particles killed by energy threshold
      // Do post-step actions on remaining particles
      // to do: group particles per process

      // Discrete processes only
      nphys = output.SortByLimitingDiscreteProcess(); // only those that kPhysics and not continous limit
      if (nphys) {
        // Do post step actions for particles suffering a given process.
        // Surviving particles are added to the output array

        // first: sample target and type of interaction for each primary tracks
        propagator->Process()->PostStepTypeOfIntrActSampling(mat, nphys, output, td);

        for (auto itr = 0; itr < nout; ++itr) {
           // input.CheckTrack(itr, "Input to PostStepFinalStateSampling.");
           input.Normalize(itr);
           // input.CheckTrack(itr, "Input to PostStepFinalStateSampling - 2nd check.");           
        }
//
// TODO: vectorized final state sampling can be inserted here through
//       a proper interface code that will:
//         1. select the necessary primary tracks based on the alrady
//            sampled interaction type
//         2. take all the member of the selected primary tracks that
//            necessary for sampling the final states
//         3. call the appropriate vector physics code that will
//            perform the physics interaction itself
//         4. update properties of the selected primary tracks,
//            insert secondary tracks into the track vector and set
//            'ntotnext' to the value of the number of secondary tracks
//            inserted to the track vector
//
#if USE_VECPHYS == 1
        propagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
#endif
        // second: sample final states (based on the inf. regarding sampled
        //         target and type of interaction above), insert them into
        //         the track vector, update primary tracks;
        propagator->Process()->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);

        // if (0 /*ntotnext*/) {
        //  Geant::Print("","============= Basket: %s\n", basket->GetName());
        //  output.PrintTracks();        

      }
      // Geant::Print("","============= WrkMgr After UsePhysics - Basket: %s\n", basket->GetName());
      // output.PrintTracks();       
    }
    
#ifdef ENERGY_BALANCE
    sumEkilled= 0.0, sumEexiting= 0.0, sumEphysics= 0.0, sumEdep= 0.0, sumEout= 0.0; // sumEremains= 0.0;

    for (int ix = 0; ix < output.GetNtracks(); ix++) {
         TrackStatus_t status= output.fStatusV[ix];
         bool       isOutside= output.fPathV[ix]->IsOutside();

         int pdg = output.fPDGV[ix];
#ifdef USE_VECGEOM_NAVIGATOR
         const Particle_t *const &partPDG = &Particle_t::GetParticle(pdg);
#else
         TParticlePDG *partPDG = TDatabasePDG::Instance()->GetParticle(pdg);
#endif
         double     restMass = partPDG->Mass();

         double   eFullTrack = output.fEV[ix];
         double       eTrack = eFullTrack ;     // Default: total energy

         // Quick and dirty rule: produced matter is 'free', antimatter requires 2x rest-mass
         //  Must exclude antiparticles which were input.
         if( (pdg ==  11) || (pdg ==  2212) || (pdg ==  2112) ) {
            eTrack -= restMass;                 // Only kinetic energy
         }
         if( (pdg == -11) || (pdg == -2212) || (pdg == -2112) ) {
            // 'Anti'-particles generated in the interaction
            eTrack += restMass;                 // + 'Creation' energy
         }
         // tracks.fStatusV[itr] == kKilled ||
         // tracks.fStatusV[itr] == kExitingSetup ||
         if( status == kKilled ) {
            sumEkilled += eTrack;
            // nKilled ++;
         } else if( status == kExitingSetup || isOutside ) {
            sumEexiting += eTrack;
            // nExiting ++;
         } else if( status == kPhysics ) {
            sumEphysics += eFullTrack;  // Un-adjusted - surviving particle!
            // nInteracting ++;
         } else {
            sumEout    += eTrack;
         }

         sumEdep += output.fEdepV[ix];
      }

      // for (int ix = 0; ix < input.GetNtracks(); ix++) {
      //  if( input.fStatusV[ix] != kKilled ) { sumEremains += input.fEV[ix]; }
      // }

      eConv = 1.0e+3; // Print in MeV
      eBalance=  sumEin - sumEdep - sumEphysics - sumEout - sumEexiting;

      // Geant::Print("WrkMgr",": ...

      if( std::fabs(eBalance) > 0.03 * sumEin ) {
        Printf("WrkMgr: #tracks: in=%4d out= %4d Ein= %13.4f  Edep = %10.4f   Eout= %13.4f  Ekilled = %9.4f  Exiting = %6.3g Physics= %14.6f - Balance= %12.6g",  // Remains= %8.6f 
             // numTracksIn,
             ninput, output.GetNtracks(),
             sumEin*eConv, sumEdep*eConv, sumEout*eConv, sumEkilled*eConv,
             sumEexiting*eConv, sumEphysics*eConv, // sumEremains*eConv,
             eBalance*eConv);
      }

#endif
       // Geant::Print("","====== WorkloadManager: Output tracks - after physics");
       // output.PrintTracks("Output Tracks: ");


    // Report problem tracks (option in comments) -- and correct them!
    if (ntotnext) {
      // Geant::Print("","============= Basket: %s\n", basket->GetName());
      for (int itr=0; itr<ntotnext; ++itr) {
        // output.CheckTrack(itr, "Product of PostStepFinalStateSampling.");
        output.Normalize(itr);
      }
    }

    if (gPropagator->fStdApplication)
      gPropagator->fStdApplication->StepManager(output.GetNtracks(), output, td);
    gPropagator->fApplication->StepManager(output.GetNtracks(), output, td);

    // WP
    if (concurrentWrite) {
      while (!(myhitFactory->fOutputsArray[tid].empty())) {
        data = myhitFactory->fOutputsArray[tid].back();
        myhitFactory->fOutputsArray[tid].pop_back();

        tree->Fill();
        // now we can recycle data memory
        myhitFactory->Recycle(data, tid);
      }
      if (tree->GetEntries() > treeSizeWriteThreshold) file->Write();
    }

    // Update geometry path for crossing tracks
    ntotnext = output.GetNtracks();
    
// Normalize directions (should be not needed)
//    for (auto itr=0; itr<ntotnext; ++itr)
//      output.Normalize(itr);

#ifdef BUG_HUNT
    for (auto itr = 0; itr < ntotnext; ++itr) {
      bool valid = output.CheckNavConsistency(itr);
      if (!valid) {
        valid = true;
      }
    }
    // First breakpoint to be set
    output.BreakOnStep(propagator->fDebugEvt, propagator->fDebugTrk, propagator->fDebugStp, propagator->fDebugRep, "EndStep");
#endif
    for (auto itr = 0; itr < ntotnext; ++itr) {
      output.fNstepsV[itr]++;
      if (output.fStatusV[itr] == kBoundary)
        *output.fPathV[itr] = *output.fNextpathV[itr];
    }
  finish:
    //      basket->Clear();
    //      Geant::Print("","======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
    //      ninjected =
    sch->AddTracks(output, ntot, nnew, nkilled, td);
    //      Geant::Print("","thread %d: injected %d baskets", tid, ninjected);
    //      wm->TransportedQueue()->push(basket);
    //    sched_locker.StartOne();
    // Make sure the basket is not recycled before gets released by basketizer
    // This should not happen vey often, just when some threads are highly
    // demoted ant the basket makes it through the whole cycle before being fully released
    while (basket->fNused.load())
      ;
    basket->Recycle(td);
  }

  // WP
  if (concurrentWrite) {
    file->Write();
  }
   
  wm->DoneQueue()->push(0);
  delete prioritizer;
  // Final reduction of counters
  propagator->fNsteps += td->fNsteps;
  propagator->fNsnext += td->fNsnext;
  propagator->fNphys += td->fNphys;
  propagator->fNmag += td->fNmag;
  propagator->fNsmall += td->fNsmall;
  
  Geant::Print("","=== Thread %d: exiting ===", tid);

  if (wm->IsStopped()) wm->MergingServer()->Finish();

  if (concurrentWrite) {
    delete file;
  }

  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracksCoprocessor(TaskBroker *broker) {
  // Thread propagating all tracks from a basket.
  //      char slist[256];
  //      TString sslist;
  //   const int max_idle = 1;
  //   int indmin, indmax;
  static std::atomic<int> counter(0);
  // int ncross;
  int nbaskets;
  int ntotransport;
  // int nextra_at_rest = 0;
  int generation = 0;
  //   int ninjected = 0;
  int nnew = 0;
  int ntot = 0;
  int nkilled = 0;
  int ngcoll = 0;
  GeantBasket *basket = 0;

  int tid = Instance()->ThreadId();
  Geant::Print("","=== Worker thread %d created for Coprocessor ===", tid);

  GeantPropagator *propagator = GeantPropagator::Instance();
  GeantTaskData *td = propagator->fThreadData[tid];
  td->fTid = tid;
  int nworkers = propagator->fNthreads;
  WorkloadManager *wm = WorkloadManager::Instance();
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  int *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = new GeantBasketMgr(sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(propagator->fNperBasket);
  prioritizer->SetFeederQueue(feederQ);
  // Start the feeder
  propagator->Feeder(td);
  // TGeoMaterial *mat = 0;
  int *waiting = wm->GetWaiting();
//  condition_locker &sched_locker = wm->GetSchLocker();
// int nprocesses = propagator->fNprocesses;
// int ninput;
// int noutput;
//   bool useDebug = propagator->fUseDebug;
//   Geant::Print("","(%d) WORKER started", tid);
#ifdef USE_VECGEOM_NAVIGATOR
// Suppose I do not need a navigator, otherwise how it would have ever worked?
#else
  // Create navigator if none serving this thread.
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
  waiting[tid] = 1;
  // int iev[500], itrack[500];
  // TGeoBranchArray *crt[500], *nxt[500];

  // broker->SetPrioritizer(prioritizer);
  while (1) {

    // Call the feeder if in priority mode
    auto feedres = wm->CheckFeederAndExit(*prioritizer, *propagator, *td);
    if (feedres == FeederResult::kFeederWork) {
       ngcoll = 0;
    } else if (feedres == FeederResult::kStopProcessing) {
       break;
    }

    // ::Info("GPU","Waiting (1) for next available stream.");
    // TaskBroker::Stream stream = broker->GetNextStream();
    // if (!stream)
    //   break;

    if (wm->FeederQueue()->empty()) {
      // There is no work to be done for now, let's just run what we have
      if (0 != broker->launchTask()) {
        // We ran something, let wait for the next free stream,
        MaybeCleanupBaskets(td);
        continue;
      } else {
        // We had nothing to run at all .... we need to wait for
        // more data ....
      }
    }
    waiting[tid] = 1;
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;
    // If prioritizers have work, just do it
    if ((basket = broker->GetBasketForTransport(*td))) {
      ngcoll = 0;
    } else {
      if (nbaskets < 1  && (!propagator->IsFeeding()) ) {
        sch->GarbageCollect(td);
        ngcoll++;
      }
      // Too many garbage collections - enter priority mode
      if ((ngcoll > 5) && (wm->GetNworking() <= 1)) {
        ngcoll = 0;
        for (int slot = 0; slot < propagator->fNevents; slot++)
          if (propagator->fEvents[slot]->Prioritize())
            propagator->fPriorityEvents++;
        while ((!sch->GarbageCollect(td, true)) && (feederQ->size_async() == 0))
          ;
      }
      if (wm->FeederQueue()->try_pop(basket)) {
         // We got a basket, let's go.
      } else {
         // No basket in the queue, let's flush our oustanding work
         if (0 != broker->launchTask(/* wait= */ true)) {
            // We ran something, new basket might be available.
            continue;
         } else {
            // We have nothing, so let's wait.
            wm->FeederQueue()->wait_and_pop(basket);
         }
      }
    }
    // Check exit condition: null basket in the queue
    if (!basket) {
      // All work should have been done, the broker's queue should
      // be empty.
      // assert( worker queue empty );
      break;
    }
    waiting[tid] = 0;
    MaybeCleanupBaskets(td,basket);
    ++counter;

    // if (!stream) {
    //   ::Info("GPU", "Waiting (2) for next available stream.");
    //   stream = broker->GetNextStream();
    // }
    // lastToClear = false;
    // if (!basket) {
    //   if (0 != broker->launchTask(/* wait= */ true)) {
    //     // We ran something, new basket might be available.
    //     continue;
    //   } else {
    //     break;
    //   }
    // }

    ntotransport = basket->GetNinput(); // all tracks to be transported
    // ninput = ntotransport;
    //    GeantTrack_v &input = basket->GetInputTracks();
    GeantTrack_v &output = *td->fTransported;
    if (!ntotransport)
      goto finish; // input list empty
    //      Geant::Print("","======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport,
    //      counter);
    //      basket->Print();
    //      Geant::Print("","==========================================");
    //      propagator->fTracksPerBasket[tid] = ntotransport;
    td->fVolume = 0;
    // mat = 0;
    if (!basket->IsMixed()) {
      td->fVolume = basket->GetVolume();
#ifdef USE_VECGEOM_NAVIGATOR
      // mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
#else
      // mat = td->fVolume->GetMaterial();
#endif
      if (ntotransport < 257)
        nvect[ntotransport] += ntotransport;
    } else {
      nvect[1] += ntotransport;
    }

    // Record tracks
    // ninput = ntotransport;
    //      if (counter==1) input.PrintTracks();
    for (int itr = 0; itr < ntotransport; itr++) {
      // iev[itr] = input.fEventV[itr];
      // itrack[itr] = input.fParticleV[itr];
      // crt[itr] = input.fPathV[itr];
      // nxt[itr] = input.fNextpathV[itr];
      //if (isnan(input.fXdirV[itr])) {
      //  Geant::Error("Before transport", "Error: track %d has NaN", itr);
      //}
    }

    // Select the discrete physics process for all particles in the basket
    // Moved to CoprocessorBroker kernel.
    // if (propagator->fUsePhysics) {
    //   propagator->ProposeStep(ntotransport, input, td);
    //   // Apply msc for charged tracks
    //   propagator->ApplyMsc(ntotransport, input, td);
    // }

    // ncross = 0;
    generation = 0;

    while (ntotransport) {
      // Interrupt condition here. Work stealing could be also implemented here...
      generation++;
      // Propagate all remaining tracks
      // NOTE: need to deal with propagator->fUsePhysics
      broker->runTask(*td, *basket); // ntotransport, basket_sch->GetNumber(), gPropagator->fTracks, particles);
      ntotransport = 0;
      // ncross += input.PropagateTracks(output);
      // ntotransport = input.GetNtracks();
    }
    // All tracks are now in the output track vector. Possible statuses:
    // kCrossing - particles crossing boundaries
    // kPhysics - particles reaching the point where the discrete physics process
    //            will happen.
    // kExitingSetup - particles exiting the geometry
    // kKilled - particles that could not advance in geometry after several tries

    {
       auto noutput = basket->GetNinput();
       for (int itr = 0; itr < noutput; itr++) {
          if (TMath::IsNaN(output.fXdirV[itr])) {
             Geant::Error("TransportTracksCoprocessor","Track %d has NaN", itr);
          }
       }
    }
    // Note: Need to apply the code if (propagator->fUsePhysics)
    // to the tracks after they came back.

  finish:
//      basket->Clear();
//      Geant::Print("","======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
    /* int ninjected = */ sch->AddTracks(output, ntot, nnew, nkilled, td);
    if (prioritizer->HasTracks()) {
       // We know that all the left over tracks injected here can not be handle
       // by the coprocessor, so we much send it for somebody else to process:

       GeantBasket *cbasket = prioritizer->GetCBasket();
       assert(cbasket->TryDispatch());
       prioritizer->SetCBasket(prioritizer->GetNextBasket(td));
       prioritizer->Push(cbasket, /*priority=*/ false, td);
       // prioritizer->BookBasket(td);
       // prioritizer->GarbageCollect(td);
    }
    //      Geant::Print("","thread %d: injected %d baskets", tid, ninjected);
    // wm->TransportedQueue()->push(basket);
    (void)ntot;
    (void)nnew;
    (void)nkilled;
//    sched_locker.StartOne();
    // Make sure the basket is not recycled before gets released by basketizer
    // This should not happen vey often, just when some threads are highly
    // demoted ant the basket makes it through the whole cycle before being fully released
    while (basket->fNused.load()) ;
    basket->Recycle(td);
  }
  // delete prioritizer;
  wm->DoneQueue()->push(0);
  delete prioritizer;
  // Final reduction of counters
  propagator->fNsteps += td->fNsteps;
  propagator->fNsnext += td->fNsnext;
  propagator->fNphys += td->fNphys;
  propagator->fNmag += td->fNmag;
  propagator->fNsmall += td->fNsmall;
  Geant::Print("","=== Coprocessor Thread %d: exiting === Processed %ld", tid, broker->GetTotalWork());
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::GarbageCollectorThread() {
  // This threads can be triggered to do garbage collection of unused baskets
  static double rsslast = 0;
  double rss;
  ProcInfo_t procInfo;
  const double MByte = 1024.;
  const double_t thr_increase = 1.05;
  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  int nthreads = propagator->fNthreads;
  double threshold = propagator->fMaxRes;
  double virtlimit = propagator->fMaxVirt;
  if (threshold == 0 && virtlimit == 0)
    return 0;
  //  condition_locker &gbc_locker = wm->GetGbcLocker();
  while (1) {
    //    gbc_locker.Wait();
    if (wm->IsStopped())
      break;
    // todo: cleaning with thread based recycle queues
    gSystem->GetProcInfo(&procInfo);
    rss = procInfo.fMemResident / MByte;
    double vsize = procInfo.fMemVirtual / MByte;
    Geant::Print("GC thread","### Checking mem: %g MBytes", rss);
    bool needClean = false;
    if (threshold && rss > threshold && (rss / rsslast > thr_increase)) {
      rsslast = rss;
      needClean = true;
    } else if (virtlimit && vsize > virtlimit) {
      needClean = true;
    }
    if (needClean) {
      for (int tid = 0; tid < nthreads; tid++) {
        GeantTaskData *td = propagator->fThreadData[tid];
        td->SetToClean(true);
      }
    }
    gSystem->Sleep(1000); // millisec
  }
  return 0;
}

//______________________________________________________________________________
int WorkloadManager::GetMonFeatures() const {
  // Get the number of monitored features
  return (fMonQueue + fMonMemory + fMonBasketsPerVol + fMonVectors + fMonConcurrency + fMonTracksPerEvent + fMonTracks);
}

//______________________________________________________________________________
bool WorkloadManager::IsMonitored(GeantPropagator::EGeantMonitoringType feature) const {
  // Check if a given feature is monitored
  switch (feature) {
  case GeantPropagator::kMonQueue:
    return fMonQueue;
  case GeantPropagator::kMonMemory:
    return fMonMemory;
  case GeantPropagator::kMonBasketsPerVol:
    return fMonBasketsPerVol;
  case GeantPropagator::kMonVectors:
    return fMonVectors;
  case GeantPropagator::kMonConcurrency:
    return fMonConcurrency;
  case GeantPropagator::kMonTracksPerEvent:
    return fMonTracksPerEvent;
  case GeantPropagator::kMonTracks:
    return fMonTracks;
  }
  return false;
}

//______________________________________________________________________________
void WorkloadManager::SetMonitored(GeantPropagator::EGeantMonitoringType feature, bool flag) {
  // Enable/disable monitoring for a feature
  int value = (int)flag;
  switch (feature) {
  case GeantPropagator::kMonQueue:
    fMonQueue = value;
    break;
  case GeantPropagator::kMonMemory:
    fMonMemory = value;
    break;
  case GeantPropagator::kMonBasketsPerVol:
    fMonBasketsPerVol = value;
    break;
  case GeantPropagator::kMonVectors:
    fMonVectors = value;
    break;
  case GeantPropagator::kMonConcurrency:
    fMonConcurrency = value;
    break;
  case GeantPropagator::kMonTracksPerEvent:
    fMonTracksPerEvent = value;
    break;
  case GeantPropagator::kMonTracks:
    fMonTracks = value;
  }
}

//______________________________________________________________________________
void *WorkloadManager::MonitoringThread() {
  // Thread providing basic monitoring for the scheduler.
  const double MByte = 1024.;
  Geant::Info("MonitoringThread","Started monitoring ...");
  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  int nmon = wm->GetMonFeatures();
  if (!nmon)
    return 0;
  int dmon = 0.5 * nmon + nmon % 2;
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  int ntotransport;
  ProcInfo_t procInfo;
  double rss;
  double nmem[101] = {0};
  GeantScheduler *sch = wm->GetScheduler();
  int nvol = sch->GetNvolumes();
  int *nvect = sch->GetNvect();
  int nthreads = wm->GetNthreads();
  int *nworking = new int[nthreads + 1];
  memset(nworking, 0, (nthreads + 1) * sizeof(int));
  int *waiting = wm->GetWaiting();
  GeantBasketMgr **bmgr = sch->GetBasketManagers();

  TCanvas *cmon = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("cscheduler");
  if (nmon == 1)
    cmon->Divide(1, 1);
  else
    cmon->Divide(2, dmon);
  TH1I *hqueue = 0;
  int nqueue[101] = {0};
  int ipad = 0;
  if (wm->IsMonitored(GeantPropagator::kMonQueue)) {
    hqueue = new TH1I("hqueue", "Work queue load", 100, 0, 100);
    hqueue->SetFillColor(kRed);
    hqueue->SetLineColor(0);
    hqueue->SetStats(false);
    cmon->cd(++ipad);
    hqueue->Draw();
  }
  TH1F *hmem = 0;
  if (wm->IsMonitored(GeantPropagator::kMonMemory)) {
    hmem = new TH1F("hmem", "Resident memory [MB]", 100, 0, 100);
    // hmem->SetFillColor(kMagenta);
    hmem->SetLineColor(kMagenta);
    hmem->SetStats(false);
    cmon->cd(++ipad);
    hmem->Draw();
  }
  TH1I *hbaskets = 0;
  TH1I *hbused = 0;
  if (wm->IsMonitored(GeantPropagator::kMonBasketsPerVol)) {
    hbaskets = new TH1I("hbaskets", "Baskets per volume", nvol, 0, nvol);
    hbaskets->SetFillColor(kBlue);
    hbaskets->SetLineColor(0);
    hbaskets->SetStats(false);
    //    hbaskets->GetYaxis()->SetRangeUser(1,10000);
    hbused = new TH1I("hbused", "Baskets per volume", nvol, 0, nvol);
    hbused->SetFillColor(kRed);
    hbused->SetFillStyle(3001);
    hbused->SetLineColor(0);
    hbused->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    hbaskets->Draw();
    hbused->Draw("SAME");
  }
  TH1I *hvectors = 0;
  if (wm->IsMonitored(GeantPropagator::kMonVectors)) {
    hvectors = new TH1I("hvectors", "Tracks in vectors of given size", 257, 0, 257);
    hvectors->SetFillColor(kBlue);
    hvectors->SetLineColor(0);
    hvectors->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    hvectors->Draw();
  }
  TH1F *hconcurrency = 0;
  TH1F *hconcavg = 0;
  if (wm->IsMonitored(GeantPropagator::kMonConcurrency)) {
    hconcurrency = new TH1F("hconcurrency", "Concurrency plot", nthreads + 1, 0, nthreads + 1);
    hconcurrency->GetYaxis()->SetRangeUser(0, 1);
    hconcurrency->GetXaxis()->SetNdivisions(nthreads + 1, true);
    hconcurrency->GetYaxis()->SetNdivisions(10, true);
    hconcurrency->SetFillColor(kGreen);
    hconcurrency->SetLineColor(0);
    hconcurrency->SetStats(false);
    hconcavg = new TH1F("hconcavg", "Concurrency plot", nthreads + 1, 0, nthreads + 1);
    hconcavg->SetFillColor(kRed);
    hconcavg->SetFillStyle(3001);
    hconcavg->SetLineColor(0);
    hconcavg->SetStats(false);
    cmon->cd(++ipad);
    hconcurrency->Draw();
    hconcavg->Draw("SAME");
  }
  TH1I *htracksmax = 0;
  TH1I *htracks = 0;
  int nbuffered = propagator->fNevents;
  if (wm->IsMonitored(GeantPropagator::kMonTracksPerEvent)) {
    htracksmax = new TH1I("htracksmax", "Tracks in flight", nbuffered, 0, nbuffered);
    htracksmax->SetFillColor(kBlue);
    htracksmax->SetLineColor(0);
    htracksmax->SetStats(false);
    htracks = new TH1I("htracks", "Tracks in flight", nbuffered, 0, nbuffered);
    htracks->SetFillColor(kRed);
    htracks->SetFillStyle(3001);
    htracks->SetLineColor(0);
    htracks->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    htracksmax->Draw();
    htracks->Draw("SAME");
  }
  TH1I *htrackstot = 0;
  int ntrackstot[101] = {0};
  if (wm->IsMonitored(GeantPropagator::kMonTracks)) {
    htrackstot = new TH1I("htrackstot", "Total number of tracks alive", 100, 0, 100);
    htrackstot->SetFillColor(kRed);
    htrackstot->SetLineColor(0);
    htrackstot->SetStats(false);
    cmon->cd(++ipad);
    htrackstot->Draw();
  }
  cmon->Update();
  double stamp = 0.;
  int i, j, bin;
  int nmaxtot;
  while (1) { // exit condition here
    i = int(stamp);
    ipad = 0;
    gSystem->Sleep(50); // millisec
    // Fill histograms
    if (stamp > 100) {
      if (hqueue) {
        ntotransport = feederQ->size_async();
        memmove(nqueue, &nqueue[1], 99 * sizeof(int));
        nqueue[99] = ntotransport;
        hqueue->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          hqueue->SetBinContent(j + 1, nqueue[j]);
      }
      if (hmem) {
        gSystem->GetProcInfo(&procInfo);
        rss = procInfo.fMemResident / MByte;
        memmove(nmem, &nmem[1], 99 * sizeof(double));
        nmem[99] = rss;
        hmem->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          hmem->SetBinContent(j + 1, nmem[j]);
      }
      if (htrackstot) {
        // Count tracks for all event slots
        int ntr = 0;
        for (int slot = 0; slot < propagator->fNevents; slot++)
          ntr += propagator->fEvents[slot]->GetNinflight();
        memmove(ntrackstot, &ntrackstot[1], 99 * sizeof(int));
        ntrackstot[99] = ntr;
        htrackstot->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          htrackstot->SetBinContent(j + 1, ntrackstot[j]);
      }
    } else {
      if (hqueue) {
        ntotransport = feederQ->size_async();
        nqueue[i] = ntotransport;
        hqueue->SetBinContent(i + 1, nqueue[i]);
      }
      if (hmem) {
        gSystem->GetProcInfo(&procInfo);
        rss = procInfo.fMemResident / MByte;
        nmem[i] = rss;
        hmem->SetBinContent(i + 1, nmem[i]);
      }
      if (htrackstot) {
        // Count tracks for all event slots
        int ntr = 0;
        for (int slot = 0; slot < propagator->fNevents; slot++)
          ntr += propagator->fEvents[slot]->GetNinflight();
        ntrackstot[i] = ntr;
        htrackstot->SetBinContent(i + 1, ntrackstot[i]);
      }
    }
    if (hbaskets) {
      for (j = 0; j < nvol; j++) {
        bin = hbaskets->GetXaxis()->FindBin(bmgr[j]->GetVolume()->GetName());
        hbaskets->SetBinContent(bin, bmgr[j]->GetNbaskets());
        hbused->SetBinContent(bin, bmgr[j]->GetNused());
      }
      int nbaskets_mixed = 0;
      int nused_mixed = 0;
      for (j = 0; j < nthreads; j++) {
        GeantTaskData *td = propagator->fThreadData[j];
        nbaskets_mixed += td->fBmgr->GetNbaskets();
        nused_mixed += td->fBmgr->GetNused();
      }
      bin = hbaskets->GetXaxis()->FindBin("MIXED");
      hbaskets->SetBinContent(bin, nbaskets_mixed);
      hbused->SetBinContent(bin, nused_mixed);
    }
    if (hvectors) {
      int vmax = 0;
      int ymax = 0;
      // mixed tracks
      // hvectors->SetBinContent(1, nvect[1]);
      for (j = 1; j < 257; ++j) {
        if (nvect[j] > 0) {
          if (j > vmax)
            vmax = j;
        }
        if (nvect[j] > ymax)
          ymax = nvect[j];
        hvectors->SetBinContent(j + 1, nvect[j]);
      }
      hvectors->GetXaxis()->SetRangeUser(-5, vmax);
      hvectors->GetYaxis()->SetRangeUser(1, 1.2 * ymax);
    }
    if (hconcurrency) {
      for (j = 0; j < nthreads + 1; j++) {
        nworking[j] += 1 - waiting[j];
        hconcurrency->SetBinContent(j + 1, 1 - waiting[j]);
        hconcavg->SetBinContent(j + 1, nworking[j] / (stamp + 1));
      }
    }
    if (htracksmax) {
      nmaxtot = 1;
      for (j = 0; j < nbuffered; j++) {
        GeantEvent *evt = propagator->fEvents[j];
        int nmax = evt->GetNmax();
        nmaxtot = max<int>(nmax, nmaxtot);
        htracksmax->SetBinContent(j + 1, nmax);
        htracks->SetBinContent(j + 1, evt->GetNinflight());
      }
      htracksmax->GetYaxis()->SetRangeUser(1, 10 * nmaxtot);
    }
    // Now draw all histograms
    if (hqueue) {
      cmon->cd(++ipad);
      hqueue->Draw();
    }
    if (hmem) {
      cmon->cd(++ipad);
      hmem->Draw();
    }
    if (hbaskets) {
      cmon->cd(++ipad);
      hbaskets->LabelsOption(">");
      hbaskets->GetXaxis()->SetRangeUser(0, 10);
      hbaskets->Draw();
      hbused->Draw("SAME");
    }
    if (hvectors) {
      cmon->cd(++ipad);
      hvectors->Draw();
    }
    if (hconcurrency) {
      cmon->cd(++ipad);
      hconcurrency->Draw();
      hconcavg->Draw("SAME");
    }
    if (htracksmax) {
      cmon->cd(++ipad);
      htracksmax->GetXaxis()->SetTitle("Event slot");
      htracksmax->Draw();
      htracks->Draw("SAME");
    }
    if (htrackstot) {
      cmon->cd(++ipad);
      htrackstot->Draw();
    }
    cmon->Modified();
    cmon->Update();
    stamp += 1;
    if (wm->IsStopped())
      break;
  }
  delete[] nworking;
  double sumvect = 0.;
  for (i = 0; i < 257; ++i)
    sumvect += nvect[i];
  Geant::Info("Monitoring stopped","=== Percent of tracks transported in single track mode: %g%%", 100. * nvect[1] / sumvect);
  // Sleep a bit to let the graphics finish
  gSystem->Sleep(100); // millisec

  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::OutputThread() {
  // Thread providing basic output for the scheduler.

  Geant::Info("OutputThread", "=== Output thread created ===");

  if (GeantPropagator::Instance()->fConcurrentWrite) {
    Printf(">>> Writing concurrently to MemoryFiles");

    TThread t;

    WorkloadManager::Instance()->MergingServer()->Listen();
  }
  else {
    TFile file("hits.root", "RECREATE");
    WorkloadManager *wm = WorkloadManager::Instance();

    GeantBlock <MyHit> *data = 0;

    TTree *tree = new TTree("Hits", "Simulation output");

    tree->Branch("hitblocks", "GeantBlock<MyHit>", &data);

    GeantFactoryStore *factoryStore = GeantFactoryStore::Instance();
    GeantFactory <MyHit> *myhitFactory = factoryStore->GetFactory<MyHit>(16);


    while (!(wm->IsStopped()) || myhitFactory->fOutputs.size() > 0) {
      // Geant::Print("","size of queue from output thread %zu", myhitFactory->fOutputs.size());

      if (myhitFactory->fOutputs.size() > 0) {
        while (myhitFactory->fOutputs.try_pop(data)) {
          // myhitFactory->fOutputs.wait_and_pop(data);
          // Geant::Print("","Popping from queue of %zu", myhitFactory->fOutputs.size() + 1);
          // if(data) std::cout << "size of the block in the queue " << data->Size() << std::endl;

          tree->Fill();
          myhitFactory->Recycle(data);
        }
      }
    }
    tree->Write();
    file.Close();
  }

  Geant::Info("OutputThread", "=== Output thread finished ===");
    
  return 0;
}
