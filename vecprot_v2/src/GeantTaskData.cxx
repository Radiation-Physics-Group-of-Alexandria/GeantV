#include "GeantTaskData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"

#include "TArrayI.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "volumes/LogicalVolume.h"
typedef vecgeom::LogicalVolume TGeoVolume;
#else
#include "TGeoVolume.h"
#endif
#include "TRandom.h"

#include "base/SOA3D.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::GeantTaskData(Int_t nthreads, Int_t maxDepth, Int_t maxPerBasket)
<<<<<<< HEAD
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false), fVolume(0), fRndm(nullptr),
      fBoolArray(0), fDblArray(0), fTrack(0, maxDepth), fPath(0), fBmgr(0), fPool() {
=======
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(0), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0,maxDepth),
      fPath(0), fBmgr(0), fPool(),
      fSOA3Dworkspace1( new vecgeom::SOA3D<vecgeom::Precision>(5*maxPerBasket) ),
      fSOA3Dworkspace2( new vecgeom::SOA3D<vecgeom::Precision>(5*maxPerBasket) ),
      fSizeInt( 5*maxPerBasket ),
      fIntArray( new int[fSizeInt] )
      {
>>>>>>> thread specific workspace data for vector navigation + minor refactoring
  // Constructor
  fNthreads = nthreads;
  fMaxDepth = maxDepth;
  fSizeBool = fSizeDbl = 5 * maxPerBasket;
  fBoolArray = new Bool_t[fSizeBool];
  fDblArray = new Double_t[fSizeDbl];
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
#ifndef GEANT_NVCC
  fRndm = new TRandom();
#endif
}

//______________________________________________________________________________
GeantTaskData::GeantTaskData()
<<<<<<< HEAD
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false), fVolume(0), fRndm(nullptr),
      fBoolArray(0), fDblArray(0), fTrack(0), fPath(0), fBmgr(0), fPool() {
=======
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(0), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0),
      fPath(0), fBmgr(0), fPool(),
      fSOA3Dworkspace1(),
      fSOA3Dworkspace2(),
      fSizeInt(0),
      fIntArray(nullptr) {
>>>>>>> thread specific workspace data for vector navigation + minor refactoring
  // Constructor
  GeantPropagator *propagator = GeantPropagator::Instance();
  fNthreads = propagator->fNthreads;
  fMaxDepth = propagator->fMaxDepth;
  fSizeBool = fSizeDbl = fSizeInt = 5 * propagator->fMaxPerBasket;
  fBoolArray = new Bool_t[fSizeBool];
  fDblArray = new Double_t[fSizeDbl];
  fIntArray = new int[fSizeInt];
  fSOA3Dworkspace1 = new vecgeom::SOA3D<double>( fSizeInt );
  fSOA3Dworkspace2 = new vecgeom::SOA3D<double>( fSizeInt );
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
  fRndm = new TRandom();
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::~GeantTaskData() {
// Destructor
//  delete fMatrix;
#ifndef GEANT_NVCC
  delete fRndm;
#endif
  delete[] fBoolArray;
  delete[] fDblArray;
  delete[] fIntArray;
  delete fSOA3Dworkspace1;
  delete fSOA3Dworkspace2;
  VolumePath_t::ReleaseInstance(fPath);
}

#ifndef GEANT_NVCC
//______________________________________________________________________________
GeantBasket *GeantTaskData::GetNextBasket() {
  // Gets next free basket from the queue.
  if (fPool.empty()) return nullptr;
  GeantBasket *basket = fPool.back();
  //  basket->Clear();
  fPool.pop_back();
  return basket;
}

//______________________________________________________________________________
void GeantTaskData::RecycleBasket(GeantBasket *b) {
  // Recycle a basket.
  fPool.push_back(b);
}

//______________________________________________________________________________
Int_t GeantTaskData::CleanBaskets(size_t ntoclean) {
  // Clean a number of recycled baskets to free some memory
  GeantBasket *b;
  Int_t ncleaned = 0;
  size_t ntodo = 0;
  if (ntoclean == 0)
    ntodo = fPool.size() / 2;
  else
    ntodo = TMath::Min(fPool.size(), ntoclean);
  for (size_t i = 0; i < ntodo; i++) {
    b = fPool.back();
    delete b;
    ncleaned++;
    fPool.pop_back();
  }
  fToClean = false;
  //  Printf("Thread %d cleaned %d baskets", fTid, ncleaned);
  return ncleaned;
}

#endif

} // GEANT_IMPL_NAMESPACE
} // geant
