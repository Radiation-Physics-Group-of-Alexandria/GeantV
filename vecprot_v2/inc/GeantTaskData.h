//===--- GeantTaskData.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTaskData.h
 * @brief Implementation of data organized per thread Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKDATA
#define GEANT_TASKDATA

#include <deque>
#include "Geant/Typedefs.h"
#include "GeantTrackVec.h"
#include "GeantPropagator.h"
#include "TrackManager.h"

namespace geantphysics {
  class PhysicsData;
}

#ifdef USE_ROOT
class TRandom;
#endif
#ifdef USE_VECGEOM_NAVIGATOR
#include "base/RNG.h"
#endif

class GUFieldPropagator;
class GUVField;

/**
 * @brief Class GeantTaskData
 * @details Class descripting data organized per thread
 *
 */
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class Basket;
class GeantBasketMgr;
class GeantBasket;
class GeantTrackGeo_v;
class StackLikeBuffer;
class TrackStat;

class GeantTaskData {
public:

  using NumaTrackBlock_t = NumaBlock<GeantTrack, true>;

  GeantPropagator *fPropagator = nullptr; /** GeantPropagator */
  int fTid = -1;         /** Thread unique id */
  int fNode = -1;        /** Locality node */
  size_t fNthreads = 0;  /** Number of transport threads */
  int fMaxDepth = 0;     /** Maximum geometry depth */
  int fSizeBool = 0;     /** Size of bool array */
  int fSizeInt = 0;      /*  Size of int array */
  int fSizeDbl = 0;      /** Size of dbl array */
  bool fToClean = false; /** Flag set when the basket queue is to be cleaned */
  Volume_t *fVolume = nullptr; /** Current volume per thread */
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::RNG *fRndm = nullptr;           /** Random generator for thread */
#elif USE_ROOT
  TRandom *fRndm = nullptr;                /** Random generator for thread */
#endif
  bool *fBoolArray = nullptr;              /** [fSizeBool] Thread array of bools */
  double *fDblArray = nullptr;             /** [fSizeDbl] Thread array of doubles */
  int *fIntArray = nullptr;                /** [fSizeInt] Thread array of ints */
  GeantTrack fTrack;                       /** Track support for this thread */
  VolumePath_t *fPath = nullptr;           /** Volume path for the thread */
  VolumePath_t **fPathV = nullptr;         /** Volume path for the thread */
  VolumePath_t **fNextpathV = nullptr;     /** Volume path for the thread */
  GeantTrackGeo_v *fGeoTrack = nullptr;    /** Geometry track SOA */
  GeantBasketMgr *fBmgr = nullptr;         /** Basket manager collecting mixed tracks */
  GeantBasket *fReused = nullptr;          /** Basket having tracks to be reused in the same volume */
  Basket *fBvector = nullptr;              /** Buffer basket used for vector API */
  Basket *fShuttleBasket = nullptr;        /** Shuttle basket from selectors to follow-up simulation stage */
  vector_t<Basket *> fStageBuffers;        /** Buffers for tracks at input of simulation stages */
  GeantBasket *fImported = nullptr;        /** Basket used to import tracks from the event server */
  StackLikeBuffer *fStackBuffer = nullptr; /** Stack buffer tor this thread */
  TrackStat *fStat = nullptr;              /** Track statictics */
  NumaTrackBlock_t *fBlock = nullptr;      /** Current track block */
  
#ifdef VECCORE_CUDA
  char fPool[sizeof(std::deque<GeantBasket *>)]; // Use the same space ...
  char fBPool[sizeof(std::deque<Basket *>)]; /** Pool of empty baskets */
#else
  std::deque<GeantBasket *> fPool; /** Pool of empty baskets */
  std::deque<Basket *> fBPool; /** Pool of empty baskets */
#endif
  GeantTrack_v  *fTransported;              // Transported tracks in current step
  vector_t<GeantTrack *> fTransported1;     // Transported tracks in current step
  int fNkeepvol = 0;     /** Number of tracks keeping the same volume */
  int fNsteps = 0;       /** Total number of steps per thread */
  int fNsnext = 0;       /** Total number of calls to getting distance to next boundary */
  int fNphys = 0;        /** Total number of steps to physics processes */
  int fNmag = 0;         /** Total number of partial steps in magnetic field */
  int fNpart = 0;        /** Total number of particles transported by the thread */
  int fNsmall = 0;       /** Total number of small steps taken */
  int fNcross = 0;       /** Total number of boundary crossings */
  int fNpushed = 0;      /** Total number of pushes with 1.E-3 */
  int fNkilled = 0;      /** Total number of tracks killed */

  geantphysics::PhysicsData  *fPhysicsData = nullptr; /** Physics data per thread */
  GUVField           *fFieldObj;             // To get value of the field!
  vecgeom::Vector3D<double>  fConstFieldValue;   // Value - if field is constant.
  double                     fBfieldMag;     // Magnitude of field for current track - or for const field
  bool                       fBfieldIsConst; // Flag - is the B field constant ?
  GUFieldPropagator       *fFieldPropagator; // For RK integration of charged particle propagation

private:
   // a helper function checking internal arrays and allocating more space if necessary
  template <typename T> static
   VECCORE_ATT_HOST_DEVICE
  void CheckSizeAndAlloc(T *&array, int &currentsize, size_t wantedsize) {
     if (wantedsize <= (size_t) currentsize)
      return;
    T *newarray = new T[wantedsize];
    memcpy(newarray,array,currentsize*sizeof(T));
    delete[] array;
    array = newarray;
    currentsize = wantedsize;
  }

  /**
   * @brief GeantTaskData constructor based on a provided single buffer.
   */
  VECCORE_ATT_DEVICE
  GeantTaskData(void *addr, size_t nTracks, int maxdepth, int maxPerBasket, GeantPropagator *prop = nullptr);

public:
  /** @brief GeantTaskData constructor */
  GeantTaskData(size_t nthreads, int maxDepth, int maxPerBasket);

  /** @brief GeantTaskData destructor */
  ~GeantTaskData();

  /**
   * @brief GeantTrack MakeInstance based on a provided single buffer.
   */
  VECCORE_ATT_DEVICE
  static GeantTaskData *MakeInstanceAt(void *addr, size_t nTracks, int maxdepth, int maxPerBasket, GeantPropagator *prop);

  /** @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth */
  VECCORE_ATT_DEVICE
  static size_t SizeOfInstance(size_t nthreads, int maxDepth, int maxPerBasket);

  /**
   * @brief Function that return double array
   *
   * @param size Size of double array
   */
  VECCORE_ATT_HOST_DEVICE
  double *GetDblArray(int size) {
    CheckSizeAndAlloc<double>(fDblArray, fSizeDbl, size);
    return fDblArray;
  }

  /**
   * @brief Function that return boolean array
   *
   * @param size Size of boolean array
   */
  VECCORE_ATT_HOST_DEVICE
  bool *GetBoolArray(int size) {
    CheckSizeAndAlloc<bool>(fBoolArray, fSizeBool, size);
    return fBoolArray;
  }

  /**
   * @brief Function that returns int array
   *
   * @param size Size of int array
   */
  VECCORE_ATT_HOST_DEVICE
  int *GetIntArray(int size) {
    CheckSizeAndAlloc<int>(fIntArray, fSizeInt, size);
    return fIntArray;
  }

  /**
   * @brief Function that returns a (per thread/task) preallocated NavigationState object
   *
   */
  VECCORE_ATT_HOST_DEVICE
  VolumePath_t *GetPath() {
    return fPath;
  }

  /**
   * @brief Get the cleared storedtrack
   */
  GeantTrack &GetTrack() {
    fTrack.Clear();
    return fTrack;
  }

  /** @brief Get new track from track manager */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack &GetNewTrack();

#ifndef VECCORE_CUDA
  /**
   * @brief Get next free basket or null if not available
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket();

  /**
   * @brief Recycles a given basket
   *
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b);

  /**
   * @brief Recycles a given basket
   *
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(Basket *b);

  /*
   * @brief Return the size of the basket pool
   *
   */
  size_t GetBasketPoolSize() const { return fPool.size(); }

  /**
   * @brief Function cleaning a number of free baskets
   *
   * @param ntoclean Number of baskets to be cleaned
   * @return Number of baskets actually cleaned
   */
  int CleanBaskets(size_t ntoclean);
#endif

  /** @brief Setter for the toclean flag */
  void SetToClean(bool flag) { fToClean = flag; }

  /** @brief Getter for the toclean flag */
  bool NeedsToClean() const { return fToClean; }
  
  /**
   * @brief Function that returns a temporary track object per task data.
   * @details Temporary track for the current caller thread
   *
   */
  GeantTrack &GetTempTrack() { fTrack.Clear(); return fTrack; }

  /** @brief  Inspect simulation stages */
  void InspectStages(int istage);

private:
  /**
   * @brief Constructor GeantTaskData
   * @todo Still not implemented
   */
  GeantTaskData(const GeantTaskData &);

  /**
   * @brief Operator &operator=
   * @todo Still not implemented
   */
  GeantTaskData &operator=(const GeantTaskData &);

  // ClassDef(GeantTaskData, 1) // Stateful data organized per thread
};
} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
