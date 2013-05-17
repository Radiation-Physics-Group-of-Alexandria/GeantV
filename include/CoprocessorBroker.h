class DevicePtrBase
{
   void *fPtr;
protected:
   void MemcpyToDevice(const void* what, unsigned long nbytes);
   void *GetPtr() const { return fPtr; }
public:
   DevicePtrBase() : fPtr(0) {}
   
   ~DevicePtrBase();

   void Malloc(unsigned long size);
};

template <typename T>
class DevicePtr : public DevicePtrBase
{
public:
   void Alloc(unsigned long nelems = 1) {
      Malloc(nelems*sizeof(nelems));
   }
   void ToDevice(const T* what, unsigned long nelems = 1) {
      MemcpyToDevice(what,nelems*sizeof(T));
   }
   operator T*() const { return reinterpret_cast<T*>(GetPtr()); }
};

class GPFieldMap;
class GPPhysicsTable;
class GPVGeometry;
class GPFieldMap;
class GPPhysicsTable;
struct GXTrack;

#ifndef __CINT__
#include <cuda.h>
#include <curand.h>
#include "random_kernel.h"
#else
class cudaStream_t;
class curandState;
#endif

class GeantTrack;

class CoprocessorBroker
{
public:
   CoprocessorBroker();
   ~CoprocessorBroker();
   
   bool UploadGeometry(GPVGeometry *geom);
   bool UploadMagneticField(GPFieldMap** fieldMap);
   
   bool Upload_eBremTable(const GPPhysicsTable &eBrem_table);
   bool Upload_eIoniTable(const GPPhysicsTable &eIoni_table);
   bool UploadMscTable(const GPPhysicsTable &msc_table);
   
   bool CudaSetup(int nblocks, int nthreads, int maxTrackPerThread);
   void prepateDataArray(long nTracks);

   void runTask(int threadid, int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin);

private:
   char       *fdGeometry; // Point to a GPGeomManager in GPU land.
   DevicePtr<GPFieldMap> fdFieldMap;
   
   DevicePtr<GPPhysicsTable> fd_eBremTable;
   DevicePtr<GPPhysicsTable> fd_eIoniTable;
   DevicePtr<GPPhysicsTable> fd_mscTable;
   
   struct StreamHelper {
      StreamHelper();
      ~StreamHelper();

      bool CudaSetup(int nblocks, int nthreads, int maxTrackPerThread);

      unsigned int TrackToDevice(int tid,
                                 GeantTrack **host_track, int *trackin,
                                 unsigned int startIdx, unsigned int basketSize,
                                 unsigned int chunkSize);

      unsigned int TrackToHost(// int tid,
                               GeantTrack **host_track, int *trackin,
                               unsigned int startIdx);
      
      GXTrack               *fTrack;
      int                   *fPhysIndex;
      int                   *fLogIndex;
      unsigned int           fNStaged;   // How many track have been copied to the scratch area so far.

      cudaStream_t  fStream;
      DevicePtr<curandState> fdRandStates;
      DevicePtr<GXTrack>     fDevTrack;
      DevicePtr<int>         fDevTrackPhysIndex;
      DevicePtr<int>         fDevTrackLogIndex;
      
      void ResetNStaged() { fNStaged = 0; }
      operator cudaStream_t() { return fStream; }
   };

   friend struct TrackToHostData;
   friend void ResetNStaged(cudaStream_t /* stream */, cudaError_t status, void *userData);
   
   StreamHelper fStream0;
   StreamHelper fStream1;
   
   int fNblocks;
   int fNthreads;
   int fNchunk;
   int fMaxTrackPerThread;
   int fKernelType;
   
  
};


