#include "trackingTest2_kernel.h"
#include "gxtracking_kernel.h"
#include "sort_kernel.h"

#include "GXTrackLiason.h"

//#include "gxtracking_kernel.cu"
//cuda/kernel/src/trackingTest2_kernel.cu

int electron_gpu(curandState* devStates,
                 size_t nSteps,
                 size_t nElectrons,
                 GXTrack *track, GXTrack * /* altTrack */,
                 int *logVolumeIndices,
                 int *physVolumeIndices,
                 GXTrack *secondaries, int *secStackSize,
                 
                 int * /* scratch */,
                 GXTrackLiason * /* trackScratch */,

                 GPGeomManager *geomManager,
                 GXFieldMap *magMap,
                 GPPhysicsTable *physicsTable,
                 GPPhysics2DVector *seltzerBergerTable,

                 int nBlocks, int nThreads,
                 cudaStream_t stream)
{
   elec_gpu(devStates, track,
            logVolumeIndices, physVolumeIndices,
            nElectrons, 
            geomManager, magMap, 
            physicsTable, seltzerBergerTable,
            secondaries, secStackSize,
            nSteps, 
            0 /* runType */,
            nBlocks, nThreads, stream);
   return 0;
}

int electron_multistage_gpu(curandState* devStates,
                            size_t nSteps,
                            size_t nElectrons,
                            GXTrack *track, GXTrack *altTrack,
                            int *logVolumeIndices,
                            int *physVolumeIndices,
                            GXTrack *secondaries, int *secStackSize,

                            int *scratch, // array of 10.
                            GXTrackLiason *trackScratch,

                            GPGeomManager *geomManager,
                            GXFieldMap *magMap,
                            GPPhysicsTable *physicsTable,
                            GPPhysics2DVector *seltzerBergerTable,
                            
                            int nBlocks, int nThreads,
                            cudaStream_t stream)
{
   int *nbrem = &(scratch[0]);
   int *nioni = &(scratch[1]);
   int *stackSize_brem = &(scratch[2]);
   int *stackSize_ioni = &(scratch[3]);


   cudaMemsetAsync( &(scratch[0]), 0, sizeof(scratch[0])*5, stream);
   cudaMemsetAsync( secStackSize, 0, sizeof(int), stream);

   elec_GPIL_gpu(devStates, track, 
                 trackScratch, 
                 logVolumeIndices, physVolumeIndices,
                 nElectrons,
                 geomManager, magMap, 
                 physicsTable, seltzerBergerTable,
                 secondaries, secStackSize,
                 nSteps,
                 0 /* runType */,
                 nBlocks, nThreads, stream);
   
   // //atomic counter for the last array position of physics processes
   // cudaThreadSynchronize();
   
   count_by_process_gpu(nElectrons, track,
                        nbrem, nioni, nBlocks, nThreads, stream);
   
   //   cudaMemcpyAsync(&nbrem_h,nbrem_d,sizeof(G4int),cudaMemcpyDeviceToHost,stream);
   //   cudaMemcpyAsync(&nioni_h,nioni_d,sizeof(G4int),cudaMemcpyDeviceToHost,stream);
   
   // This copies the track from track to altTrack
   sort_by_process_gpu(nElectrons, track, altTrack,
                       nbrem, stackSize_brem,
                       nioni, stackSize_ioni,
                       nBlocks, nThreads, stream);
   
   // cudaThreadSynchronize();
   
   elec_doit_gpu(devStates, altTrack, 
                 trackScratch, 
                 nElectrons, 
                 geomManager, magMap, 
                 physicsTable,seltzerBergerTable,
                 secondaries, secStackSize,
                 nSteps, 
                 0 /* runType */,
                 nBlocks, nThreads, stream);

   return 1; // The real data is in altTrack
}