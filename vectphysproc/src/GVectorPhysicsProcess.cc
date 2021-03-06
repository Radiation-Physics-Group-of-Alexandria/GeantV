
#include "GVectorPhysicsProcess.h"

#include <iostream>

#include "TPartIndex.h"

// Vector physics model related
#include "globals.h"
#include "GeantTrack.h"

#include "GeantTaskData.h"

// Vector physics process related
#include "GVComptonProcess.h"

//------------------------------------------------------------------------------
GVectorPhysicsProcess::GVectorPhysicsProcess()
      : PhysicsProcess(),
        fVComptonProcess(0), fEnergyLimit(0), fNumThreads(0)
{
   // Default type is already discrete.
   // SetType(kDiscrete);
}

//------------------------------------------------------------------------------
GVectorPhysicsProcess::GVectorPhysicsProcess(double energyLimit, int numThreads)
      : PhysicsProcess(),
        fVComptonProcess(0), fEnergyLimit(energyLimit), fNumThreads(numThreads)
{
   // Default type is already discrete.
   // SetType(kDiscrete);
}

//------------------------------------------------------------------------------
GVectorPhysicsProcess::~GVectorPhysicsProcess()
{
  if(fVComptonProcess) 
    for(int i=0; i<fNumThreads; ++i)
     delete fVComptonProcess[i];
  delete [] fVComptonProcess;
}

//------------------------------------------------------------------------------
void GVectorPhysicsProcess::Initialize()
{
// Initialize physics.
  if (fVComptonProcess) return;
  fVComptonProcess = new GVComptonProcess*[fNumThreads];
  for(int i=0; i<fNumThreads; ++i) 
   fVComptonProcess[i] = new GVComptonProcess(TPartIndex::I()->ProcIndex("Compton"), 
                                              fEnergyLimit, 16);
   // 16 is the initial cpacity; cpacity is self adopted during the run 
}

//------------------------------------------------------------------------------
void GVectorPhysicsProcess::PostStepFinalStateSampling( Material_t * /*mat*/,
                                                        int ntracks, 
                                                        GeantTrack_v &tracks,
                                                        int &nout, 
                                                        GeantTaskData *td)
{
    if(!fVComptonProcess) {
       std::cerr<<"\n"
                <<"\t|*********************************************************|\n"
                <<"\t| ERROR in :                                              |\n"
                <<"\t|   GVectorPhysicsProcess::PostStepFinalStateSampling     |\n"
                <<"\t| ------------------------------------------------------- |\n"
                <<"\t| GVectorPhysicsProcess has not been initialized properly:|\n"
                <<"\t|    the GVectorPhysicsProcess::fVComptonProcess member   |\n"
                <<"\t|    is still NULL                                        |\n"
                <<"\t***********************************************************\n"
                <<std::endl;
      exit(EXIT_FAILURE);       
    }

    nout = fVComptonProcess[td->fTid]->ApplyPostStepProcess(tracks, ntracks, td);
}

