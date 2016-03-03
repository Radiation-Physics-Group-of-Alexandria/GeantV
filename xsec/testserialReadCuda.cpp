
/*
#include "TTabPhysMgr.h"
#include "GeantPropagator.h"
#ifdef USE_ROOT
#include "TGeoManager.h"
#endif
*/

#ifdef USE_VECGEOM_NAVIGATOR
#include "base/RNG.h"
using vecgeom::RNG;
#define UNIFORM() RNG::Instance().uniform()
#elif USE_ROOT
#include <TRandom.h>
#define UNIFORM() gRandom->Uniform()
#else
#define UNIFORM() ((double)rand())/RAND_MAX
#endif

#include <iostream>
#include <fstream>
#include "backend/cuda/Interface.h"
//#include "TEXsec.h"
//#include "TEFstate.h"
//#include "TPDecay.h"
//void launchExpandPhysicsOnDevice(vecgeom::DevicePtr<char>&, int nBlocks, int nThreads, int nrep, vecgeom::DevicePtr<double> devIPart, vecgeom::DevicePtr<double> devIEnergy);
void launchExpandPhysicsOnDevice(vecgeom::DevicePtr<char>&, int nBlocks, int nThreads, double* iSampled,int* devIPart, float* devIEnergy);
/*
void expandPhysicsLocal(char *buf) {
   std::cout << "Rebuilding TPartIndex store" << std::endl;
   TPartIndex::I()->RebuildClass(buf);
   int sizet = TPartIndex::I()->SizeOf();
   std::cout << "Number of bytes for TPartIndex " << sizet << std::endl;
   buf += sizet;
   std::cout << "Rebuilding x-sec store" << std::endl;
   TEXsec::RebuildStore(buf);
   int sizex = TEXsec::SizeOfStore();
   std::cout << "Number of bytes for x-sec " << sizex << std::endl;
   buf += sizex;
   std::cout << "Rebuilding decay store" << std::endl;
   TPDecay *dec = (TPDecay *) buf;
   dec->RebuildClass();
   TEFstate::SetDecayTable(dec);
   int sized = dec->SizeOf();
   std::cout << "Number of bytes for decay " << sized << std::endl;
   buf += sized;
   std::cout << "Rebuilding final state store" << std::endl;
   TEFstate::RebuildStore(buf);
   int sizef = TEFstate::SizeOfStore();
   std::cout << "Number of bytes for final state -- HOST --" << sizef << std::endl;
}
*/
int main()
{
   char *hostBuf=nullptr;
   int totsize;
   // read from file
   std::ifstream fin("xfphys.bin", std::ios::binary);
   fin.read(reinterpret_cast<char*>(&totsize), sizeof(totsize));
   //buf = new char[totsize];
   hostBuf = (char*)_mm_malloc(totsize,sizeof(double));
   fin.read(reinterpret_cast<char*>(hostBuf), totsize);
   fin.close();
   vecgeom::DevicePtr<char> devBuf;
   devBuf.Allocate(totsize);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR ALLOC buffer\n");
      return 0;
   } 
   devBuf.ToDevice(hostBuf,totsize);
   if (cudaSuccess!=cudaGetLastError()) {
      printf("ERROR MEMCPY buffer\n");
      return 0;
   }

   printf("Total size of store %d\n", totsize);

   const unsigned int nrep = 50;
   double *iSampled;
   float *iEnergy; 
   int *iPart; 
   cudaMallocHost(&iSampled, nrep*sizeof(double));
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating iSample on Host: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMallocHost(&iEnergy, nrep*sizeof(float));
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating iEnergy on Host: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMallocHost(&iPart, nrep*sizeof(int));
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating iPart on Host: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }

   printf("Number of repetitions %d\n", nrep);
   

   for(auto irep=0; irep<nrep; irep++) {
      iSampled[irep] = (double) UNIFORM();
   }
/*
   vecgeom::DevicePtr<double> devIPart;
   vecgeom::DevicePtr<double> devIEnergy;

   devIPart.Allocate(nrep);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR ALLOC buffer\n");
      return 0;
   } 

   devIEnergy.Allocate(nrep);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR ALLOC buffer\n");
      return 0;
   }

 
   devIPart.ToDevice(iSampled,nrep);
 */
    double  *devISampled;
    int  *devIPart;
    float  *devIEnergy;

    cudaMalloc((void**) &devISampled, nrep*sizeof(double)); 
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating devIPart to Device: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
    cudaMalloc((void**) &devIPart, nrep*sizeof(int)); 
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating devIPart to Device: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMalloc((void**) &devIEnergy, nrep*sizeof(float)); 
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating devIEnergy to Device: %s\n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
 
    cudaMemcpy(devISampled,iSampled, nrep*sizeof(double),cudaMemcpyHostToDevice);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR copying iSampled to devIPart on Device:%s \n", cudaGetErrorString(cudaGetLastError()));
      return 0;
   }

  launchExpandPhysicsOnDevice(devBuf, 1, 1,devISampled, devIPart,devIEnergy);



cudaThreadSynchronize();
    cudaMemcpy(iSampled,devISampled, nrep*sizeof(double),cudaMemcpyDeviceToHost);
   cudaError_t error=cudaGetLastError(); 
   if (error != cudaSuccess) {
      printf(" ERROR copy iSampled from Device: %s\n", cudaGetErrorString(error));
      return 0;
   }

    cudaMemcpy(iPart,devIPart, nrep*sizeof(int),cudaMemcpyDeviceToHost);
    error=cudaGetLastError(); 
   if (error != cudaSuccess) {
      printf(" ERROR copy iPart from Device: %s\n", cudaGetErrorString(error));
      return 0;
   }
    cudaMemcpy(iEnergy,devIEnergy, nrep*sizeof(float),cudaMemcpyDeviceToHost);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR copy iEnergy from Device\n");
      return 0;
   }
    

   
/*
   CudaCopyFromDevice(iEnergy,devIEnergy, devIEnergy.SizeOf());
   CudaCopyFromDevice(iPart,devIPart, devIPart.SizeOf());
*/

   printf("============ BACK TO THE HOST ==============\n");
   for(auto irep=0; irep<nrep; irep++) 
      printf(" energy for iPart %d is %f\n",iPart[irep], iEnergy[irep]);

/*
   expandPhysicsLocal(hostBuf);
   const char *fxsec = "/dev/null";
   const char *ffins = "/dev/null";
   #ifdef USE_ROOT
   GeantPropagator::Instance(1,1,1);
   TGeoManager *geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");

   #endif
   TTabPhysMgr::Instance(fxsec, ffins );

   constexpr int nrep = 1000;

   #ifndef USE_VECGEOM_NAVIGATOR
   #ifdef USE_ROOT
   gRandom->SetSeed(12345);
   #else
   srand(12345);
   #endif
   #endif
*/

   std::ofstream fftest("xphysR.txt");

//   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
   int iel =1;
      for(auto irep=0; irep<nrep; ++irep) {
/*
	 int ipart = iPart[irep] * TPartIndex::I()->NPartReac();
         int ireac = iPart[irep] * FNPROC;
	 float en =  iPart[irep] * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
          //cout<<"using RNG "<<ipart<<endl;
         float xs = TEXsec::Element(iel)->XS(ipart, ireac, en);
 	 if(xs < 0) continue;
	 int npart=0;
	 float weight=0;
	 float kerma=0;
	 float enr=0;
	 const int *pid=0;
	 const float *mom=0;
	 int ebinindx=0;
	 TEFstate::Element(iel)->SampleReac(ipart, ireac, en, npart, weight, kerma, enr, pid, mom, ebinindx);
	 if(npart <= 0) continue;
         printf(" energy HOST %d\n",enr);
	 fftest <<  iel << ":" << TPartIndex::I()->PartName(ipart) << ":" << ireac << ":" << en
		<< ":" << xs << ":" << npart << ":" << weight << ":" << kerma << ":" << enr << ":";
	 for(auto i=0; i<npart; ++i)
	    fftest << pid[i] << ":" << mom[i*3] << ":" << mom[i*3+1] << ":" << mom[i*3+2];
	 fftest <<":" << ebinindx << std::endl;
*/
	 fftest <<  iPart[irep] << ":" << iEnergy[irep] <<std::endl;
      }
  /// }
   fftest.close();
/*
   #ifdef USE_ROOT
   delete geom;
   #endif
*/
   return 0;
}
