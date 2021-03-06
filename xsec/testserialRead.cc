#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"
#include <iostream>
#include <fstream>

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

void expandPhysics(char *buf);

using std::cout;
using std::endl;

int main()
{
   char *buf=nullptr;
   int totsize;
   // read from file
   std::ifstream fin("xfphys.bin", std::ios::binary);
   if ( ! fin.is_open() ) {
     std::cerr << "Could not open file xfphys.bin in current directory\n";
     return 1;
   }
   fin.read(reinterpret_cast<char*>(&totsize), sizeof(totsize));
   //buf = new char[totsize];
   buf = (char*)_mm_malloc(totsize,sizeof(double));
   if (!buf) {
     std::cerr << "Failed to allocate memory, needed: " << totsize*sizeof(double) << '\n';
     return 2;
   }
   fin.read(reinterpret_cast<char*>(buf), totsize);
   fin.close();
   std::cout << "Total size of store " << totsize << std::endl;

   expandPhysics(buf);

   constexpr int nrep = 1000;

   #ifndef USE_VECGEOM_NAVIGATOR
   #ifdef USE_ROOT
   gRandom->SetSeed(12345);
   #else
   srand(12345);
   #endif
   #endif

   std::ofstream fftest("xphysR.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {

	 int ipart = UNIFORM() * TPartIndex::I()->NPartReac();
         int ireac = UNIFORM() * FNPROC;
	 float en =  UNIFORM() * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
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
	 fftest <<  iel << ":" << TPartIndex::I()->PartName(ipart) << ":" << ireac << ":" << en
		<< ":" << xs << ":" << npart << ":" << weight << ":" << kerma << ":" << enr << ":";
	 for(auto i=0; i<npart; ++i)
	    fftest << pid[i] << ":" << mom[i*3] << ":" << mom[i*3+1] << ":" << mom[i*3+2];
	 fftest <<":" << ebinindx << std::endl;
      }
   }
   fftest.close();
   return 0;
}

void expandPhysics(char *buf) {
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
   std::cout << "Number of bytes for final state " << sizef << std::endl;
}
