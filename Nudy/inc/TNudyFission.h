#ifndef TNudyFission_H
#define TNudyFission_H
#include "TNudyEndfRecoPoint.h"
#include "Particle.h"
#include "ElementProp.h"
#include "TNudyEndfMat.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#include <string>
#include <fstream>
#include "TNudyElement.h"
#define PI acos(-1.0)
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

class TNudyFission {
public:
   TNudyFission();
   TNudyFission(int ElemId, const char*);
   TNudyFission(ElementProp *, TNudyEndfRecoPoint *recoPoint);
   virtual ~TNudyFission();
   void nFissionXsec(double, TNudyElement * );
   double GetFissionXsec();
   std::vector<double> GetKineticEnergyNeutron();
   std::vector<double> GetcosAngNeutron();
   double GetPromptneutron();
   std::vector<int> GetFissionFragmentsmass();
   std::vector<int> GetFissionFragmentscharge();
   std::string GetsecParticleName( );
   std::string GetFissionProcessName();
   std::string GetParticleFission();
   void processReset();
private:
  double kineticE;
  double cosCMF = 0, cosLab = 0, secEnergyCMF = 0, secEnergyLab = 0;
  double sigmaPartial;
  double sigmaTotal;
  double fissionfragmentsCharge = 0;
  
  std::vector<int> massFissionFragments;
  std::vector<int> chargeFissionFragments;
  
  std::vector<double> secEnergyLabF;
  std::vector<double> cosLabF;
  
  double fissionFragmentmass1 = 0.0, fissionFragmentmass2 = 0.0;
  
  double totalkineticEnergyFF = 0.0;
  double kineticEnergyFragment1 = 0.0, kineticEnergyFragment2=0;
  
  int mtValues;
  int elemId;
  int ielemId;
  int LCT, MF, MT, MF4, MF5, MF6;
  int mt;
  double sigmainElastic=0.0;
  std::string  secParticleFiss;
  std::string  procNameFis;
  std::string  nameF;
  double mass;
  double charge;
  double sigmaFission;
  int mt1;
  double nup;
  double nud;
  double fissHeat;
  
  //TNudyElement *elementProp;
  TNudyEndfRecoPoint *rp;
 // ofstream fout;
 #ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyFission, 1) // class for TNudyFission
#endif
};
#endif
