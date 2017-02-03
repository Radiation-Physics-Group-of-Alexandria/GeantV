#ifndef TNudyElastic_H
#define TNudyElastic_H
#include "TNudyEndfRecoPoint.h"
#include "Particle.h"
#include "ElementProp.h"
#include "TNudyEndfMat.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#include <string>
#include "TNudyElement.h"
#define PI acos(-1.0)
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

class TNudyElastic {
public:
   TNudyElastic();
   TNudyElastic(int ElemId, const char*);
   TNudyElastic(ElementProp *, TNudyEndfRecoPoint *recoPoint);
   virtual ~TNudyElastic();
   void nElasticXsec(double, TNudyElement *);
   double GetElasticXsec();
   std::vector<double> GetsecParticleKiEn();
   std::vector<double> GetsecParticlecosAng();
   std::string GetElasticProcessName();
   std::vector<std::string> GetParticleElastic();
   void nElasticProcessReset();
private:
  void GetSecParameter(TNudyElement *, TNudyEndfRecoPoint *recoPoint);
  std::string residueNameElastic(int , int );
  static const char fkElNam[119][4];
  double kineticE;
  double cosCME = 0, cosLabE = 0, secEnergyCME = 0, secEnergyLabE = 0;
  double secondarycosAng = 0;
  double secondarycosAng2 = 0;
  double sigmaPartial;
  double sigmaTotal;
  int residueAEl, residueZEl;
  int mtValues;
  int elemId;
  int ielemId;
  int LCT, MF, MT, MF4, MF5, MF6;
  int mt;
  std::string  secParticleElas;
  std::string  processName;
  std::string  secParticleName;
  std::vector<std::string>  productsName;
  std::vector<double> energyProductsLab;
  std::vector<double> cosAngProductsLab;
  std::string residueType;
  double mass;
  double charge;
  double sigmaElastic;
  
  //TNudyElement *elementProp;
  TNudyEndfRecoPoint *rp;
 // ofstream fout;
 #ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyElastic, 1) // class for TNudyElastic
#endif
};
#endif
