#ifndef TNudyCapture_H
#define TNudyCapture_H
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


class TNudyCapture {

public:
  TNudyCapture();
  TNudyCapture(int ElemId, const char*);
  TNudyCapture(ElementProp *, TNudyEndfRecoPoint *recoPoint);
  virtual ~TNudyCapture();
   void nCaptureXsec(double, TNudyElement *);
   double GetCaptureXsec();
   double GetKiEnCapture();
   double GetcosAngCapture();
   std::string GetCaptureProcessName();
   std::string GetParticleCapture();
   
private:
  void GetSecParameter(TNudyElement *, TNudyEndfRecoPoint *recoPoint);
  void FillHisto(double icosLab, double isecEnergyLab);
  std::vector<double> crs;
  double kineticE;
  double cosCM = 0, cosLab = 0, secEnergyCM = 0, secEnergyLab = 0;
  double cosCMC = 0, cosLabC = 0, secEnergyCMC = 0, secEnergyLabC = 0;
  
  static const int n1=1000000;
  double sigmaPartial;
  double sigmaTotal;
  double ratio;
  double sum1;
  double residueACa, residueZCa;
  
  int mtValues;
  int elemId;
  int ielemId;
  int isel     = 0;
  int counter  = 0;
  int ecounter = 0;
  int LCT, MF, MT, MF4, MF5, MF6;
  int events;
  int mZA;
  
 int mt;
 double sigmainElastic=0.0;
  rowd Eout;
  rowd costheta;
  std::string  secParticleElas;
  std::string  secParticleFiss;
  std::string  secParticleCapt;
  std::string  secParticleInel;
  std::string  name;
  
  std::string  procNameEla;
  std::string  procNameFis;
  std::string  procNameCap;
  std::string  procNameIne;
  std::string  nameE;
  std::string  nameF;
  std::string  nameC;
  std::string  nameI;
  std::vector<std::string>ProcName;
  double mass;
  double charge;
  double sigmaElastic,sigmaFission,sigmaCapture,sigmaInelastic;
 
  
   int mt1;
   std::vector<int> mtTemp;
   std::vector<double> crssIel;
   std::vector<double> crss;
  
  
 
  //TNudyElement *elementProp;
  TNudyEndfRecoPoint *rp;
 // ofstream fout;
 #ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyCapture, 1) // class for TNudyCapture
#endif
};
#endif
