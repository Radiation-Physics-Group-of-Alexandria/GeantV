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
   void GetSecParameters(TNudyElement *targetElement, TNudyEndfRecoPoint *recoPoint);
   double GetCaptureXsec();
   std::vector<double> nCaptureGetEnergy();
   std::vector<double> nCaptureGetcosAng();
   std::string GetCaptureProcessName();
   std::vector<std::string> nCaptureGetsecParticles();
   std::string residueName(int , int );
   void nCaptureProcessReset();
private:
  void GetSecParameter(TNudyElement *, TNudyEndfRecoPoint *recoPoint);
  void FillHisto(double icosLab, double isecEnergyLab);
  std::vector<double> crs;
  double kineticE;
  double cosCM = 0, cosLab = 0, secEnergyCM = 0, secEnergyLab = 0;
  double cosCMC = 0, cosLabC = 0, secEnergyCMC = 0, secEnergyLabC = 0;
  double residueKineticEnergy = 0.0;
  double residuecosAng = 0.0;
  std::vector<double> nCaptureProductsEnergy;
  std::vector<double> nCaptureProductscosAng;
  static const int n1=1000000;
  double sigmaPartial;
  double sigmaTotal;
  double ratio;
  double sum1;
  int residueACa = 0, residueZCa = 0;
  
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
  std::string  secParticleCapt;
  std::string  name;
  std::string  procNameCap;
  std::string  nameC;
  std::string residueType;
  std::vector<std::string>ProcName;
  std::vector<std::string>productsName;
  static const char fkElNam[119][4];
  double mass;
  double charge;
  double sigmaCapture;
  int mt1;
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
