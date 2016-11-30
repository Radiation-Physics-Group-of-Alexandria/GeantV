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


typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowint> matrixint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;
typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
typedef std::vector<std::vector<std::vector<std::vector<rowd>>>> matrixd5;


class TNudyElastic {

public:
  TNudyElastic();
  TNudyElastic(int ElemId, const char*);
  
  TNudyElastic(ElementProp *, TNudyEndfRecoPoint *recoPoint);
  virtual ~TNudyElastic();
   
   
   void Neutron_Process(TNudyElement *);
   void nElasticXsec(double, TNudyElement *);
  
   double GetElasticXsec();
   double GetFissionXsec();
   double GetCaptureXsec();
   double GetInElasticXsec();
   double GetsecParticleKiEn();
   double GetsecParticlecosAng();
   
   std::vector<double> GetXsec();
   
   std::string GetsecParticleName( );
   std::string GetElasticProcessName();
   std::string GetFissionProcessName();
   std::string GetCaptureProcessName();
   std::string GetInelasticProcessName();
   std::vector<std::string>GetProcessNames();
   
   std::string GetParticleFission();
   std::string GetParticleElastic();
   std::string GetParticleCapture();
   std::string  GetParticleInelastic();
   
   
   
   
private:
  void GetSecParameter(TNudyElement *, TNudyEndfRecoPoint *recoPoint);
  void FillHisto(double icosLab, double isecEnergyLab);
  std::vector<double> crs;
  double kineticE;
  double cosCM = 0, cosLab = 0, secEnergyCM = 0, secEnergyLab = 0;
  double cosCME = 0, cosLabE = 0, secEnergyCME = 0, secEnergyLabE = 0;
  
  static const int n1=1000000;
  double sigmaPartial;
  double sigmaTotal;
  double ratio;
  double sum1;
  double residueAEl, residueZEl;
  
  int mtValues;
  int elemId;
  int pId;
  int ielemId;
  int isel     = 0;
  int counter  = 0;
  int ecounter = 0;
  int LCT, MF, MT, MF4, MF5, MF6;
  int events;
  int mZA;
  std::string elementsymb;
  
 int mt;
 //std::string secParticle;
 double sigmainElastic=0.0;
 double En;
 double costhlab;
 double secE; 
 double secpcosAng; 
 double sigma; 
  double Ein;
  std::vector<double> xec;
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
  ClassDef(TNudyElastic, 1) // class for TNudyElastic
#endif
};
#endif
