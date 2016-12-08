#ifndef TNudyInelastic_H
#define TNudyInelastic_H
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

class TNudyInelastic {
public:
  TNudyInelastic();
  TNudyInelastic(int ElemId, const char*);
  TNudyInelastic(ElementProp *, TNudyEndfRecoPoint *recoPoint);
  virtual ~TNudyInelastic();
  void Neutron_Process(TNudyElement *);
  void nInelasticXsec(double, TNudyElement * );
  double GetInelasticXsec();
  double GetKiEnInelastic();
  double GetcosAngInelastic();
  std::string GetInelasticProcessName();
  std::string  GetParticleInelastic();
  void GetInelasticParameters(double &xsec1, double &E);
private:
  void GetSecParameter(TNudyElement *, TNudyEndfRecoPoint *recoPoint);
  std::vector<double> crs;
  double kineticE;
  double cosCM = 0, cosLab = 0, secEnergyCM = 0, secEnergyLab = 0;
  double cosCMI = 0, cosLabI = 0, secEnergyCMI = 0, secEnergyLabI = 0;
  double sigmaPartial;
  double sigmaTotal;
  double ratio;
  double sum1;
  double residueAIn, residueZIn;
  double residueA=0, residueZ=0;
  double sigmainElastic=0.0;
  double En;
  double costhlab;
  double secE; 
  double secpcosAng; 
  double sigma; 
  double Ein;
  double mass;
  double charge;
  double sigmaInelastic;
  
  int mtValues;
  int nMTs;
  int nMTinelastic[200][200];
  int elemId;
  int pId;
  int ielemId;
  int isel     = 0;
  int counter  = 0;
  int ecounter = 0;
  int LCT, MF, MT, MF4, MF5, MF6;
  int events;
  int mZA;
  int mt,mt1;
 //std::string secParticle;
  std::string elementsymb;
  std::vector<double> xec;
  rowd Eout;
  rowd costheta;
  std::string  secParticleInel;
  std::string  procNameIne;
  std::string  name;
  std::vector<std::string>ProcName;
  std::vector<int> mtTemp;
  std::vector<double> crssIel;
  std::vector<double> crss;
  //TNudyElement *elementProp;
  TNudyEndfRecoPoint *rp;
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyInelastic, 1) // class for TNudyInelastic
#endif
};
#endif
