#ifndef TNudyInelastic_H
#define TNudyInelastic_H
#include "TNudyEndfRecoPoint.h"
#include "TNudyEndfMat.h"
#include <string>
#include "TNudyElement.h"
#define PI acos(-1.0)
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
class TNudyEndfFile;
#endif


class TNudyInelastic {
public:
  TNudyInelastic();
  TNudyInelastic(int ElemId, const char*);
  virtual ~TNudyInelastic();
  void nInelasticXsec(double, TNudyElement * );
  double GetInelasticXsec();
  std::vector<double> GetKiEnInelastic();
  std::vector<double> GetcosAngInelastic();
  std::string GetInelasticProcessName();
  std::vector<std::string>  GetParticleInelastic();
  int reactionChannelNumber();
  std::string residueName(int, int);
  void GetInelasticParameters(double &xsec1, double &E);
  double qValue[1000];
  void Reset();
private:
  void GetSecParameter(TNudyElement *, TNudyEndfRecoPoint *recoPoint);
  double kineticE;
  double cosCM = 0, cosLab = 0, secEnergyCM = 0, secEnergyLab = 0;
  double residueKineticEnergy=0;
  double residueCosang=0;
  double sigmaPartial;
  double sigmaTotal;
  double ratio;
  double sum1;
  double projectileMass = 0, projectileCharge = 0;
  int residueA = 0;
  int residueZ = 0;
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
  int neutronsN;
  static const char fkElNam[119][4];
  int nMTinelastic[200][200];
  int mtValues;
  int nMTs;
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
  int nsecParticles;
  std::vector<double> xec;
  std::string  secParticleInel;
  std::string  procNameInelastic;
  std::string  name;
  std::string residueType;
  std::vector<std::string>ProcName;
  std::vector<double> crs;
  std::vector<double> secEnergyinLab;
  std::vector<double> seccosAnginLab;
  std::vector<std::string> productsName;
  std::vector<std::string> tempproductsName;
  std::vector<double> productMass; //except neutrons
  std::vector<double>residueEnergy;
  std::vector<double>residueAngle;
  std::vector<int> mtTemp;
  std::vector<double> crssIel;
  std::vector<double> crss;
  std::vector<int> zP;
  TNudyEndfEnergyAng *recoEnergyAng;
  TNudyEndfRecoPoint *rp;

#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyInelastic, 1) // class for TNudyInelastic
#endif
};
#endif
