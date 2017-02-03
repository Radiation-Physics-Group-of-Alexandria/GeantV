#ifndef __TNudyEndfThermal__
#define __TNudyEndfThermal__

#include <vector>
#include <fstream>
class TNudyEndfFile;
class TNudyEndfList;
class TList;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

#define PI acos(-1.0)
#define x2(x) (x * x)
#define x3(x) (x * x2(x))
#define x4(x) (x2(x) * x2(x))
#define x5(x) (x * x4(x))
#define x6(x) (x4(x) * x2(x))
#define x8(x) (x6(x) * x2(x))
#define kconst 2.196771e-3
#define Mn 939.565378e+6       // in eV/c^2   1.00866491588   // in u
#define hcross 6.5821192815e-6 // eV.s  1.054571800e-34 // in J.s
#define Fac1(x) (1.0 + x2(x))
#define Fac2(x) (9.0 + 3.0 * x2(x) + x4(x))
#define Fac3(x) 225.0 + 50445.0 * x2(x) - 13500.0 * x3(x) + 1386.0 * x4(x) - 60.0 * x5(x) + x6(x)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowint> matrixint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;
typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
typedef std::vector<std::vector<std::vector<std::vector<rowd>>>> matrixd5;


class TNudyEndfThermal {

public:
  TNudyEndfThermal();
  TNudyEndfThermal(const char *irENDF, double isigDiff);
  virtual ~TNudyEndfThermal();
  void GetData(const char *irENDF, double isigDiff);
  //double nThermalElasticXsecion();
  double nThermalElasticXsecion(double inputEnergy, double inputTemprature);
  void   storeElasticXsecion();
  
  private:
 void ReadFile7(TNudyEndfFile *file);
 double LinearInterpolFile7(double x1, double x2, double sig1, double sig2, rowd x3, rowd x4);
 
  double fun(double, double, double, double);
  double simpson(double, double, int);
  
  
  const char *rENDF;             // precision for cross-section reconstruction
  double sigDiff;                // precision for cross-section reconstruction
  rowint nbt1, int1;
  rowint nbt2, int2;
  rowint nbt3, int3;
  rowd eLinearFile3;
  rowd xLinearFile3;
  int NR, NP, NE; // standard ENDF parameters for range and interpolation
  int flagRead = -1;
  rowint MtLct,Nbeta;                // LCT numbers
  rowint MtNumbers, MtNum4, MtNum5, MtNum6; // MT numbers
  rowint mtLTHR;
  
  
  rowd eintFile1, nutFile1, einFile1, nuFile1;
  rowd eindFile1, nudFile1, einphFile1, phFile1;
  rowd einfFile1, heatFile1;
  
  rowd eLinElastic;
  rowd sLinElastic;
  
   
  
  rowd dummy_eElastic;
  rowd dummy_sElastic;
  rowd dummy2_eElastic;
  rowd dummy2_sElastic;
  
  rowd eIncInelastic;
  rowd sIncInelastic;
  
  rowd tempsIncInelastic;
  
  
  matrixd2 alpint;
  matrixd2 sint;
  
  
  matrixd2 tempSet;//S(E,T) for each values of T
  
  matrixd2 incIneS;
  
  matrixd2 eneUni,nThermalElasticXsec;
  rowd tempThermalXsec;
  rowd tempEthermal;
  double beta[200];
  //static const int n1=5000;
  int nEnergygrid =0 ;
  
  double tSet[1000][1000];
  double sigb;   // bound x-section
  double nXsec=0; // coherent elastic x-section
  double sig_ies = 0; // coherent inelastic x-section
 // double DW;    // DebyeWaller factor
  rowd temp;      // temperature for incoherent elastic scattering
  rowd DW; // DebyeWaller factor
  //double En_temp; // temporary energy of neutron
  //TRandom3 r1;//=new TRandom3();
  rowd BN;
  rowd temp_T;
  matrixd3 tempalpS;// temperature,alpha,S
  rowd cnc, Set;
  matrixd2 tempIncoElastic, DWIncoElastic;                     // total incident energy and nu,  all elements
  matrixd2 einp, nup;                     // prompt incident energy and nu,  all elements
  matrixd2 eind, nud, lambdaD;            // delayed incident energy and nu,  all elements
  matrixd2 einFissHeat, fissHeat;         // fission incident energy and heat,  all elements
  matrixd2 einfId, qvalue;                // incident energy for fission yield
  matrixd3 zafId, pdfYieldId, cdfYieldId; // za and yield fission
  rowd ein;
  TRandom3 *r1;
  TRandom3 *r2;
     // dimensions
  int N = 10000;
  int M = 10000;
   
  // dynamic allocation
  double** E; 
  double** X;  
  double* tempE; 
  double* tempX; 
  double* Temprature;
  
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyEndfThermal, 1) // class for an ENDF reconstruction
#endif
};
#endif
