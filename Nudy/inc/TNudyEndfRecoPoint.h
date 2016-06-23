#ifndef __TNudyEndfRecoPoint__
#define __TNudyEndfRecoPoint__

#include <vector>
#include <fstream>
class TNudyEndfNuPh;
class TNudyEndfFissionYield;
class TNudyEndfEnergy;
class TNudyEndfEnergyAng;
class TNudyEndfAng;
class TNudyEndfDoppler;
class TNudyEndfFile;
class TNudyEndfList;
class TList;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom;
#endif

#define PI acos(-1.0)
#define x2(x) (x * x)
#define x3(x) (x * x2(x))
#define x4(x) (x2(x) * x2(x))
#define x5(x) (x * x4(x))
#define x6(x) (x4(x) * x2(x))
#define x8(x) (x6(x) * x2(x))
#define kconst 2.196771e-3
#define Mn 939.565378e+6    // in eV/c^2   1.00866491588   // in u
#define hcross 6.5821192815e-6   // eV.s  1.054571800e-34 // in J.s
#define Fac1(x) (1.0 + x2(x))
#define Fac2(x) (9.0 + 3.0 * x2(x) + x4(x))
#define Fac3(x) 225.0 + 50445.0*x2(x)  - 13500.0*x3(x)  + 1386.0*x4(x)  - 60.0*x5(x)  + x6(x)

class  TNudyEndfRecoPoint {

public: 
  TNudyEndfRecoPoint();
  virtual ~TNudyEndfRecoPoint();
  void GetData(const char *rENDF);
  double SetsigPrecision(double x1){return sigDiff = x1;}
  void broadSigma(std::vector<double> &x1, std::vector<double> &x2, std::vector<double> &x3);
  void fixupTotal(std::vector<double> &x1, std::vector<double> &x2);
  double GetSigmaTotal(double energyK);
  double GetSigmaPartial(int i, double energyK);
  double GetCos4(int mt, double energyK);
  double GetEnergy5(int mt, double energyK);
  std::fstream out,outtotal;
  std::string outstring,outstringTotal;
  std::vector<double> energyUni,sigmaUniTotal;		// unionization of energy and total cross-section
  std::vector<std::vector<double> > sigmaOfMts;         // sigma for each reaction
  std::vector<std::vector<double> > sigmaUniOfMts;      // sigma for each reaction afte unionization of energy
  std::vector<std::vector<int> > MtValues;              // MT values for which cross-section/ heating values are given 
  std::vector<int> energyLocationMts;			// MT wise starting energy for cross-section
  int NoOfElements = 0;
  double sigDiff;					// precision for cross-section reconstruction
protected:
  std::vector<std::vector<std::vector<double> > >cosPdf4OfMts;        // cosine and pdf from file 4 for each reaction
  std::vector<std::vector<std::vector<double> > >cosCdf4OfMts;        // cosine and cdf from file 4 for each reaction
  std::vector<std::vector<double> > energy4OfMts;       // incident energy in file 4 for each reaction
  std::vector<std::vector<int> > Mt4Values;             // MT values for which angular distributions are given in file 4
  std::vector<std::vector<int> > Mt4Lct;                // CM and Lab flag for angular distributions as given in file 4
  std::vector<std::vector<std::vector<double> > >energyPdf5OfMts;        // cosine and pdf from file 4 for each reaction
  std::vector<std::vector<std::vector<double> > >energyCdf5OfMts;        // cosine and cdf from file 4 for each reaction
  std::vector<std::vector<double> > energy5OfMts;       // incident energy in file 4 for each reaction
  std::vector<std::vector<int> > Mt5Values;             // MT values for which angular distributions are given in file 4
  double AWRI;
  double QValue[999];
private:
  void ReadFile2(TNudyEndfFile *file);
  void ReadFile3(TNudyEndfFile *file);
  double recursionLinearFile3(double x1, double x2, double sig1, double sig2, std::vector<double> x3, std::vector<double> x4);
  void ReadResDat4l(int l1, int mm, TNudyEndfList *theList, TList *rec);
  void recoPlusBroad(int flagNer);
  void GetSigma(int lrfp, double x, double &siga, double &sigb, double &sigc);
  void Linearize(int flagNer);  
  void GetSigmaRMP(double x, double &siga, double &sigb, double &sigc);
  void InverseMatrix();
  double backCrsAdler(double x, int l1);
  double calcPhi(double x, int l);
  double calcShift(double x, int l);
  double calcPene(double x, int l);
  double GetERP(double x, int r, int lVal);
  double GetRho(double x, int lVal);
  double GetRhoC(double x, int isDiff, int lVal);
  double Gamma_reduced(double y, int ii, int lval);
  double K_wnum(double x);
  double Gamma_nrE(double x, int ii, int lval);
  double Gamma_xrE(int ii, int lrx);
  double Gamma_rE(double x, int ii, int lval, int lrx);
  void additionalSigma(int LRF, double x);
  double recursionLinear(double x1, double x2, double sig1, double sig2);
  double recursionLinear(double x1, double x2, double sig1, double sig2, double sig3, double sig4, double sig5, double sig6);
  int widthFluctuation(double gnr, double gx, double gg, double gf, int jval);
  double Thinning(std::vector<double> &x1, std::vector<double> &x2);
  double addFile3Resonance(double &x1, double &x2, std::vector<double> &x3, std::vector<double> &x4);
  double insertFile3(std::vector<double> &x1, std::vector<double> &x2);
  double insertFile3High(std::vector<double> &x1, std::vector<double> &x2);

  int Z, ZA, ZAI, LFW, NER, LRU, LRF, NRO, NAPS, NLS, LSSF, NLS2, NJS, INT,NIS,intLinLru1=0; // standard ENDF parameters
  int LRX,cueMat=0;                                     // flag for inelastic reaction, number of J
  int dopplerBroad=0, flagResolve=0, flagUnResolve=0;   // flag for Doppler broadening for thinning, Resolve  and URR parameter exist
  double eLo1=0, eLo2=0, eHi1=0,  eHi2=0, eLo=0, eHi=0; // Resonance energy range Low, High
  double SPI, AP, APL[10], rad_a;                       // Target Spin (I), Scattering Radius (AP), L-dependent AP, Channel radius (a)
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  double factor_k;                                      // factor for wave vector
  double JMIN, JMAX;                                    // J values
  double RN,RG,RF,RX;                                   // standard ENDF parameters
  int totalAdler,crsAdler[4];                           // adler adler cross-section
  double R[3][3],S[3][3];                               // matrix for RM formalism
  double RI[3][3],SI[3][3];                             // matrix for RM formalism
  double MissingJ[5][50], MisGj[5]; int NJValue[5];     // J values in sorted form
  int NR, NP, NE;                         // standard ENDF parameters for range and interpolation
  std::vector<int> MtNumbers;				// MT numbers
  std::vector<double> sigmaMts;				// MT numbers for sigma in file3
  std::vector<double> eLinElastic,eLinCapture,eLinFission;
  std::vector<double> xLinElastic,xLinCapture,xLinFission;
  std::vector<double> xBroadElastic,xBroadCapture,xBroadFission;
  std::vector<int> nbt1,int1;
  std::vector<double> eLinearFile3;
  std::vector<double> xLinearFile3;
  std::vector<double> sigma; 
  std::vector<int> l;					// l values
  std::vector<int> NRS;              			// no. of resolved resonances
  std::vector<int> NRJ;              			// no. of URR J
  std::vector<int> JSM;              			// URR J
  std::vector<double> Er;            			// resolved resonance energy
  std::vector<double> J;            			// associated J
  std::vector<double> GJ;				// spin multiplication factor
  std::vector<double> Gamma_r;       			// total width = Gamma_n + Gamma_g + Gamma_f
  std::vector<double> Gamma_n;       			// neutron scattering width
  std::vector<double> Gamma_g;       			// Capture width
  std::vector<double> Gamma_f;       			// fission width
  std::vector<double> Gamma_x;       			// Inelastic width
  std::vector<double> Gamma_fa,Gamma_fasq;       	// fission width 1
  std::vector<double> Gamma_fb,Gamma_fbsq;       	// fission width 2
  std::vector<double> at1;       			// 1 background constant (Reich-Moore)
  std::vector<double> at2;       			// 2 background constant (Reich-Moore)
  std::vector<double> at3;       			// 3 background constant (Reich-Moore)
  std::vector<double> at4;       			// 4 background constant (Reich-Moore)
  std::vector<double> bt1;       			// 5 background constant (Reich-Moore)
  std::vector<double> bt2;       			// 6 background constant (Reich-Moore)
  std::vector<double> det1;       			// 1 resonance energy (Reich-Moore)
  std::vector<double> dwt1;       			// 2 half width (Reich-Moore)
  std::vector<double> grt1;       			// 3 symmetrical cross-section parameter G (Reich-Moore)
  std::vector<double> git1;       			// 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  std::vector<double> def1;       			// 5 background constant (Reich-Moore)
  std::vector<double> dwf1;       			// 6 background constant (Reich-Moore)
  std::vector<double> grf1;       			// 3 symmetrical cross-section parameter G (Reich-Moore)
  std::vector<double> gif1;       			// 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  std::vector<double> dec1;       			// 5 background constant (Reich-Moore)
  std::vector<double> dwc1;       			// 6 background constant (Reich-Moore)
  std::vector<double> grc1;       			// 3 symmetrical cross-section parameter G (Reich-Moore)
  std::vector<double> gic1;       			// 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  std::vector<double> amux, amun, amug, amuf;		// standard ENDF parameters
  std::vector<double> Es;				// energy URR
  std::vector<double> D, GX, GNO, GG, GF;		// URR parameters
  std::vector<double> PhiEr,ShiftEr;			// penetration and shift factors
  std::vector<double> eneTemp,sigTemp;			// temporary vectors to store energy and sigma
  TNudyEndfDoppler *doppler;
  TNudyEndfAng *recoAng;
  TNudyEndfEnergy *recoEnergy;
  TNudyEndfEnergyAng *recoEnergyAng;
  TNudyEndfNuPh *recoNuPh;
  TNudyEndfFissionYield *recoFissY;
#ifdef USE_ROOT
  TRandom *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyEndfRecoPoint, 1) // class for an ENDF reconstruction
#endif
};
#endif
