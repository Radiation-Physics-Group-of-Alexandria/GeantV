// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TEXsec_H
#define TEXsec_H
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEXSec                                                               //
//                                                                      //
// X-section for GV per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPartIndex.h"
#include "TPXsec.h"
#include "Geant/Error.h"

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
#include "Rtypes.h"
class TGHorizontalFrame;
class TGListBox;
class TGMainFrame;
class TGraph;
class TRootEmbeddedCanvas;
#endif
class TFile;
#else
class TEXsec;
extern VECCORE_ATT_DEVICE int fNLdElemsDev;            //! number of loaded elements
extern  int fNLdElemsHost;            //! number of loaded elements
extern VECCORE_ATT_DEVICE TEXsec *fElementsDev[NELEM]; //! databases of elements
extern  TEXsec *fElementsHost[NELEM]; //! databases of elements
#endif

class TEXsec {
public:
  enum { kCutGamma, kCutElectron, kCutPositron, kCutProton };
VECCORE_ATT_HOST_DEVICE
  TEXsec();
  TEXsec(int z, int a, float dens, int np);
VECCORE_ATT_HOST_DEVICE
  TEXsec(const TEXsec &other);
  virtual ~TEXsec();
  static const char *ClassName() { return "TEXsec"; }
  bool AddPart(int kpart, int pdg, int nxsec);
  bool AddPartXS(int kpart, const float xsec[], const int dict[]);
  bool AddPartIon(int kpart, const float dedx[]);
  bool AddPartMS(int kpart, const float angle[], const float ansig[], const float length[], const float lensig[]);

VECCORE_ATT_HOST_DEVICE
  int Ele() const { return fEle; }
VECCORE_ATT_HOST_DEVICE
  int Index() const { return fIndex; }
  void SetIndex(int index) { fIndex = index; }
VECCORE_ATT_HOST_DEVICE
  double Emin() const { return fEmin; }
VECCORE_ATT_HOST_DEVICE
  double Emax() const { return fEmax; }
VECCORE_ATT_HOST_DEVICE
  int NEbins() const { return fNEbins; }
VECCORE_ATT_HOST_DEVICE
  double EilDelta() const { return fEilDelta; }
VECCORE_ATT_HOST_DEVICE
  float XS(int pindex, int rindex, float en) const;
  float DEdx(int pindex, float en) const;
  bool MS(int index, float en, float &ang, float &asig, float &len, float &lsig) const;
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
  TGraph *XSGraph(const char *part, const char *reac, float emin, float emax, int nbin) const;
  TGraph *DEdxGraph(const char *part, float emin, float emax, int nbin) const;
  TGraph *MSGraph(const char *part, const char *what, float emin, float emax, int nbin) const;
#endif
#endif

  float Lambda(int pindex, double en) const;
  bool Lambda_v(int npart, const int pindex[], const double en[], double lam[]) const;
  bool Lambda_v(int npart, int pindex, const double en[], double lam[]) const;
  int SampleReac(int pindex, double en) const;
  int SampleReac(int pindex, double en, double randn) const;

  static bool FloatDiff(double a, double b, double prec) { return fabs(a - b) > 0.5 * fabs(a + b) * prec; }

  const float *Cuts() const { return fCuts; }
  bool SetCuts(const double cuts[4]) {
    for (int jc = 0; jc < 4; ++jc)
      fCuts[jc] = cuts[jc];
    return true;
  }

  void DumpPointers() const;
#ifndef VECCORE_CUDA
  void Draw(const char *option);
  void Viewer(); // *MENU*
  void UpdateReactions();
  void SelectAll();
  void DeselectAll();
  void PreDraw();
  void ResetFrame();
#endif

VECCORE_ATT_HOST_DEVICE
  int SizeOf() const;
VECCORE_ATT_HOST_DEVICE
  void Compact();
VECCORE_ATT_HOST_DEVICE
  void RebuildClass();
VECCORE_ATT_HOST_DEVICE
  static int SizeOfStore();
VECCORE_ATT_HOST_DEVICE
  static size_t MakeCompactBuffer(char* &b);
VECCORE_ATT_HOST_DEVICE
  static void RebuildStore(char *b);
#ifdef MAGIC_DEBUG
VECCORE_ATT_HOST_DEVICE
  int GetMagic() const {return fMagic;}
#endif
#ifdef VECCORE_CUDA
VECCORE_ATT_HOST_DEVICE
char *strncpy(char *dest, const char *src, size_t n);
#endif

VECCORE_ATT_HOST_DEVICE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) fEGrid) % sizeof(fEGrid[0]) != 0) { Geant::Error("TEXsec::CheckAlign","fEGrid","misaligned\n");isaligned=false;}
  if(((unsigned long) &fAtcm3) % sizeof(fAtcm3) != 0) { Geant::Error("TEXsec::CheckAlign","fAtcm3","misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmin) % sizeof(fEmin) != 0) { Geant::Error("TEXsec::CheckAlign","fEmin","misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmax) % sizeof(fEmax) != 0) { Geant::Error("TEXsec::CheckAlign","fEmax","misaligned\n");isaligned=false;}
  if(((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) { Geant::Error("TEXsec::CheckAlign","fEilDelta","misaligned\n");isaligned=false;}
  if(((unsigned long) fCuts) % sizeof(fCuts[0]) != 0) { Geant::Error("TEXsec::CheckAlign","fCuts","misaligned\n");isaligned=false;}
  if(((unsigned long) &fEle) % sizeof(fEle) != 0) { Geant::Error("TEXsec::CheckAlign","fEle","misaligned\n");isaligned=false;}
  if(((unsigned long) &fIndex) % sizeof(fIndex) != 0) { Geant::Error("TEXsec::CheckAlign","fIndex","misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) { Geant::Error("TEXsec::CheckAlign","fNEbins","misaligned\n");isaligned=false;}
  if(((unsigned long) &fNRpart) % sizeof(fNRpart) != 0) { Geant::Error("TEXsec::CheckAlign","fNRpart","misaligned\n");isaligned=false;}
  for(auto i=0; i< fNRpart; ++i)
   if(((unsigned long) fPXsecP[i]) % sizeof(double) != 0) { Geant::Error("TEXsec::CheckAlign","fPXsecP[%d] misaligned\n",i);isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TEXsec::CheckAlign","fMagic","misaligned\n");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TEXsec::CheckAlign","fStore","misaligned\n");isaligned=false;}
  return isaligned;
}

  bool Resample();

  bool Prune();
#ifndef VECCORE_CUDA
  static int NLdElems() { return fNLdElems; }
  static TEXsec *Element(int i) {
    if (i < 0 || i >= fNLdElems)
      return 0;
    return fElements[i];
  }
 #else

#ifdef VECCORE_CUDA_DEVICE_COMPILATION
 VECCORE_ATT_DEVICE 
 static int NLdElems() { return fNLdElemsDev; }
 VECCORE_ATT_DEVICE 
 static TEXsec *Element(int i) {
    if (i < 0 || i >= fNLdElemsDev)
      return 0;
    return fElementsDev[i];
  }
 #else
  static int NLdElems() { return fNLdElemsHost; }
  static TEXsec *Element(int i) {
    if (i < 0 || i >= fNLdElemsHost)
      return 0;
    return fElementsHost[i];
  }
#endif
#endif
  const char *GetName() const { return fName; }
  const char *GetTitle() const { return fTitle; }

#if !defined(USE_ROOT) || defined(VECCORE_CUDA)
  VECCORE_ATT_HOST_DEVICE
  static TEXsec *GetElement(int z, int a = 0);
#else
  static TEXsec *GetElement(int z, int a = 0, TFile *f = 0);
#endif

  
#ifdef VECCORE_CUDA
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  VECCORE_ATT_DEVICE TEXsec **GetElements() { return fElementsDev; }
#else
  TEXsec **GetElements() { return fElementsHost; }
#endif
#else
  static TEXsec **GetElements() { return fElements; }
#endif

private:
#ifndef VECCORE_CUDA
  TEXsec &operator=(const TEXsec &); // Not implemented
#endif
  char fName[32];   // Name
  char fTitle[128]; // Title

  const double *fEGrid; //! Common energy grid
  double fAtcm3;        // Atoms per cubic cm unit density
  double fEmin;         // Minimum of the energy Grid
  double fEmax;         // Maximum of the energy Grid
  double fEilDelta;     // Inverse log energy step
  float fCuts[4];       // Production cuts "a la G4"
  int fEle;             // Element code Z*10000+A*10+metastable level
  int fIndex;           // Index of this in TTabPhysMgr::fElemXsec
  int fNEbins;          // Number of log steps in energy
  int fNRpart;          // Number of particles with reaction

  TPXsec *fPXsec;       // [fNRpart] Cross section table per particle
  TPXsec **fPXsecP;     // [fNRpart] Cross section table per particle
#ifndef VECCORE_CUDA
  static int fNLdElems;            //! number of loaded elements
  static TEXsec *fElements[NELEM]; //! databases of elements
#endif
#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
  static TGMainFrame *fMain;           //! Main window
  static TGHorizontalFrame *fSecond;   //! Window for the graph and the bar on left
  static TRootEmbeddedCanvas *fCanvas; //! For the graphs
  static TGListBox *fReactionBox;      //! Reaction list
  static TGListBox *fParticleBox;      //! Particle list

  ClassDefNV(TEXsec, 7) // Element X-secs
#endif
#endif

private:
  alignas(sizeof(double)) char fStore[1];              //! Pointer to the compact store part of the class
};

#endif
