#include "TEXsec.h"
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TText.h>
#include <TFrame.h>
#include <TGFrame.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#endif
#endif
#include "base/Global.h"
#ifndef VECCORE_CUDA
#include "base/MessageLogger.h"
#endif
using vecgeom::kAvogadro;
using std::numeric_limits;
using std::string;

#ifndef VECCORE_CUDA

using std::max;

#else

template<class T>
VECCORE_ATT_HOST_DEVICE
const T& max(const T&a,const T& b) {
    return (a<b)?b:a;
}
#endif

#ifndef VECCORE_CUDA
int TEXsec::fNLdElems = 0;
TEXsec *TEXsec::fElements[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#else
int fNLdElemsDev = 0;            //! number of loaded elements
TEXsec *fElementsDev[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int fNLdElemsHost = 0;            //! number of loaded elements
TEXsec *fElementsHost[NELEM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#endif

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
TGMainFrame *TEXsec::fMain = 0;
TGHorizontalFrame *TEXsec::fSecond = 0;
TRootEmbeddedCanvas *TEXsec::fCanvas = 0;
TGListBox *TEXsec::fReactionBox = 0;
TGListBox *TEXsec::fParticleBox = 0;
#endif
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
TEXsec::TEXsec()
    : fEGrid(TPartIndex::I()->EGrid()), fAtcm3(0), fEmin(0), fEmax(0), fEilDelta(0),
      fEle(0), fIndex(-1), fNEbins(0), fNRpart(0), fPXsec(nullptr), fPXsecP(nullptr) {
  fName[0] = '\0';
  fTitle[0] = '\0';
}

//___________________________________________________________________
TEXsec::TEXsec(int z, int a, float dens, int np)
    : fEGrid(TPartIndex::I()->EGrid()), fAtcm3(dens * kAvogadro * 1e-24 / TPartIndex::I()->WEle(z)),
      fEmin(TPartIndex::I()->Emin()), fEmax(TPartIndex::I()->Emax()),
      fEilDelta(TPartIndex::I()->EilDelta()), fEle(z * 10000 + a * 10), fIndex(-1),
      fNEbins(TPartIndex::I()->NEbins()), fNRpart(np),
      fPXsec(new TPXsec[fNRpart]),
      fPXsecP(new TPXsec*[fNRpart])
{
  strncpy(fName, TPartIndex::I()->EleSymb(z), 31);
  strncpy(fTitle, TPartIndex::I()->EleName(z), 127);
  memset(fCuts, 0, 4 * sizeof(float));
  for(auto i=0; i<fNRpart; ++i) fPXsecP[i] = &fPXsec[i];
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
TEXsec::TEXsec(const TEXsec &other): fEGrid(other.fEGrid), fAtcm3(other.fAtcm3),
				     fEmin(other.fEmin), fEmax(other.fEmax), fEilDelta(other.fEilDelta),
				     fEle(other.fEle), fIndex(other.fIndex), fNEbins(other.fNEbins),
				     fNRpart(other.fNRpart), fPXsecP(other.fPXsecP)
 {
   #ifndef VECCORE_CUDA
   strncpy(fName,other.fName,32);
   strncpy(fTitle,other.fTitle,128);
   #else
   this->strncpy(fName,other.fName,32);
   this->strncpy(fTitle,other.fTitle,128);
   #endif
   memcpy(fCuts,other.fCuts,4*sizeof(float));
 }

//___________________________________________________________________
TEXsec::~TEXsec() {
   if(fPXsecP)
      for(auto ipart=0; ipart<fNRpart; ++ipart) delete fPXsecP[ipart];
   delete [] fPXsecP;
   delete [] fPXsec;
}

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
//______________________________________________________________________________
void TEXsec::Streamer(TBuffer &R__b)
{
   // Stream an object of class TEXsec.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TEXsec::Class(),this);

      if(fPXsecP != nullptr)
	 for(auto ipart=0; ipart<fNRpart; ++ipart) delete fPXsecP[ipart];
      delete [] fPXsecP;
      fPXsecP = new TPXsec*[fNRpart];
      for(auto i=0; i<fNRpart; ++i) fPXsecP[i] = &fPXsec[i];

   } else {
      R__b.WriteClassBuffer(TEXsec::Class(),this);
   }
}
#endif
#endif

//___________________________________________________________________
bool TEXsec::AddPart(int kpart, int pdg, int nxsec) {
   if(fPXsecP[kpart] == nullptr) fPXsecP[kpart] = new TPXsec();
   return fPXsecP[kpart]->SetPart(pdg, nxsec); }

//___________________________________________________________________
bool TEXsec::AddPartMS(int kpart, const float angle[], const float ansig[], const float length[],
                       const float lensig[]) {
   if(fPXsecP[kpart] == nullptr) fPXsecP[kpart] = new TPXsec();
   return fPXsecP[kpart]->SetPartMS(angle, ansig, length, lensig);
}

//___________________________________________________________________
bool TEXsec::AddPartXS(int kpart, const float xsec[], const int dict[]) {
   if(fPXsecP[kpart] == nullptr) fPXsecP[kpart] = new TPXsec();
   return fPXsecP[kpart]->SetPartXS(xsec, dict); }

//___________________________________________________________________
bool TEXsec::AddPartIon(int kpart, const float dedx[]) {
   if(fPXsecP[kpart] == nullptr) fPXsecP[kpart] = new TPXsec();
   return fPXsecP[kpart]->SetPartIon(dedx); }

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
float TEXsec::XS(int pindex, int rindex, float en) const {
   return fPXsecP[pindex]->XS(rindex, en); }

//___________________________________________________________________
float TEXsec::DEdx(int pindex, float en) const {
   return fPXsecP[pindex]->DEdx(en); }

//___________________________________________________________________
bool TEXsec::MS(int pindex, float en, float &ang, float &asig, float &len, float &lsig) const {
  return fPXsecP[pindex]->MS(en, ang, asig, len, lsig);
}

//___________________________________________________________________
void TEXsec::DumpPointers() const {
  Geant::Printf("Material %d emin %f emax %f NEbins %d EilDelta %f Npart %d\n", fEle, fEmin, fEmax, fNEbins, fEilDelta,
                fNRpart);
  for (int i = 0; i < fNRpart; ++i)
    fPXsecP[i]->Dump();
}

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
//___________________________________________________________________
TGraph *TEXsec::XSGraph(const char *part, const char *reac, float emin, float emax, int nbin) const {
  char title[200];
  const double delta = exp(log(emax / emin) / (nbin - 1));
  double en = emin;
  int pindex = TPartIndex::I()->PartIndex(part);
  if (pindex < 0) {
    Error("XSGraph", "Unknown particle %s\n", part);
    return 0;
  }
  float *xsec = new float[nbin];
  float *energy = new float[nbin];
  int proc = TPartIndex::I()->ProcIndex(reac);
  for (int i = 0; i < nbin; ++i) {
    energy[i] = en;
    xsec[i] = XS(pindex, proc, en);
    en *= delta;
  }
  TGraph *tg = new TGraph(nbin, energy, xsec);
  memset(title, 0, 200);
  snprintf(title, 199, "%s %s on %s", part, reac, GetTitle());
  tg->SetTitle(title);
  delete[] xsec;
  delete[] energy;
  return tg;
}

//___________________________________________________________________
TGraph *TEXsec::DEdxGraph(const char *part, float emin, float emax, int nbin) const {
  char title[200];
  const double delta = exp(log(emax / emin) / (nbin - 1));
  float *dedx = new float[nbin];
  float *energy = new float[nbin];
  double en = emin;
  int pindex = TPartIndex::I()->PartIndex(part);
  if (pindex < 0) {
    Error("DEdxGraph", "Unknown particle %s\n", part);
    delete[] dedx;
    delete[] energy;
    return 0;
  }
  for (int i = 0; i < nbin; ++i) {
    energy[i] = en;
    dedx[i] = DEdx(pindex, en);
    en *= delta;
  }
  TGraph *tg = new TGraph(nbin, energy, dedx);
  memset(title, 0, 200);
  snprintf(title, 199, "%s dEdx on %s", part, GetTitle());
  tg->SetTitle(title);
  delete[] dedx;
  delete[] energy;
  return tg;
}

//___________________________________________________________________
TGraph *TEXsec::MSGraph(const char *part, const char *what, float emin, float emax, int nbin) const {
  char title[200];
  const char *whatname[4] = {"MSangle", "MSangle_sig", "MSCorr", "MSCorr_sig"};
  const double delta = exp(log(emax / emin) / (nbin - 1));
  double en = emin;
  int pindex = TPartIndex::I()->PartIndex(part);
  if (pindex < 0) {
    Error("MSGraph", "Unknown particle %s\n", part);
    return 0;
  }
  int iopt = 4;
  while (iopt--)
    if (!strcmp(whatname[iopt], what))
      break;
  if (iopt < 0) {
    Error("MSGraph", "Unknown parameter %s\nShould be one of %s %s %s %s\n", what, whatname[0], whatname[1],
          whatname[2], whatname[3]);
    return 0;
  }
  float *mscat = new float[nbin];
  float *energy = new float[nbin];
  for (int i = 0; i < nbin; ++i) {
    float answ[4];
    energy[i] = en;
    MS(pindex, en, answ[0], answ[1], answ[2], answ[3]);
    mscat[i] = answ[iopt];
    en *= delta;
  }
  TGraph *tg = new TGraph(nbin, energy, mscat);
  memset(title, 0, 200);
  snprintf(title, 199, "%s %s on %s", part, whatname[iopt], GetTitle());
  tg->SetTitle(title);
  delete[] mscat;
  delete[] energy;
  return tg;
}

//___________________________________________________________________
void TEXsec::Draw(const char *option) // mode=0->terminal, mode=1->viewer
{
  // Draw cross sections and other physics quantities for this material
  //
  string help("Method to draw the cross sections. The format is:\n");
  help += "particle,reaction|reaction|reaction,emin,emax,nbin:<options for draw>\n";
  help += "Available reactions are:\n";
  help += "Total,Transport,MultScatt,Ionisation,Decay,inElastic,Elastic,Capture,Brehms\n";
  help += "PairProd,dEdx,MSang,MSang_sig,MSCorr,MSCorr_sig\n";
  help += "the option All can be given to draw all reaction cross sections for a given particle\n";
  help += "the option \"same\" superimposes the plot on the previous one";

  // CoulombScatt,Photoel,Compton,Conversion,Capture,Killer,dEdx,MSangle,MSlength
  const EColor col[14] = {kBlack,  kGray,   kRed,  kGreen, kBlue,   kMagenta, kCyan,
                          kOrange, kSpring, kTeal, kAzure, kViolet, kPink};
  char title[200];
  char gtitle[200];

  TString opt = option;

  static int isame = 0;

  TObjArray *sections = opt.Tokenize(":");
  TString sec1 = ((TObjString *)sections->At(0))->GetName();
  TString sec2;
  if (sections->GetEntries() > 1)
    sec2 = ((TObjString *)sections->At(1))->GetName();
  bool same = sec2.Contains("same");

  TObjArray *token = sec1.Tokenize(",");
  int narg = token->GetEntries();
  if (narg < 2) {
    Info("Draw", "%s", help.c_str());
    return;
  }
  const char *part = ((TObjString *)token->At(0))->GetName();
  TString reactions = ((TObjString *)token->At(1))->GetName();

  float emin = TPartIndex::I()->Emin();
  if (narg > 2)
    sscanf(((TObjString *)token->At(2))->GetName(), "%f", &emin);
  float emax = TPartIndex::I()->Emax();
  if (narg > 3)
    sscanf(((TObjString *)token->At(3))->GetName(), "%f", &emax);
  int nbin = 100;
  if (narg > 4)
    sscanf(((TObjString *)token->At(4))->GetName(), "%d", &nbin);
  int mode = 0;
  if (narg > 5)
    sscanf(((TObjString *)token->At(5))->GetName(), "%d", &mode);
  if (gFile)
    gFile->Get("PartIndex");
  int pindex = TPartIndex::I()->PartIndex(part);
  if (pindex < 0) {
    Error("Draw", "Unknown particle %s\n", part);
    return;
  } else if (pindex > TPartIndex::I()->NPartReac()) {
    Error("Draw", "No reaction for particle %s\n", part);
    return;
  }
  if (reactions.Contains("All") || reactions.Contains("*")) {
    TString allrea = "";
    for (int i = 0; i < TPartIndex::I()->NProc() - 1; ++i) {
      if (XS(pindex, i, emin) >= 0)
        allrea = allrea + TPartIndex::I()->ProcName(i) + "|";
    }
    allrea += TPartIndex::I()->ProcName(TPartIndex::I()->NProc() - 1);
    reactions.ReplaceAll("All", allrea);
    reactions.ReplaceAll("*", allrea);
  }
  TObjArray *rnames = reactions.Tokenize("|");
  int nreac = rnames->GetEntries();
  snprintf(gtitle, 199, "%s %s on %s", part, reactions.ReplaceAll("|", ",").Data(), GetTitle());
  TMultiGraph *tmg = new TMultiGraph("GV", gtitle);
  TLine **line = new TLine *[nreac];
  TText **text = new TText *[nreac];
  float lstartx = 0.7;

  TCanvas *tc;
  if (mode) {
    tc = fCanvas->GetCanvas();
    tc->SetLogx();
    tc->SetLogy();
    tc->SetGrid();
    same = false;
  } else {
    tc = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("GVcanvas");
    //   if(!tc) tc = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
    if (!tc) {
      tc = new TCanvas("GVcanvas", gtitle, 600, 400);
      tc->SetLogx();
      tc->SetLogy();
      tc->SetGrid();
      same = false;
      sec2.ReplaceAll("same", "");
    }
  }

  if (same) {
    sec2.ReplaceAll("same", "C");
    ++isame;
  } else {
    isame = 0;
    tc->Clear();
    tc->SetTitle(gtitle);
  }

  if (isame == 2 || isame == 3)
    lstartx = 0.3;
  const float lendx = lstartx + 0.05;
  const float tstartx = lstartx + 0.07;
  float lstarty = 0.85;
  if (isame == 1 || isame == 3)
    lstarty = 0.4;
  const float lstepy = 0.03;
  int coff = isame;
  char ytitle[50];
  double ymin = 1e10;
  for (int j = 0; j < nreac; ++j) {
    const char *reac = ((TObjString *)rnames->At(j))->GetName();
    TGraph *tg;
    if (TString(reac).BeginsWith("MS")) {
      tg = MSGraph(part, reac, emin, emax, nbin);
      for (int i = 0; i < tg->GetN(); ++i) {
        double y = tg->GetY()[i];
        if (y > 0 && y < ymin)
          ymin = y;
      }
      snprintf(title, 199, "%s %s on %s", part, reac, GetTitle());
      tg->SetName(title);
      tg->SetTitle(title);
      tg->SetLineColor(col[(coff + j) % 14]);
      tmg->Add(tg, sec2.Data());
      line[j] = new TLine(lstartx, lstarty - lstepy * j, lendx, lstarty - lstepy * j);
      line[j]->SetLineColor(col[(coff + j) % 14]);
      text[j] = new TText(tstartx, lstarty - lstepy * j, reac);

      if (TString(reac).BeginsWith("MSangle"))
        snprintf(ytitle, 49, "Radians");
      else
        snprintf(ytitle, 49, "Relative Step Correction");
    } else if (!strcmp(reac, "dEdx")) {
      snprintf(ytitle, 49, "GeV/cm / g/cm^3");
      tg = DEdxGraph(part, emin, emax, nbin);
      for (int i = 0; i < tg->GetN(); ++i) {
        double y = tg->GetY()[i];
        if (y > 0 && y < ymin)
          ymin = y;
      }
      snprintf(title, 199, "%s dEdx on %s", part, GetTitle());
      tg->SetName(title);
      tg->SetTitle(title);
      tg->SetLineColor(col[(coff + j) % 14]);
      tmg->Add(tg, sec2.Data());
      line[j] = new TLine(lstartx, lstarty - lstepy * j, lendx, lstarty - lstepy * j);
      line[j]->SetLineColor(col[(coff + j) % 14]);
      text[j] = new TText(tstartx, lstarty - lstepy * j, reac);
    } else {
      snprintf(ytitle, 49, "barn");
      if (TPartIndex::I()->ProcIndex(reac) < 0) {
        Error("Draw", "Reaction %s does not exist\n", reac);
        TPartIndex::I()->Print("reactions");
        Geant::Printf("dEdx, MSangle, MSangle_sig, MSCorr, MSCorr_sig\n");
        delete[] text;
        delete[] line;
        return;
      }
      tg = XSGraph(part, reac, emin, emax, nbin);
      for (int i = 0; i < tg->GetN(); ++i) {
        double y = tg->GetY()[i];
        if (y > 0 && y < ymin)
          ymin = y;
      }
      // a x-sec less than 1nb makes little sense...
      ymin = max<double>(ymin, 2e-9);
      snprintf(title, 199, "%s %s on %s", part, reac, GetTitle());
      tg->SetName(title);
      tg->SetTitle(title);
      tg->SetLineColor(col[(coff + j) % 14]);
      tmg->Add(tg, sec2.Data());
      line[j] = new TLine(lstartx, lstarty - lstepy * j, lendx, lstarty - lstepy * j);
      line[j]->SetLineColor(col[(coff + j) % 14]);
      text[j] = new TText(tstartx, lstarty - lstepy * j, reac);
    }
  }
  const char *gopt = 0;
  if (strlen(sec2.Data()))
    gopt = sec2.Data();
  else
    gopt = "AC";
  tmg->SetMinimum(0.5 * ymin);
  tmg->Draw(gopt);
  if (!same) {
    ((TAxis *)tmg->GetHistogram()->GetXaxis())->SetTitle("Energy (GeV)");
    ((TAxis *)tmg->GetHistogram()->GetYaxis())->SetTitle(ytitle);
    ((TAxis *)tmg->GetHistogram()->GetXaxis())->CenterTitle();
    ((TAxis *)tmg->GetHistogram()->GetYaxis())->CenterTitle();
  }
  TText **ptext = new TText *[nreac];
  char string[100] = {"\0"};
  if (same) {
    snprintf(string, 99, "%s on %s", part, GetTitle());
  }
  for (int j = 0; j < nreac; ++j) {
    line[j]->SetBit(TLine::kLineNDC);
    line[j]->Draw();
    line[j]->SetLineColor(col[(coff + j) % 14]);
    text[j]->SetNDC(true);
    text[j]->SetTextSize(0.03);
    text[j]->SetTextAlign(12);
    text[j]->Draw();
    if (same) {
      ptext[j] = new TText(lstartx * 0.95, lstarty - lstepy * j, string);
      ptext[j]->SetTextSize(0.03);
      ptext[j]->SetNDC(true);
      ptext[j]->SetTextAlign(32);
      ptext[j]->Draw();
    } else
      ptext[j] = 0;
  }
  tc->Update();
  for (int i = 0; i < nreac; ++i) {
    delete ptext[i];
    delete text[i];
    delete line[i];
  }
  delete[] ptext;
  delete[] text;
  delete[] line;
}

//___________________________________________________________________

//  All the methods below are used by the Viewer
void TEXsec::Viewer() {
  if (fMain == NULL) { // if the viewer ins't already active
    fMain = new TGMainFrame(gClient->GetRoot(), 1050, 520);
    fMain->Connect("CloseWindow()", "TEXsec", this, "ResetFrame()");
    fSecond = new TGHorizontalFrame(fMain, 200, 40);         // contains everything
    fCanvas = new TRootEmbeddedCanvas(0, fSecond, 900, 520); // contains the graphs

    TGVerticalFrame *fVertical =
        new TGVerticalFrame(fSecond, 200, 520); // Create a vertical frame containing the lists & the buttons
    TGHorizontalFrame *fSelect =
        new TGHorizontalFrame(fVertical, 200, 40); // contains the select all & deselect all buttons

    TGLabel *lParticle = new TGLabel(fVertical, "Particles : "); // label for the particle selector
    lParticle->SetTextFont("verdana");
    fVertical->AddFrame(lParticle, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));

    fParticleBox = new TGListBox(fVertical);
    fParticleBox->Connect("Selected(int)", "TEXsec", this,
                          "UpdateReactions()"); // when the user changes the particle with the mouse

    for (int i = 0; i < TPartIndex::I()->NPartReac(); i++) {
      fParticleBox->AddEntry(TPartIndex::I()->PartName(i), i); // adding all the particles of the list
    }

    fParticleBox->SortByName();
    fParticleBox->Select(8);
    fParticleBox->Resize(150, 160);
    fVertical->AddFrame(fParticleBox, new TGLayoutHints(kLHintsCenterX, 5, 5, 5, 5));

    TGLabel *lSpace = new TGLabel(fVertical, "         ");
    fVertical->AddFrame(lSpace, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));

    TGLabel *lReaction = new TGLabel(fVertical, "Reactions : "); // label for the reaction selector
    lReaction->SetTextFont("verdana");
    fVertical->AddFrame(lReaction, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));

    fReactionBox = new TGListBox(fVertical, 90);
    fReactionBox->Connect("Selected(int)", "TEXsec", this, "PreDraw()");
    this->UpdateReactions();
    fReactionBox->Resize(150, 160);
    fReactionBox->SetMultipleSelections(true);
    fVertical->AddFrame(fReactionBox, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));

    TGTextButton *select = new TGTextButton(fSelect, "&Select All"); // adding the 3 buttons
    select->Connect("Clicked()", "TEXsec", this, "SelectAll()");
    fSelect->AddFrame(select, new TGLayoutHints(kLHintsLeft, 0, 5, 3, 4));
    TGTextButton *deselect = new TGTextButton(fSelect, "Deselect &All");
    deselect->Connect("Clicked()", "TEXsec", this, "DeselectAll()");
    fSelect->AddFrame(deselect, new TGLayoutHints(kLHintsRight, 5, 5, 3, 4));
    fVertical->AddFrame(fSelect, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));

    TGTextButton *exit = new TGTextButton(fVertical, "&Exit", "gApplication->Terminate(0)");
    fVertical->AddFrame(exit, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));

    fSecond->AddFrame(fVertical, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
    fSecond->AddFrame(fCanvas, new TGLayoutHints(kLHintsRight | kLHintsExpandY | kLHintsExpandX, 0, 0, 1, 1));
    fMain->AddFrame(fSecond, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 0, 0, 1, 1));

    string name = "ROOT XSec Viewer - ";
    name += GetTitle();
    fMain->SetWindowName(name.c_str());
    fMain->MapSubwindows();
    fMain->MapWindow();
    fMain->Layout();
  }
}
//___________________________________________________________________
void TEXsec::SelectAll() {

  for (int i = 0; i < fReactionBox->GetNumberOfEntries(); i++) {
    if (!fReactionBox->GetEntry(i)->IsActive()) {
      fReactionBox->GetEntry(i)->Toggle();
    }
  }
  fReactionBox->Layout();
  PreDraw();
}

//___________________________________________________________________
void TEXsec::DeselectAll() {
  for (int i = 0; i < fReactionBox->GetNumberOfEntries(); i++) {
    if (fReactionBox->GetEntry(i)->IsActive()) {
      fReactionBox->GetEntry(i)->Toggle();
    }
  }
  fReactionBox->Layout();
  PreDraw();
}

//___________________________________________________________________
void TEXsec::UpdateReactions() // when the user changes the selected particle
{
  fReactionBox->RemoveAll();
  if (fParticleBox->GetSelected() > -1) { // if a particle is selected
    int nb = 0;
    for (int i = 0; i < TPartIndex::I()->NProc() - 1; ++i) {
      if (XS(TPartIndex::I()->PartIndex(((TGTextLBEntry *)fParticleBox->GetSelectedEntry())->GetTitle()), i,
             0.0000000001) >= 0) {
        fReactionBox->AddEntry(TPartIndex::I()->ProcName(i),
                               nb); // adding the possible reactions for the chosen particle
        fReactionBox->GetEntry(nb)->Toggle();
        nb++;
      }
    }
  }

  fReactionBox->Layout();

  PreDraw();
}

//___________________________________________________________________
void TEXsec::PreDraw() // preparation of Draw() in the viewer, generation of the option
{
  if (fParticleBox->GetSelected() > -1) { // if a particle is selected
    TString opt = ((TGTextLBEntry *)fParticleBox->GetSelectedEntry())->GetTitle();

    opt += ",Total|";
    for (int i = 0; i < fReactionBox->GetNumberOfEntries(); i++) {
      if (fReactionBox->GetEntry(i)->IsActive()) {
        opt += ((TGTextLBEntry *)fReactionBox->GetEntry(i))->GetTitle();
        opt += "|";
      }
    }
    opt += ",0.000000001,1000000000,100,1";
    Draw(opt);
  }
}

//___________________________________________________________________
void TEXsec::ResetFrame() {
  delete fMain;
  fMain = NULL;
}
#endif
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TEXsec::SizeOf() const {
   size_t size = sizeof(*this);
   for(auto i=0; i<fNRpart; ++i)
      size += fPXsecP[i]->SizeOf();
  size -= sizeof(char); // fStore already holds one char
  size = sizeof(double)*((size-1)/sizeof(double)+1);
  return (int) size;
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TEXsec::Compact() {
   char *start = fStore;
   for(auto i=0; i<fNRpart; ++i) {
      TPXsec *px = new(start) TPXsec(*fPXsecP[i]);
      px->Compact();
      // This line can be reactivated when we remove the back compat array
      //      fPXsecP[i]->~TPXsec();
      start += px->SizeOf();
      fPXsecP[i]=px;
   }
}

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TEXsec::RebuildClass() {
  if(((unsigned long) this) % sizeof(double) != 0) {
    Geant::Fatal("aTEXsec::RebuildClass","%s","the class is misaligned\n");
    return;
  }
   char *start = fStore;
   // We consider that the fPXsecP pointer is bad and we do not delete it
   // this is the case if we read back the array.
   fPXsecP = new TPXsec*[fNRpart];
   for(auto i=0; i<fNRpart; ++i) {
#ifdef MAGIC_DEBUG
      if(((TPXsec*) start)->GetMagic() != -777777) {
        Geant::Fatal("TEXsec::RebuildClass","Broken magic %d\n",((TPXsec*) start)->GetMagic());
        return;
      }
#endif
      ((TPXsec *) start)->RebuildClass();
      fPXsecP[i] = (TPXsec *) start;
      if(!fPXsecP[i]->CheckAlign()) 
#ifndef VECCORE_CUDA
      exit(1);;
#else
      return;
#endif 
      start += ((TPXsec*) start)->SizeOf();
   }
}

//___________________________________________________________________
#ifndef VECCORE_CUDA
int TEXsec::SizeOfStore() {
   // First calculate how much we need
   int totsize = 0;
   for(auto i=0; i<fNLdElems; ++i) totsize += fElements[i]->SizeOf();
   return totsize + 2 * sizeof(int);
}
#else
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
VECCORE_ATT_DEVICE
int TEXsec::SizeOfStore() {
   // First calculate how much we need
   int totsize = 0;
   for(auto i=0; i<fNLdElemsDev; ++i) totsize += fElementsDev[i]->SizeOf();
   return totsize + 2 * sizeof(int);
}
#else 
int TEXsec::SizeOfStore() {
   // First calculate how much we need
   int totsize = 0;
   for(auto i=0; i<fNLdElemsHost; ++i) totsize += fElementsHost[i]->SizeOf();
   return totsize + 2 * sizeof(int);
}
#endif
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
size_t TEXsec::MakeCompactBuffer(char* &b) {
   // First calculate how much we need
   size_t totsize = SizeOfStore();
   if(b == nullptr) {
   // Now allocate buffer
     b = (char*) malloc(totsize);
     memset(b,0,totsize);
   }
   char* start = b;
   memcpy(start,&totsize,sizeof(int));
   start += sizeof(int);
  #ifndef VECCORE_CUDA
   memcpy(start,&fNLdElems,sizeof(int));
   start += sizeof(int);
   // now copy and compact
   for(auto i=0; i<fNLdElems; ++i) {
      TEXsec *el = new(start) TEXsec(*fElements[i]);
      el->Compact();
      start += el->SizeOf();
   }
   #else
  #ifdef VECCORE_CUDA_DEVICE_COMPILATION
   memcpy(start,&fNLdElemsDev,sizeof(int));
   start += sizeof(int);
   // now copy and compact
   for(auto i=0; i<fNLdElemsDev; ++i) {
      TEXsec *el = new(start) TEXsec(*fElementsDev[i]);
      el->Compact();
      start += el->SizeOf();
   }
   #else
   memcpy(start,&fNLdElemsHost,sizeof(int));
   start += sizeof(int);
   // now copy and compact
   for(auto i=0; i<fNLdElemsHost; ++i) {
      TEXsec *el = new(start) TEXsec(*fElementsHost[i]);
      el->Compact();
      start += el->SizeOf();
   }
   #endif
   #endif
   return totsize;
}

#ifndef VECCORE_CUDA
#define EXIT_OR_RETURN(val) exit(val);
#else
#define EXIT_OR_RETURN(val) return;
#endif

//___________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TEXsec::RebuildStore(char *b) {
   size_t size = 0;
   typedef TEXsec *ElementArray_t[NELEM];
#ifndef VECCORE_CUDA
   int &nElems(fNLdElems);
   ElementArray_t &elements(fElements);
#elif defined(VECCORE_CUDA_DEVICE_COMPILATION)
   int &nElems(fNLdElemsDev);
   ElementArray_t &elements(fElementsDev);
#else // NVCC Host
   int &nElems(fNLdElemsHost);
   ElementArray_t &elements(fElementsHost);
#endif
   char *start = b;
   memcpy(&size,start,sizeof(int));
   start += sizeof(int);
   memcpy(&nElems,start,sizeof(int));
   start += sizeof(int);
   for(auto i=0; i<nElems; ++i) {
      TEXsec *current = (TEXsec *) start;
#ifdef MAGIC_DEBUG
      if(current->GetMagic() != -777777) {
        Geant::Fatal("TEXsec::RebuildStore:","Broken Magic@ %d %d\n", i, current->GetMagic());
        EXIT_OR_RETURN(1);
      }
#endif
      current->RebuildClass();
      elements[i] = current;
      if(!elements[i]->CheckAlign())
        EXIT_OR_RETURN(1);
      start += current->SizeOf();
   }
   if(size_t (start - b) != size) {
     Geant::Fatal("TEXsec::RebuildStore","expected size %ld found size %ld", size, start - b);
     EXIT_OR_RETURN(1);
   }
}

//___________________________________________________________________
float TEXsec::Lambda(int pindex, double en) const {
  double xs = fPXsecP[pindex]->XS(TPartIndex::I()->NProc() - 1, en);
  return xs ? 1. / (fAtcm3 * xs) : numeric_limits<float>::max();
}

//___________________________________________________________________
bool TEXsec::Lambda_v(int npart, const int pindex[], const double en[], double lam[]) const {
  const int itot = TPartIndex::I()->NProc() - 1;
  for (int ip = 0; ip < npart; ++ip) {
    double xs = fPXsecP[pindex[ip]]->XS(itot, en[ip]);
    lam[ip] = xs ? 1. / (fAtcm3 * xs) : numeric_limits<float>::max();
  }
  return true;
}

//___________________________________________________________________
bool TEXsec::Lambda_v(int npart, int pindex, const double en[], double lam[]) const {
  const int itot = TPartIndex::I()->NProc() - 1;
  return fPXsecP[pindex]->XS_v(npart, itot, en, lam);
}

//___________________________________________________________________
int TEXsec::SampleReac(int pindex, double en) const {
   return fPXsecP[pindex]->SampleReac(en); }

//___________________________________________________________________
int TEXsec::SampleReac(int pindex, double en, double randn) const {
   return fPXsecP[pindex]->SampleReac(en, randn); }

//___________________________________________________________________
#if !defined(USE_ROOT) || defined(VECCORE_CUDA)
VECCORE_ATT_HOST_DEVICE
TEXsec *TEXsec::GetElement(int z, int a) {
  int ecode = z * 10000 + a * 10;
#ifdef VECCORE_CUDA
#ifdef VECCORE_CUDA_DEVICE_COMPILATION 
  for (int el = 0; el < fNLdElemsDev; ++el)
    if (ecode == fElementsDev[el]->Ele())
      return fElementsDev[el];
  return 0;
#else
  for (int el = 0; el < fNLdElemsHost; ++el)
    if (ecode == fElementsHost[el]->Ele())
      return fElementsHost[el];
  return 0;
#endif   
#else
  for (int el = 0; el < fNLdElems; ++el)
    if (ecode == fElements[el]->Ele())
      return fElements[el];
  return 0;
#endif   
}
#else
TEXsec *TEXsec::GetElement(int z, int a, TFile *f) {
  //   printf("Getting Element %d %d %d\n",z,a,fNLdElems);
  int ecode = z * 10000 + a * 10;
  for (int el = 0; el < fNLdElems; ++el)
    if (ecode == fElements[el]->Ele())
      return fElements[el];
  // Element not found in memory, getting it from file
  TFile *ff = gFile;
  if (f)
    ff = f;
  if (!ff)
    ::Fatal("TEXsec::GetElement", "No file open!");
  fElements[fNLdElems] = (TEXsec *)ff->Get(TPartIndex::I()->EleSymb(z));
  if (!fElements[fNLdElems]) {
    ::Fatal("GetElement", "Element z %d a %d not found", z, a);
    return 0; // just to make the compiler happy
  } else {
    // We loaded the element, however we have to see whether
    // the energy grid is the right one
    // NO, don't need to check. It will be loaded from xsec.root
    //        if(FloatDiff(TPartIndex::I()->Emin(),fElements[fNLdElems]->Emin(),1e-7) ||
    //           FloatDiff(TPartIndex::I()->Emax(),fElements[fNLdElems]->Emax(),1e-7) ||
    //           TPartIndex::I()->NEbins() != fElements[fNLdElems]->NEbins())
    // we have to resize the energy grid of the element
    //            fElements[fNLdElems]->Resample();
    return fElements[fNLdElems++];
  }
}
#endif
//___________________________________________________________________
bool TEXsec::Prune() {
  for (int ip = 0; ip < fNRpart; ++ip)
    fPXsecP[ip]->Prune();
  return true;
}

//___________________________________________________________________
bool TEXsec::Resample() {
  for (int ip = 0; ip < fNRpart; ++ip)
    fPXsecP[ip]->Resample();
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fNEbins = TPartIndex::I()->NEbins();
  fEGrid = TPartIndex::I()->EGrid();
  return true;
}

//___________________________________________________________________
#ifdef VECCORE_CUDA
VECCORE_ATT_HOST_DEVICE
char *TEXsec::strncpy(char *dest, const char *src, size_t n)
{
    char *ret = dest;
    do {
        if (!n--)
            return ret;
    } while (*dest++ = *src++);
    while (n--)
        *dest++ = 0;
    return ret;
}
#endif

