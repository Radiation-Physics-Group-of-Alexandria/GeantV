// 	This class is reconstructing probability tables for Energy distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyCore.h"
#include "TNudyEndfPhEnergy.h"
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"
using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfPhEnergy)
#include "TRandom3.h"
#endif

    TNudyEndfPhEnergy::TNudyEndfPhEnergy()
{
}

//______________________________________________________________________________
TNudyEndfPhEnergy::TNudyEndfPhEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < sec->GetN1(); k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT              = sec->GetMT();
      int NC              = sec->GetN1();
      MtNumbers.push_back(MT);
      int LF = tab1->GetL2();
      std::cout << " LF = " << LF << " MT " << MT << "  NC " << NC << std::endl;
      NR = tab1->GetN1();
      NP = tab1->GetN2();
      //****************************************************************************
      // arbitrary tabulated function
      if (LF == 1) {
        for (int cr = 0; cr < NR; cr++) {
          nbt1.push_back(tab1->GetNBT(cr));
          int1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < NP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        nr2                 = tab2->GetN1();
        np2                 = tab2->GetN2();

        for (int cr = 0; cr < nr2; cr++) {
          nbt2.push_back(tab2->GetNBT(cr));
          int2.push_back(tab2->GetINT(cr));
        }
        for (int cr = 0; cr < np2; cr++) {
          TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
          ein.push_back(tab12->GetC2());
          nr3 = tab12->GetNR();
          np3 = tab12->GetNP();
          for (int i = 0; i < nr3; i++) {
            nbt3.push_back(tab12->GetNBT(i));
            int3.push_back(tab12->GetINT(i));
          }
          for (int crs = 0; crs < np3; crs++) {
            energyFile5.push_back(tab12->GetX(crs));
            energyPdfFile5.push_back(tab12->GetY(crs));
// 	    std::cout << " E = " << tab12->GetX(crs) <<"  "<< tab12->GetY(crs) << std::endl ;
          }
          for (int cr = 0; cr < np3 - 1; cr++) {
            // std::cout << energyFile5[cr] <<"  "<< energyPdfFile5[cr] << std::endl;
//             recursionLinearFile5Prob(energyFile5[cr], energyFile5[cr + 1], energyPdfFile5[cr], energyPdfFile5[cr + 1]);
          }
          fillPdf1d();
          nbt3.clear();
          int3.clear();
        }
        nbt1.clear();
        int1.clear();
        nbt2.clear();
        int2.clear();
        fE1.clear();
        //****************************************************************************
        // general evaporation spectrum
      } 
    }
    ein2d.push_back(ein);
    frac2d.push_back(fP1);
    ene3d.push_back(ene2d);
    pdf3d.push_back(pdf2d);
    cdf3d.push_back(cdf2d);
    ein.clear();
    ene2d.clear();
    pdf2d.clear();
    cdf2d.clear();
    fP1.clear();
  }
  Mt5Values.push_back(MtNumbers);
  energy5OfMts.push_back(ein2d);
  fraction5OfMts.push_back(frac2d);
  energyOut5OfMts.push_back(ene3d);
  energyPdf5OfMts.push_back(pdf3d);
  energyCdf5OfMts.push_back(cdf3d);
  MtNumbers.clear();
  ein2d.clear();
  ene3d.clear();
  pdf3d.clear();
  cdf3d.clear();
  frac2d.clear();
}

TNudyEndfPhEnergy::~TNudyEndfPhEnergy()
{
  MtNumbers.shrink_to_fit();
  fE1.shrink_to_fit();
  fP1.shrink_to_fit();
  fE2.shrink_to_fit();
  fP2.shrink_to_fit();
  fE3.shrink_to_fit();
  fP3.shrink_to_fit();
  INorm.shrink_to_fit();
  nbt1.shrink_to_fit();
  int1.shrink_to_fit();
  nbt2.shrink_to_fit();
  int2.shrink_to_fit();
  nbt3.shrink_to_fit();
  int3.shrink_to_fit();
  energyFile5.shrink_to_fit();
  energyPdfFile5.shrink_to_fit();
  energyCdfFile5.shrink_to_fit();
  ein.shrink_to_fit();
  eneE.shrink_to_fit();
  cdf.shrink_to_fit();
  pdf.shrink_to_fit();
  ene2d.shrink_to_fit();
  frac2d.shrink_to_fit();
  cdf2d.shrink_to_fit();
  pdf2d.shrink_to_fit();
  ein2d.shrink_to_fit();
  ene3d.shrink_to_fit();
  cdf3d.shrink_to_fit();
  pdf3d.shrink_to_fit();
}

double TNudyEndfPhEnergy::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(nbt3, int3, nr3, energyFile5, energyPdfFile5, np3, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  energyFile5.push_back(mid);
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5Prob(x1, mid, pdf1, pdf);
  recursionLinearFile5Prob(mid, x2, pdf, pdf2);
  return 0;
}

double TNudyEndfPhEnergy::recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double gx            = 0.0;
  double pe            = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, fE1, fP1, NP, energy);
  double thetae        = TNudyCore::Instance()->Interpolate(nbt2, int2, nr2, fE2, fP2, np2, energy);
  if (thetae > 0.0) gx = TNudyCore::Instance()->Interpolate(nbt3, int3, nr3, fE3, fP3, np3, mid / thetae);
  pdf                  = 0.0;
  pdf                  = pe * gx;
  double pdfmid1       = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  energyFile5.push_back(mid);
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5GenEva(x1, mid, pdf1, pdf, energy);
  recursionLinearFile5GenEva(mid, x2, pdf, pdf2, energy);
  return 0;
}

double TNudyEndfPhEnergy::recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, fE1, fP1, NP, energy);
  double I         = TNudyCore::Instance()->Interpolate(nbt2, int2, nr2, fE2, INorm, np2, energy);
  double thetae    = TNudyCore::Instance()->Interpolate(nbt2, int2, nr2, fE2, fP2, np2, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sqrt(mid) * exp(-mid / thetae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  energyFile5.push_back(mid);
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5Maxwell(x1, mid, pdf1, pdf, energy);
  recursionLinearFile5Maxwell(mid, x2, pdf, pdf2, energy);
  return 0;
}

double TNudyEndfPhEnergy::recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, fE1, fP1, NP, energy);
  double I         = TNudyCore::Instance()->Interpolate(nbt2, int2, nr2, fE2, INorm, np2, energy);
  double ae        = TNudyCore::Instance()->Interpolate(nbt2, int2, nr2, fE2, fP2, np2, energy);
  double be        = TNudyCore::Instance()->Interpolate(nbt2, int2, nr2, fE3, fP3, np3, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sinh(sqrt(be * mid)) * exp(-mid / ae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  energyFile5.push_back(mid);
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5Watt(x1, mid, pdf1, pdf, energy);
  recursionLinearFile5Watt(mid, x2, pdf, pdf2, energy);
  return 0;
}
void TNudyEndfPhEnergy::fillPdf1d()
{
  TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
  TNudyCore::Instance()->ThinningDuplicate(energyFile5, energyPdfFile5);
  TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
  for (unsigned long i = 0; i < energyFile5.size(); i++) {
    if (energyPdfFile5[i] > 1E-15) {
      eneE.push_back(energyFile5[i]);
      pdf.push_back(energyPdfFile5[i]);
      cdf.push_back(energyCdfFile5[i]);
    }
  }
  ene2d.push_back(eneE);
  pdf2d.push_back(pdf);
  cdf2d.push_back(cdf);
  energyFile5.clear();
  energyPdfFile5.clear();
  energyCdfFile5.clear();
  eneE.clear();
  pdf.clear();
  cdf.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfPhEnergy::GetEnergy5(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  // std::cout<<"mt "<< mt <<"  "<< Mt5Values[ielemId].size() <<"  "<<energyK << std::endl;
  for (unsigned int l = 0; l < Mt5Values[ielemId].size(); l++) {
    if (Mt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i <<"  "<< Mt5Values[ielemId][i] <<"  "<<energy5OfMts[ielemId][i].size() << std::endl;
  if (i < 0) return 99;
  int min = 0;
  int max = energy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= energy5OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= energy5OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < energy5OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  if (min < 0) min = 0;
  // for(unsigned int kk = 0; kk < energy5OfMts[ielemId][i].size(); kk++)
  // std::cout << min <<"  "<< energy5OfMts[ielemId][i][kk] <<"  "<< energy5OfMts[ielemId][i].size() << std::endl;
  double fraction =
      (energyK - energy5OfMts[ielemId][i][min]) / (energy5OfMts[ielemId][i][min + 1] - energy5OfMts[ielemId][i][min]);
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  // std::cout<<"min "<< min <<"  "<< energy5OfMts[ielemId][i][min] << std::endl;
  int k    = 0;
  int size = energyCdf5OfMts[ielemId][i][min].size();
  for (unsigned int j = 1; j < energyPdf5OfMts[ielemId][i][min].size(); j++) {
    if (rnd1 <= energyCdf5OfMts[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  // for (unsigned int j1 = 0; j1 < energyOut5OfMts[ielemId][i][min].size(); j1++){
  // std::cout<< energyOut5OfMts[ielemId][i][min][j1] <<"  "<< energyPdf5OfMts[ielemId][i][min][j1] <<std::endl;
  //}
  double plk = (energyPdf5OfMts[ielemId][i][min][k + 1] - energyPdf5OfMts[ielemId][i][min][k]) /
               (energyOut5OfMts[ielemId][i][min][k + 1] - energyOut5OfMts[ielemId][i][min][k]);
  double plk2 = energyPdf5OfMts[ielemId][i][min][k] * energyPdf5OfMts[ielemId][i][min][k];

  double edes = 0;
  if (plk != 0)
    edes = energyOut5OfMts[ielemId][i][min][k] +
           (sqrt(plk2 + 2 * plk * (rnd1 - energyCdf5OfMts[ielemId][i][min][k])) - energyPdf5OfMts[ielemId][i][min][k]) /
               plk;
  return edes;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfPhEnergy::GetDelayedFraction(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  // std::cout<<"mt "<< mt <<"  "<< Mt5Values[ielemId].size() <<"  "<<energyK << std::endl;
  for (unsigned int l = 0; l < Mt5Values[ielemId].size(); l++) {
    if (Mt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i <<"  "<< Mt5Values[ielemId][i] << std::endl;
  if (i < 0) return 99;
  int min = 0;
  int max = energy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= energy5OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= energy5OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < energy5OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fraction5OfMts[ielemId][i][min];
}
