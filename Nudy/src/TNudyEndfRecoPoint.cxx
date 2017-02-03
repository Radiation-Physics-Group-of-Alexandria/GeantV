// This class is reconstructing ENDF data and creating probability tables for the angle
// and energy distributions of the secondatries
// Author: Dr. Harphool Kumawat
// Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// date of creation: March 22, 2016
//Original-> Provides information of neutron only
//SPB, date of modification: December 12, 2016
//Functions are added which return information of secondary products including neutron
#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfList.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfNuPh.h"
#include "TNudyEndfAng.h"
#include "TNudyEndfEnergy.h"
#include "TNudyEndfEnergyAng.h"
#include "TNudyEndfFissionYield.h"
#include "TNudyCore.h"
#include "TNudyEndfPhYield.h"
#include "TNudyEndfPhProd.h"
#include "TNudyEndfPhAng.h"
#include "TNudyEndfPhEnergy.h"
#include "TNudyEndfRecoPoint.h"
#ifdef USE_ROOT
ClassImp(TNudyEndfRecoPoint)
#endif

#ifdef USE_ROOT
#include "TRandom3.h"
#endif

    TNudyEndfRecoPoint::TNudyEndfRecoPoint()
    : elemId(0), rENDF()
{
}
TNudyEndfRecoPoint::TNudyEndfRecoPoint(int ielemId, const char *irENDF)
    : elemId(ielemId), rENDF(irENDF)
{
  GetData(ielemId, irENDF);
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile3(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  int mt1multi = -1;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int MT = sec->GetMT();
    if (MT != 1 && MT != 3 && MT != 4 && MT != 27 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 &&
        (MT < 120 || MT >= 600) ) {
      MtNumbers.push_back(MT);
      TIter recIter(sec->GetRecords());
      TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
      double QI            = header->GetC2();
      QValue[sec->GetMT()] = QI;
      qvaluetemp.push_back(QI);
      qvaluetemp.push_back(MT);
      NR = header->GetN1();
      NP = header->GetN2();
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
      //std::cout<<"ReactionMT: "<<MT<<"   Q-value: "<<QI<<std::endl;
      //std::cout<<"Recopoint: "<<header->GetL1()<<" "<<header->GetL2()<<std::endl;
      for (int crs = 0; crs < NP; crs++) {
        eLinearFile3.push_back(tab1->GetX(crs));
        xLinearFile3.push_back(tab1->GetY(crs));
      }
      eneTemp.insert(std::end(eneTemp),std::begin(eLinearFile3),std::end(eLinearFile3));
      eneTemp.insert(std::end(eneTemp),std::begin(xLinearFile3),std::end(xLinearFile3));

   //SPB
    //auto ai=std::end(eneTemp);
    //auto ai1=std::begin(eLinearFile3);
    //auto ai2=std::begin(xLinearFile3);
    //cout<<"end of eneTemp:::::::::::::\t"<<*ai<<"\t"<<*ai1<<"\t"<<*ai2<<endl;
    //SPB, here in the x-section, last two numbers are half of the input x-section data, has to be corrected
      eLinearFile3.clear();
      xLinearFile3.clear();
      sigmaOfMts.push_back (eneTemp);
      eneTemp.clear();
    } if (MT==1 && mt1multi==-1){
      mt1multi = 0;
      TIter recIter(sec->GetRecords());
      TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
      NR = header->GetN1();
      NP = header->GetN2();
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
      for (int crs = 0; crs < NP; crs++) {
        eLinearFile3.push_back(tab1->GetX(crs));
        xLinearFile3.push_back(tab1->GetY(crs));
      }
      energyMts.insert(std::end(energyMts),std::begin(eLinearFile3),std::end(eLinearFile3));
      sigmaMts.insert(std::end(sigmaMts),std::begin(xLinearFile3),std::end(xLinearFile3));
      eLinearFile3.clear();
      xLinearFile3.clear();     
    }
  }
  MtValues.push_back(MtNumbers);
  qvalue.push_back(qvaluetemp);
  qvaluetemp.clear();
  MtNumbers.clear();
}
void TNudyEndfRecoPoint::GetData(int ielemId, const char *rENDF)
{
  elemId      = ielemId;
  TFile *rEND = TFile::Open(rENDF);
  //  TFile *rEND = TFile::Open(rENDF,"UPDATE");
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey              = (TKey *)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape *)rkey->ReadObj();
  TNudyEndfMat *tMat      = 0;
  TList *mats             = (TList *)rENDFVol->GetMats();
  int nmats               = mats->GetEntries();
  int mt455               = 1000;

  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat *)mats->At(iMat);
    for (int inf = 0; inf < tMat->GetNXC(); inf++) {
      switch (tMat->GetMFn(inf)) {
      case 4:
        MtNum4.push_back(tMat->GetMFn(inf));
        MtNum4.push_back(tMat->GetMTn(inf));
        break;
      case 5:
        MtNum5.push_back(tMat->GetMFn(inf));
        if (tMat->GetMTn(inf) == 455) {
          MtNum5.push_back(mt455);
          mt455++;
        } else {
          MtNum5.push_back(tMat->GetMTn(inf));
        }
        break;
      case 6:
        MtNum6.push_back(tMat->GetMFn(inf));
        MtNum6.push_back(tMat->GetMTn(inf));
        //std::cout<<"MF: "<<tMat->GetMFn(inf)<<" mtn: "<<tMat->GetMTn(inf)<<endl;
        break;
      }
    }
    Mt4.push_back(MtNum4);
    Mt5.push_back(MtNum5);
    Mt6.push_back(MtNum6);
    std::vector<int>().swap(MtNum4);
    std::vector<int>().swap(MtNum5);
    std::vector<int>().swap(MtNum6);
    //    std::cout<<" MAT number "<< tMat->GetMAT() <<"  "<< nmats << std::endl;
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file;
    while ((file = (TNudyEndfFile *)iter.Next())) {
      // std::cout <<" mf "<< file->GetMF() << std::endl;
      // Read File data into class structures
      switch (file->GetMF()) {
      case 1:
        if (flagRead != 1) {
          recoNuPh = new TNudyEndfNuPh(file);
          flagRead = 1;
          std::cout << "file 1 OK: Should be printed only one time " << std::endl;
        }
        break;
      case 2:
        break;
      case 3:
        { 
          ReadFile3(file);
          TIter secIter(file->GetSections());
          TNudyEndfSec *sec;
          while ((sec = (TNudyEndfSec *)secIter.Next())) {
	    int MT = sec->GetMT();
	    if (MT == 1) {
		eneUni.push_back (energyMts);
		sigUniT.push_back (sigmaMts);
	      break;
	    }
	  }
	  fixupTotal(energyMts);
	  energyMts.clear();
	  sigmaMts.clear();
        }break;
      case 4:
        std::cout << "before file 4 " << std::endl;
        recoAng = new TNudyEndfAng(file);
        std::cout << "file 4 OK " << std::endl;
        break;
      case 5:
        std::cout << "before file 5 " << std::endl;
        recoEnergy = new TNudyEndfEnergy(file);
        std::cout << "file 5 OK " << std::endl;
        break;
      case 6:
        std::cout << "before file 6 " << std::endl;
        recoEnergyAng = new TNudyEndfEnergyAng(file, QValue);
        std::cout << "file 6 OK " << std::endl;
        break;
      case 8:
        std::cout << "before file 8 " << std::endl;
        recoFissY = new TNudyEndfFissionYield(file);
        std::cout << "file 8 OK " << std::endl;
        break;
// //       case 12:
// //         std::cout << "before file 12 " << std::endl;
// // 	recoPhYield = new TNudyEndfPhYield(file);
// // 	std::cout<<"file 12 OK "<<std::endl;
// // 	break;
// //       case 13:
// //         std::cout << "before file 13 " << std::endl;
// // 	recoPhProd = new TNudyEndfPhProd(file);
// // 	std::cout<<"file 13 OK "<<std::endl;
// // 	break;
       case 14:
         std::cout << "before file 14 " << std::endl;
 	recoPhAng = new TNudyEndfPhAng(file);
 	std::cout<<"file 14 OK "<<std::endl;
 	break;
// //       case 15:
// //         std::cout << "before file 15 " << std::endl;
// // 	recoPhEnergy = new TNudyEndfPhEnergy(file);
// // 	std::cout<<"file 15 OK "<<std::endl;
// // 	break;
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::fixupTotal(std::vector<double> &x1)
{
  for (unsigned long i = 0; i < sigmaOfMts.size(); i++) {
    int size = sigmaOfMts[i].size()/ 2;
    int flagLoc = -1 ;
    for (unsigned long k = 0; k < x1.size(); k++) {
      int min = 0;
      int max = size - 1;
      int mid = 0;
      //cout<<" x1[k]: "<<x1[k]<<" "<<sigmaOfMts[i][min]<<endl;
      if (x1[k] <= sigmaOfMts[i][min])
        min = 0;
      else if (x1[k] >= sigmaOfMts[i][max])
        min = max;
      else {
        while (max - min > 1) {
          mid = (min + max) / 2;
          if (x1[k] < sigmaOfMts[i][mid])
            max = mid;
          else
            min = mid;
        }
      }
      if (x1[k] == sigmaOfMts[i][min] && sigmaOfMts[i][size + min] > 0) {
        eneTemp.push_back(x1[k]);
        sigTemp.push_back(sigmaOfMts[i][size + min]);
      }
      if (eneTemp.size() == 1 && flagLoc == -1) {
        energyLocationMts.push_back(k);
	flagLoc = 0;
      }
    }
    sigmaUniOfMts.push_back(sigTemp);
    eneTemp.clear();
    sigTemp.clear();
  }
  //SPB
   // for(int i=0;i<sigmaUniOfMts.size();i++){
     // for(int j=0; j<sigmaUniOfMts[i].size(); j++){
    //if(i==0 && j==0)cout<<"Temp energy and x-sections:\t"<<sigmaUniOfMts.size()<<"\t"<<sigmaUniOfMts[i].size()<<endl;
   // } 
   //}
  sigmaOfMts.clear();
  sigUniOfMt.push_back(sigmaUniOfMts);
  energyLocMtId.push_back(energyLocationMts);
  sigmaUniOfMts.clear();
  energyLocationMts.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaTotal(int ielemId, double energyK)
{
  int min = 0;
  int max = eneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= eneUni[ielemId][min]){
    min = 0;
  }else if (energyK >= eneUni[ielemId][max]){
    min = max - 1;
  }else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eneUni[ielemId][mid]){
        max = mid;
      }else{
        min = mid;
       } 
    }
  }
  return sigUniT[ielemId][min] +
         (sigUniT[ielemId][min + 1] - sigUniT[ielemId][min]) * (energyK - eneUni[ielemId][min]) /
             (eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaPartial(int ielemId, int i, double energyK)
{
  //std::cout<< eneUni[ielemId].size() <<std::endl;
  int min = 0;
  int max = eneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= eneUni[ielemId][min])
    min = 0;
  else if (energyK >= eneUni[ielemId][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eneUni[ielemId][mid])
        max = mid;
      else
        min = mid;
    }
  }
  int minp = min - energyLocMtId[ielemId][i];
//     for (unsigned int j1 = 0; j1 < eneUni[ielemId].size(); j1++)
//       std::cout << eneUni[ielemId][j1] <<"  "<< sigUniOfMt[ielemId][i][j1] <<"  "<< sigUniOfMt[ielemId][i].size() << std::endl;
  if (minp <= 0) return 0;
  if (minp >= (int)sigUniOfMt[ielemId][i].size()) return 0;
  return sigUniOfMt[ielemId][i][minp] +
         (sigUniOfMt[ielemId][i][minp + 1] - sigUniOfMt[ielemId][i][minp]) * (energyK - eneUni[ielemId][min]) /
             (eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos6(int ielemId, int mt, double energyK)
{
  return recoEnergyAng->GetCos6(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy6(int ielemId, int mt, double energyK)
{
  return recoEnergyAng->GetEnergy6(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos64(int ielemId, int mt, double energyK)
{
  return recoEnergyAng->GetCos64(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos4(int ielemId, int mt, double energyK)
{
//std::cout<<"recoAng->GetCos4(ielemId, mt, energyK): "<<recoAng->GetCos4(ielemId, mt, energyK)<<std::endl;
  return recoAng->GetCos4(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy5(int ielemId, int mt, double energyK)
{
  return recoEnergy->GetEnergy5(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetAd6(int ielemId, int mt)
{
  return recoEnergyAng->GetAd6(ielemId, mt);
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetZd6(int ielemId, int mt)
{
//  std::cout<<" recoEnergyAng->GetZd6(ielemId, mt): " <<recoEnergyAng->GetZd6(ielemId, mt)<<std::endl;
 int productZ=0;
     productZ = recoEnergyAng->GetZd6(ielemId, mt);
 return productZ;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetLaw6(int ielemId, int mt)
{
  return recoEnergyAng->GetLaw6(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetQValue(int ielemId, int mt)
{
  int size = qvalue[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (qvalue[ielemId][2 * i + 1] == mt) {
      return qvalue[ielemId][2 * i];
    return qvalue[ielemId][2 * i];
    }
  }
  return 0;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt4(int ielemId, int mt)
{
 //std::cout<<"TNudyEndfRecoPoint::GetMt4: "<<Mt4.size()<<std::endl;
  if (Mt4.size() <= 0) return 99;
  int size = Mt4[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (Mt4[ielemId][2 * i + 1] == mt) {
     //std::cout<<"Mt4[ielemId][2 * i + 1] " <<Mt4[ielemId][2 * i + 1]<<std::endl;
      return Mt4[ielemId][2 * i];
    }
  }
  return 99;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt5(int ielemId, int mt)
{
  if (Mt5.size() <= 0) return -1;
  int size = Mt5[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (Mt5[ielemId][2 * i + 1] == mt) {
      return Mt5[ielemId][2 * i];
    }
  }
  return -1;
}
//-------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetMt6Neutron(int ielemId, int mt)
{
  return recoEnergyAng->GetMt6Neutron(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt6(int ielemId, int mt)
{ 
  //std::cout<<"TNudyEndfRecoPoint::GetMt6: "<<Mt6.size()<<std::endl;
  if (Mt6.size() <= 0) return 99;
  int size = Mt6[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (Mt6[ielemId][2 * i + 1] == mt) {
    //std::cout<<"file: "<<Mt6[ielemId][2 * i]<<std::endl;
      return Mt6[ielemId][2 * i];
    }
  }
  return 99;
}
//-------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetCos4Lct(int ielemId, int mt)
{
  return recoAng->GetCos4Lct(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuTotal(int elemid, double energyK)
{
  return recoNuPh->GetNuTotal(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuPrompt(int elemid, double energyK)
{
  return recoNuPh->GetNuPrompt(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuDelayed(int elemid, double energyK)
{
  return recoNuPh->GetNuDelayed(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetFissHeat(int elemid, double energyK)
{
  return recoNuPh->GetFissHeat(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetFisYield(int elemid, double energyK)
{
  return recoFissY->GetFisYield(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetLambdaD(int ielemId, int time)
{
  return recoNuPh->GetLambdaD(ielemId, time);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetDelayedFraction(int ielemId, int mt, double energyK)
{
  return recoEnergy->GetDelayedFraction(ielemId, mt, energyK);
}
//SPB
//______________________________________________________________________________
double TNudyEndfRecoPoint::GetkineticEnergyFF(int elemid, double energyK)
{
 return recoNuPh->GetkineticEnergyFF(elemid, energyK);
} 
//______________________________________________________________________________
double TNudyEndfRecoPoint::GetPromptGammaEnergyFission(int elemid, double energyK)
{
 return recoNuPh->GetPromptGammaEnergyFission(elemid, energyK);
} 
//______________________________________________________________________________
double TNudyEndfRecoPoint::GetDelayedGammaEnergyFission(int elemid, double energyK)
{
 return recoNuPh->GetDelayedGammaEnergyFission(elemid, energyK);
} 
//______________________________________________________________________________
double TNudyEndfRecoPoint::GetDelayedBetaEnergyFission(int elemid, double energyK)
{
 return recoNuPh->GetDelayedBetaEnergyFission(elemid, energyK);
} 
//______________________________________________________________________________
double TNudyEndfRecoPoint::GetNeutrinoEnergyFission(int elemid, double energyK)
{
 return recoNuPh->GetNeutrinoEnergyFission(elemid, energyK);
} 
//______________________________________________________________________________
double TNudyEndfRecoPoint::GetTotalEnergyFission(int elemid, double energyK)
{
 return recoNuPh->GetTotalEnergyFission(elemid, energyK);
} 
//______________________________________________________________________________
//SPB
//products law
std::vector<int> TNudyEndfRecoPoint::GetLawProducts(int ielemId, int mt)
{
 recoEnergyAng->GetLaw6V(ielemId,  mt);
 //std::cout<<" from recopoint : "<<(recoEnergyAng->GetInformation()).lawZA.size()<<endl;
 return (recoEnergyAng->GetInformation()).lawZA;
}
//______________________________________________________________________________
//SPB
//products charge, Z
std::vector<int> TNudyEndfRecoPoint::GetProductsZ(int ielemId, int mt)
{
 recoEnergyAng->GetZd6V(ielemId,  mt);
 //std::cout<<" from recopoint : "<<(recoEnergyAng->GetInformation()).pZ.size()<<endl;
 return (recoEnergyAng->GetInformation()).pZ;
}
//______________________________________________________________________________
//SPB
//Products mass number, A
std::vector<int> TNudyEndfRecoPoint::GetProductsA(int ielemId, int mt)
{
 recoEnergyAng->GetAd6V(ielemId,  mt);
 //std::cout<<" from recopoint : "<<(recoEnergyAng->GetInformation()).pZ.size()<<endl;
 return (recoEnergyAng->GetInformation()).pA;
}
//______________________________________________________________________________
//SPB
//Products Energy
std::vector<double> TNudyEndfRecoPoint::GetProductsEnergy(int ielemId, int mt, double energyK)
{
 recoEnergyAng->GetEnergy6V(ielemId,  mt, energyK);
 //std::cout<<" from recopoint : "<<(recoEnergyAng->GetInformation()).pZ.size()<<endl;
 return (recoEnergyAng->GetInformation()).Energy6;
}
//______________________________________________________________________________
//SPB
std::vector<double> TNudyEndfRecoPoint::GetProductscosAng6(int ielemId, int mt, double energyK)
{
 recoEnergyAng-> GetCos6V(ielemId,  mt, energyK);
 //std::cout<<" from recopoint : "<<(recoEnergyAng->GetInformation()).pZ.size()<<endl;
 return (recoEnergyAng->GetInformation()).cosang6;
}
//______________________________________________________________________________
//SPB
std::vector<double> TNudyEndfRecoPoint::GetProductscosAng64(int ielemId, int mt, double energyK)
{
 recoEnergyAng-> GetCos64V(ielemId,  mt, energyK);
 //std::cout<<" from recopoint : "<<(recoEnergyAng->GetInformation()).pZ.size()<<endl;
 return (recoEnergyAng->GetInformation()).cosang64;
}
//______________________________________________________________________________
TNudyEndfRecoPoint::~TNudyEndfRecoPoint()
{
  sigmaOfMts.shrink_to_fit();
  sigmaUniOfMts.shrink_to_fit();
  energyLocationMts.shrink_to_fit();
  MtNumbers.shrink_to_fit();
  MtNum4.shrink_to_fit();
  MtNum5.shrink_to_fit();
  MtNum6.shrink_to_fit();
  sigmaMts.shrink_to_fit();
  qvaluetemp.shrink_to_fit();
  eLinearFile3.shrink_to_fit();
  xLinearFile3.shrink_to_fit();
  eneTemp.shrink_to_fit();
  sigTemp.shrink_to_fit();
}
//-------------------------------------------------------------------------------------------------------
//SPB
void TNudyEndfRecoPoint::ResetRecoEnergyAng()
{
  recoEnergyAng->Reset();
}
