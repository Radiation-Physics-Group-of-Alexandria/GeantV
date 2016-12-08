// This class is reconstructing ENDF thermal scattering cross-section data and rewrite to ROOT file
#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfList.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfDoppler.h"
#include "TNudyCore.h"
#include "TNudyEndfThermal.h"
#include "Math/SpecFuncMathMore.h"
using namespace std;

#ifdef USE_ROOT
ClassImp(TNudyEndfThermal)
#endif

#ifdef USE_ROOT
#include "TRandom3.h"
#endif

TNudyEndfThermal::TNudyEndfThermal():rENDF(), sigDiff(0)
{}

TNudyEndfThermal::TNudyEndfThermal(const char *irENDF, double isigDiff) : rENDF(irENDF), sigDiff(isigDiff)
{
  GetData(irENDF, isigDiff);
}
void TNudyEndfThermal::GetData(const char *rENDF, double isigDiff)
{
  sigDiff = isigDiff;
  //  TFile *rEND = TFile::Open(rENDF);
  TFile *rEND = TFile::Open(rENDF, "UPDATE");
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey              = (TKey *)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape *)rkey->ReadObj();
  TNudyEndfMat *tMat      = 0;
  TList *mats             = (TList *)rENDFVol->GetMats();
  int nmats               = mats->GetEntries();
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat *)mats->At(iMat);
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file;
    while ((file = (TNudyEndfFile *)iter.Next())) {
      switch (file->GetMF()){
        case 7:
        if(flagRead != 7) {
          ReadFile7(file);
          flagRead = 7;
          std::cout << "file 7 OK: Should be printed only one time " << std::endl;
        }
        break;
      }
    }
  }
  rENDFVol->Write();
  rEND->Close();
 }
//______________________________________________________________________________

double TNudyEndfThermal::recursionLinearFile7(double x1, double x2, double sig1, double sig2,
                                              std::vector<double> ene, std::vector<double> sig)
{
  //x:alpha, siga:S(alpha,T)
  double siga;
  double mid = 0.5 * (x1 + x2);
  //cout<<"x1::\t"<<x1<<"\t"<<x2<<"\t"<<mid<<endl;
  if (sig1 == 0.0) return 0;
  if (sig1 == 0.0 && sig2 == 0.0) return 0;
  if (x1 == x2 || x1 < 1E-5 || x2 < 1E-5) {
    return 0;
  }
  siga           = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, ene, sig, NP, mid);
  // cout<<mid<<"\t"<<siga<<endl;
  if (siga == 0.0) return 0;
  double sigmid1 = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  //cout<<"sigmid1:\t"<<sig1<<"\t"<<sig2<<endl;
  if (fabs((siga - sigmid1) / sigmid1) <= sigDiff) {
    return 0;
  }else{
    eLinElastic.push_back(mid);
    sLinElastic.push_back(siga);
   }
   //recursionLinearFile7(x1, mid, sig1, siga, ene, sig);
   //recursionLinearFile7(mid, x2, siga, sig2, ene, sig);
   return 0;
}
//______________________________________________________________________________
void TNudyEndfThermal::ReadFile7(TNudyEndfFile *file)
{
  for (int i=0;i<1000;i++){
    for (int j=0;j<1000;j++){
      tSet[i][j]=0.0;
    }
  }
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  int thLI;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int MT                = sec->GetMT();
    NR                    = header->GetN1();
    //NP = header->GetN2();
    int LI  = header->GetL1();
    int LCT = header->GetL2();
    int LTT = sec->GetL2();
    int LTHR = sec->GetL1();
    MtNumbers.push_back(sec->GetMT());
    mtLTHR.push_back(LTHR);
    //cout<<header->GetC1()<<"\t"<<header->GetC2()<<"\t"<<header->GetL1()<<"\t"<<header->GetL2()<<"\t"<<header->GetN1()<<"\t"<<header->GetN2()<<endl;
    cout<<sec->GetC1()<<"\t"<<sec->GetC2()<<"\t"<<sec->GetL1()<<"\t"<<sec->GetL2()<<"\t"<<sec->GetN1()<<"\t"<<sec->GetN2()<<"\t"<<sec->GetMAT()<<"\t"<<sec->GetMF()<<"\t"<<sec->GetMT()<<endl;
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
    if(MT==2){//Elastic scattering process
      //std::cout<<"energy "<< tab1->GetC2() <<"\t"<<tab1->GetNR()<<"\t"<<tab1->GetNP()<<std::endl;
      //ein.push_back(tab->GetC2());
      NR = tab1->GetNR();
      NP = tab1->GetNP();
      //if(sec->GetL1()==1){ //GetL1()==LTHR, 1->for coherent elastic scattering
      switch(LTHR) {
        case 1 :{ //GetL1()==LTHR, 1->for coherent elastic scattering
          for (int i = 0; i < tab1->GetNR(); i++) {
            nbt1.push_back(tab1->GetNBT(i));
            int1.push_back(tab1->GetINT(i));
            //cout<<tab1->GetNBT(i)<<"\t"<<tab1->GetINT(i)<<endl;
          }
          //cout<<"temp:\t"<<tab1->GetC1()<<"\t"<<endl;
          temp_T.push_back(tab1->GetC1());//filling temp value from tab1
          for (int j = 0; j < NP; j++) {
            eLinElastic.push_back(tab1->GetX(j));//energy of tab1
            sLinElastic.push_back(tab1->GetY(j));//S(E,T) in tab1
            dummy_eElastic.push_back(tab1->GetX(j));
            dummy_sElastic.push_back(tab1->GetX(j));
            //cout<<"E:\t"<<eLinElastic[j]<<"\t"<<sLinElastic[j]<<endl;
          }
          //cout<<tab1->GetL1()<<"\t"<<endl;
          Nbeta.push_back(tab1->GetL1());// No. of temperature
          for(int k=0; k<tab1->GetL1(); k++){
            TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
            cout<<tab->GetC1()<<"\t"<<endl;
            temp_T.push_back(tab->GetC1()); //filling temp value from list
            for(int l=0; l<tab->GetNPL(); l++){
              tSet[k][l]= tab->GetLIST(l);
              //sLinElastic.push_back(tab->GetLIST(l));//S(E,T) in tab1
              //cout<<"List:\t"<<tab->GetLIST(l)<<endl;
            }
          }
          for (int is = 0; is<NP-1; is++) {// NP - 1
            //std::cout<<"is = "<<is<<"  "<<eLinElastic[is] <<"  "<< sLinElastic[is] << std::endl;
            recursionLinearFile7(eLinElastic[is], eLinElastic[is + 1], sLinElastic[is], sLinElastic[is + 1],
                         eLinElastic, sLinElastic);
            //
           //cout<<"interpolated x-section:\t"<<  recursionLinearFile7(eLinElastic[is], eLinElastic[is + 1], sLinElastic[is], sLinElastic[is + 1],eLinElastic, sLinElastic)<<endl;
          }
          TNudyCore::Instance()->Sort(eLinElastic, sLinElastic);
          TNudyCore::Instance()->ThinningDuplicate(eLinElastic, eLinElastic);
          nbt1.clear();
          int1.clear();
        }break;
        case 2 : { //GetL1()==LTHR, 2->for incoherent elastic scattering
          for (int i = 0; i < tab1->GetNR(); i++) {
            nbt1.push_back(tab1->GetNBT(i));
            int1.push_back(tab1->GetINT(i));
            //cout<<"case2:\t"<<tab1->GetNBT(i)<<"\t"<<tab1->GetINT(i)<<endl;
          }
          for(int j = 0; j < NP; j++) {
            sigb=tab1->GetC1(); //bound x-section in barns
            temp.push_back(tab1->GetX(j));
            DW.push_back(tab1->GetY(j));
            //cout<<"LTHR:2\t"<<tab1->GetX(j)<<"\t"<<tab1->GetY(j)<<endl;
          }
        }//LTHR, 2
      }//switch for LTHR
    }//MT=2
//______________________________________________________________________________
//Next we are looking for incoherent inelastic scattering
    if (MT==4){
      // cout<<"NR:\t"<<NR<<"\t"<<NP<<endl;
      TNudyEndfList *list = (TNudyEndfList *)(sec->GetRecords()->At(0));
      //B(0) = total free stom x-section, M0sigf0
      //B(1) = E/kT
      //B(2) = A0, used to estimate the values of alpha
      //B(3) = E_max , upper energy limit for the constant sigma_f0
      //B(4) = not used
      //B(5) = M0, the number of principal scattering atoms in the material
      int LAT   = sec->GetL2();  //flag for temperature used for estimating alpha & beta
      int LASYM = sec->GetN1();  //flag for symmetric(0) or asymmetric(1) S(alpha,beta)
      int LLN   = list->GetL1(); //flag for how S(alpha,beta) stored, directly(0), stored in ln(S)->1
      switch (LLN) {
        case 0:{ //S(alpha,beta) stored, directly(0)
          for(int ib=0;ib<list->GetNPL(); ib++){
            //cout<<"list:\t"<<list->GetLIST(ib)<<endl;// B(N) values,
            BN.push_back(list->GetLIST(ib));//list of constants
          }
        }break;
        case 1:{ // S(alpha,beta) stored in ln(S)
          for(int ib=0;ib<list->GetNPL(); ib++){
            //cout<<"list:\t"<<list->GetLIST(ib)<<endl;// B(N) values,
            BN.push_back(list->GetLIST(ib));//list of constants
          }
        }break;
      }
      for(int ibn = 0; ibn < BN.size(); ibn++){
        // cout<<"printing the values of BN:\t"<<BN[ibn]<<endl;
      }
      switch (LAT) {
        case 0:{
          //if(LAT==1){
          TNudyEndfTab2 *Tab2 = (TNudyEndfTab2 *)recIter.Next();
          //cout<<"NB:number of beta values:\t"<<Tab2->GetN2()<<"\t"<<Tab2->GetNR()<<"\t"<<Tab2->GetNZ()<<"\tnpl from list:\t"<<list->GetNPL()<<endl;
          int NB=Tab2->GetN2(); //NB(number of beta values)
          for (int i = 0; i < NB; ++i) {//provides 96
            //cout<<"value of NR:\t"<<NR<<endl;
            TNudyEndfTab1 *Tab11 = (TNudyEndfTab1 *)recIter.Next();
            NR = Tab11->GetNR();
            //NP = Tab11->GetNP();
            for (int i = 0; i < Tab11->GetNR(); i++) {
              nbt1.push_back(Tab11->GetNBT(i));
              int1.push_back(Tab11->GetINT(i));
              //cout<<"MT:4\t"<<Tab11->GetNBT(i)<<"\t"<<Tab11->GetINT(i)<<endl;
            }
            //cout<<"Temp:\t"<<Tab11->GetC1()<<endl;
            switch(LASYM){
              case 0:{ // S(alpha,beta) is symmetric
                for (int it1 = 0; it1 < Tab11->GetNP(); it1++) {
                  eIncInelastic.push_back(Tab11->GetX(it1));
                  sIncInelastic.push_back(Tab11->GetY(it1));
                  //cout<<"T and beta:\t"<<Tab11->GetC1()<<"\t"<<Tab11->GetC2()<<endl;//
                  // cout<<Tab11->GetX(it1)<<"  "<<Tab11->GetY(it1)<<endl;
                }
              }break;
              case 1:{ //S(alpha,beta) is asymmetric
                NR = Tab11->GetNR();
                //NP = Tab11->GetNP();
                for (int i = 0; i < Tab11->GetNR(); i++) {
                  nbt2.push_back(Tab11->GetNBT(i));
                  int2.push_back(Tab11->GetINT(i));
                  //cout<<"MT:4\t"<<Tab11->GetNBT(i)<<"\t"<<Tab11->GetINT(i)<<endl;
                }
                for(int it1 = 0; it1 < Tab11->GetNP(); it1++) {
                  //cout<<"NB:number of beta values:\t"<<Tab11->GetL1()<<"\t"<<Tab11->GetNP()<<"\t"<<Tab11->GetN2()<<endl;
                  //cout<<"values of temperaure for whcih S(alpha,beta1,T1) are given:\t"<<Tab11->GetC1()<<endl;
                  //cout<<"get energy:\t"<<"\t"<<Tab11->GetX(it1)<<"\t"<<Tab11->GetY(it1)<<"\t"<<endl;
                  eIncInelastic.push_back(Tab11->GetX(it1));
                  sIncInelastic.push_back(Tab11->GetY(it1));
                }// alpha vs S(alpha,beta)
              }break; //1->asymmetric
            }//LASYM
            for(int j=0; j<Tab11->GetL1(); j++){//loop for T
              TNudyEndfList *listT1 = (TNudyEndfList *)recIter.Next();
              for (int it1 = 0; it1 < Tab11->GetNP(); it1++) {
                //cout<<"list values case0:\t"<<listT1->GetLIST(it1)<<endl;
                //cout<<"T and beta:\t"<<listT1->GetC1()<<"\t"<<listT1->GetC2()<<endl;//
              }
            }//T
          }//NB
        }break; //LAT, case 0
        case 1: {
          TNudyEndfTab2 *Tab2 = (TNudyEndfTab2 *)recIter.Next();
          //cout<<"NB:number of beta values:\t"<<Tab2->GetN2()<<"\t"<<Tab2->GetNR()<<"\t"<<Tab2->GetNZ()<<"\tnpl from list:\t"<<list->GetNPL()<<endl;
          int NB=Tab2->GetN2(); //NB(number of beta values)
          for (int i = 0; i < NB; ++i) {//provides 96
            TNudyEndfTab1 *Tab11 = (TNudyEndfTab1 *)recIter.Next();
            //cout<<"NB:number of beta values:\t"<<Tab11->GetL1()<<"\t"<<Tab11->GetNP()<<"\t"<<Tab11->GetN2()<<endl;//
            //cout<<Tab11->GetC1()<<"\t"<<Tab11->GetC2()<<endl;
            beta[i]=Tab11->GetC2();
            NR = Tab11->GetNR();
            //NP = Tab11->GetNP();
            for(int i = 0; i < Tab11->GetNR(); i++){
              nbt3.push_back(Tab11->GetNBT(i));
              int3.push_back(Tab11->GetINT(i));
              //cout<<"LAT:1, Mt:4\t"<<Tab11->GetNBT(i)<<"\t"<<Tab11->GetINT(i)<<endl;
            }
            //std::cout<<"energy "<< Tab11->GetC1() <<"\t"<<Tab11->GetNR()<<"\t"<<Tab11->GetNP()<<std::endl;
            switch(LASYM) {
              case 0:{ // S(alpha,beta) is symmetric
                for (int i = 0; i < Tab11->GetNR(); i++){
                  nbt1.push_back(Tab11->GetNBT(i));
                  int1.push_back(Tab11->GetINT(i));
                  //cout<<Tab11->GetNBT(i)<<"\t"<<Tab11->GetINT(i)<<endl;
                }
                for(int it1 = 0; it1 < Tab11->GetNP();it1++){
                  //cout<<"NB:number of beta values:\t"<<Tab11->GetL1()<<"\t"<<Tab11->GetNP()<<"\t"<<Tab11->GetN2()<<endl;
                  // cout<<"values of temperaure for whcih S(alpha,beta1,T1) are given:\t"<<Tab11->GetC1()<<endl;
                  //cout<<Tab11->GetX(it1)<<"\t"<<Tab11->GetY(it1)<<endl;
                  //cout<<"T and beta:\t"<<Tab11->GetC1()<<"\t"<<Tab11->GetC2()<<endl;//
                  eIncInelastic.push_back(Tab11->GetX(it1));
                  sIncInelastic.push_back(Tab11->GetY(it1));
                  tempsIncInelastic.push_back(Tab11->GetY(it1));
                }
              }break;
             case 1:{ //S(alpha,beta) is asymmetric
              NR = Tab11->GetNR();
              //NP = Tab11->GetNP();
              for (int i = 0; i < Tab11->GetNR(); i++) {
                nbt2.push_back(Tab11->GetNBT(i));
                int2.push_back(Tab11->GetINT(i));
                //cout<<Tab11->GetNBT(i)<<"\t"<<Tab11->GetINT(i)<<endl;
              }
              for (int it1 = 0; it1 < Tab11->GetNP(); it1++) {
                eIncInelastic.push_back(Tab11->GetX(it1));
                sIncInelastic.push_back(Tab11->GetY(it1));
              }// alpha vs S(alpha,beta)
             }break; //1->asymmetric
            }//LASYM
            for(int j=0; j<Tab11->GetL1(); j++){//loop for T
              TNudyEndfList *listT1 = (TNudyEndfList *)recIter.Next();
              for (int it1 = 0; it1 < Tab11->GetNP(); it1++) {
                //cout<<"list values case1:\t"<<listT1->GetLIST(it1)<<endl;
                //cout<<"T and beta:\t"<<listT1->GetC1()<<"\t"<<listT1->GetC2()<<endl;//
                tempsIncInelastic.push_back(listT1->GetLIST(it1));
              }
            }//T
          }//NB
        }break; //LAT, case 1
      }
      TNudyEndfTab1 *Teff = (TNudyEndfTab1 *)recIter.Next(); // tab1 for effective temperature
      for (int i = 0; i < Teff->GetNR(); i++) {
        nbt3.push_back(Teff->GetNBT(i));
        int3.push_back(Teff->GetINT(i));
        //cout<<Teff->GetNBT(i)<<"\t"<<Teff->GetINT(i)<<endl;
      }
      for (int ief = 0; ief < Teff->GetNP(); ief++) {
        //cout<<"get energy:\t"<<"\t"<<Teff->GetX(ief)<<"\t"<<Teff->GetY(ief)<<"\t"<<endl;
      }
    }//LAT
 }//loop over sec
   int jc=0;
   int count=0;
   double out[100];
   double sap[200];
   for(int k=0; k <200;k++){
     sap[k]=0.0;
   }
   ofstream fout;
   fout.open("HinZrH.txt",ios::out);
   int j;
   for(int i =0; i<eIncInelastic.size();i++){
     fout<<eIncInelastic[i]<<"\t"<<sIncInelastic[i]<<endl;
     //jc=i;
     count=count+1;
     if(count%48 == 0)jc=jc+1;
     // cout<<" size: "<< jc<<"  "<<eIncInelastic[i]<<" "<<sIncInelastic[i]<<endl;//.size()/(48*8)<<endl;
     out[0] = 0;
     for(int ii=0;ii<200;ii++){
       j=ii*48;
       for (j = 0; j < 48; j++) {
         out[j] = out[j] + (eIncInelastic[j+1] - eIncInelastic[j]) * (sIncInelastic[j+1] + sIncInelastic[j]) * 0.5;
       }
       sap[ii]= out[j];
      // if(i== ii)cout<<i<<" "<<sap[ii]<<" "<<endl;
     }
   }
}//readfile7 loop
//______________________________________________________________________________
double TNudyEndfThermal::nThermalElasticXsecion(double En)
{
  double nXsec=0;
  double energyK = En;
  int ielemId=0;
  int min = 0;
  int mid = 0;
  int max;
  for(int mlthr=0; mlthr<mtLTHR.size(); mlthr++){
    // cout<<mtLTHR[mlthr]<<"  "<< MtNumbers.size()<< " "<<MtNumbers[mlthr]<<endl;
    if(MtNumbers[mlthr] == 2 && mtLTHR[mlthr] == 1){
    cout<<"*****Thermal neutron coherent elastic scattering cross-section available****"<<endl;
      for(int j=0; j<1;j++){ //temp_T.size()-1
        for(int k=0; k<eLinElastic.size();k++){
          if(j==0){
            tempEthermal.push_back(eLinElastic[k]);
            tempThermalXsec.push_back(sLinElastic[k]/eLinElastic[k]);
          }
        }
      }
      eneUni.push_back(tempEthermal);
      nThermalElasticXsec.push_back(tempThermalXsec);
      //tempThermalXsec.clear();
      //cout<<"values of NP:\t"<<NP<<endl;
      max = eneUni[ielemId].size() - 1;
      mid = 0;
      if (energyK <= eneUni[ielemId][min]){
        min = 0;
      }else if (energyK >= eneUni[ielemId][max]){
        min = max - 1;
       }else {
         while (max - min > 1) {
           mid = (min + max) / 2;
           //cout<<min<<"\t"<<max<<"\tmid:\t"<<mid<<endl;
           //cout<<"energy: "<<eneUni[ielemId][mid]<<endl;
           //cout<<energyK <<"  "<<eneUni[ielemId][mid]<<endl;
           if (energyK < eneUni[ielemId][mid]){
             max = mid;
             //cout<<"max: "<<max<<endl;
           }else{
             min = mid;
             //cout<<"min: "<<min<<endl;
            }
         }
       }
       nXsec= nThermalElasticXsec[ielemId][min] +
       (nThermalElasticXsec[ielemId][min + 1] - nThermalElasticXsec[ielemId][min])*(energyK - eneUni    [ielemId][min]) /(eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
    }//if loop MT and LTHR
    if(MtNumbers[mlthr] == 2 && mtLTHR[mlthr] == 2){
      cout<<"*****Thermal neutron incoherent elastic scattering cross-section available****"<<endl;  
      double En_temp = energyK ; // temporary energy of neutron
      double fact;
      double mu0;
      int nb=1000;
      double fact1,fact2,fact3,mui;
      double muav;
      double termc,term2,term3,term4;
      r1=new TRandom3();
      for(int j=0; j<temp.size();j++){ // temperature loop
         for(int k=0; k<DW.size();k++){ //DW loop
            if(j==0 && k==0){
              mu0=-1.0;
              fact= En_temp*DW[k];
              //cout<<"Debye Waller factor:\t"<<DW[k]<<endl;
              sig_ies=(sigb/2)*(1-exp(-4*fact))/(2*fact);
              //cout<<En_temp<<"\t"<<(1-exp(-4*fact))/(2*fact)<<"\t"<<sig_ies<<endl;
              fact1=(1-exp(-4*fact))/nb;
              fact2=exp(-2*fact*(1-mu0));
              fact3=1+(1/(2*fact))*log(fact1+fact2);//gives mui
              termc=nb/(2*fact);
              term2=exp(-2*fact*(1-fact3))*(2*fact*fact3-1);
              term3=exp(-2*fact*(1-mu0))*(2*fact*mu0-1);
              term4=1-exp(-4*fact);
              muav=termc*(term2-term3)/term4;//provides avg mui
              mu0=fact3;
              // cout<<En_temp<<"\t"<<fact3<<"\t"<<muav<<"\t"<<sig_ies<<endl;
            }
         }
      }
    }
  } // for loop MT and LTHR
//cout<<"the cross section: "<<nXsec<<endl;
return nXsec + sig_ies;
}
//______________________________________________________________________________

double TNudyEndfThermal::fun(double a, double E, double x,double b)
{

double f;
f=a*sqrt(x/E)*exp(-(x-E)/b);
return f;
}
//______________________________________________________________________________

TNudyEndfThermal::~TNudyEndfThermal()
{
//delete r1;
//delete r2;
}
