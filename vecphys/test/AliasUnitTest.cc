/*
 Author: Marilena Bandieramonte
 Email: marilena.bandieramonte@cern.ch
 
 Description:
 This is an unit test to check generated alias tables in VecPhys
 This test uses the Alias table manager to create tables for well-known distributions such as the normal and poisson
 distributions, then check that a histogram of many generated values actually follows the theoretical curve using a simple Chi squared test.
 
 How to run:
 
 */

#include <iostream>
#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include <TLatex.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TApplication.h>

//
#include <TFrame.h>
#include <TF1.h>
#include <TFile.h>
#include <TPaveLabel.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include <sys/time.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TPaveText.h>

#define ENTRIES 1000000

typedef int Random_t;
using namespace vecphys;

//--------------------------------------------------
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//--------------------------------------------------
int sampleBin(int size)
{
    //double r1=fRand(0, size-1);
    TRandom3 *r = new TRandom3(0);
    double r1=r->Uniform(0, size);
    int bin = floor(r1);
    return bin;
}


//--------------------------------------------------
// interactAlias method: void interactAlias(GUAliasTable* fAliasTable, TH1F *histo, int size)
// This method starts from an Alias Table and try to use the Alias method to reconstruct the original pdf
void interactAlias(GUAliasTable* fAliasTable, TH1F *histo, int size, double lowEdge, double upEdge){
    
    TH1F* reconstructedPdf=new TH1F("recHist", "ReconstructedHisto", size, lowEdge, upEdge);
    //TH1F* randomHisto=new TH1F("randomHisto", "randomHisto", 100, -3, 3);
    
    int selectedBin;
    bool takeNonAlias=false;
    int sumAlias=0;
    double randomNumber=0.;
    
    TRandom3 *randomNew = new TRandom3(0);
    
    int entriesRecPdf=0;
    //iterate to sample the bin, and sample from the aliasTable the corresponding entries
    for(int i = 0; i < ENTRIES; i++){
        selectedBin=sampleBin(size);
        if(selectedBin>size || selectedBin<0) {std::cout<<"Error!\n"; return;}
        //std::cout<<"+++++ fAliasTable->fProbQ["<<selectedBin<<"]: "<<fAliasTable->fProbQ[selectedBin]<<"\n";
        
        //randomNumber=fRand(0,1);
        randomNumber=randomNew->Uniform(0, 1);
        //std::cout<<"rand(): "<<randomNumber<<"\n";
        //randomHisto->Fill(randomNumber);
        
        takeNonAlias=(randomNumber <= fAliasTable->fProbQ[selectedBin]); //Non-alias probability, ir random<fProbQ, take the non-Alias
        if(takeNonAlias){
            sumAlias++;
            //std::cout<<" -----> fAliasTable->fProbQ["<<selectedBin<<"]="<<fAliasTable->fProbQ[selectedBin]<<"\n";
            reconstructedPdf->AddBinContent(selectedBin+1);
            entriesRecPdf++;
            takeNonAlias=false;
        }
        else{
            entriesRecPdf++;
            //std::cout<<" -----> EV PROB fAliasTable->fAlias["<<selectedBin<<"]="<<fAliasTable->fAlias[selectedBin]<<"\n";
            reconstructedPdf->AddBinContent(fAliasTable->fAlias[selectedBin]+1);
        }
    }
    
    //std::cout<<"+++++ sumAlias: "<<sumAlias<<"\n";
    TCanvas *c2 = new TCanvas("c2","RecHisto",200,10,700,500);
    reconstructedPdf->SetLineColor(kMagenta);
    
    ////try
    reconstructedPdf->SetEntries(entriesRecPdf);
    Double_t norm = reconstructedPdf->GetEntries();
    std::cout<<"entries: "<<norm<<"\n"; if(norm==0) { std::cout<<"ERROR! no entries in the rec histo\n"; return; }
    reconstructedPdf->Scale(1/norm);
    reconstructedPdf->Draw("same");
    
    //want to know the bin with the maximum value
    //int binmax = reconstructedPdf->GetMaximumBin();
    //double x = reconstructedPdf->GetXaxis()->GetBinCenter(binmax);
    //std::cout<<"BinMax: "<<binmax<<" getBinCenter: "<<x<<"\n";
    
    histo->SetLineColor(kGreen);
    histo->Draw("same");
    
    
    //Fit with a Gaussian
    histo->Fit("gaus","0");
    TF1 *fit1 = (TF1*)histo->GetFunction("gaus");
    fit1->SetLineColor(kRed);
    //fit1->Draw("same");
    
    
    //Fit with a Gaussian
    reconstructedPdf->Fit("gaus","0");
    TF1 *fit2 = (TF1*)reconstructedPdf->GetFunction("gaus");
    fit2->SetLineColor(kBlue);
    //fit2->Draw("same");
    
    std::vector <double> resLo(histo->GetNbinsX()); //empty vector with number of elements equal to the number of bins of histo
    double pvalue=histo->Chi2Test(reconstructedPdf, "WW NORM P", resLo.data()); //calculate the  and store the results in resLo
    
    //Draw the chiSquare value on the Canvas
    /*auto t=new TLatex(0.1, 0.95, Form("p-value: %g ", chi2));
     t->SetTextFont(42);
     t->SetNDC();
     t->SetTextColor(6);
     t->SetTextSize(0.03);
     t->Draw();
     */
    
    //pValueP=eOutG4->Chi2Test(eOutScalar,"UU P",resA1);
    TPaveText* pt3 = new TPaveText(0.13, 0.55, 0.43, 0.75, "brNDC");
    Char_t str[50];
    sprintf(str, "p-value from Pearson chi2 test: %g", pvalue);
    pt3->AddText(str);
    pt3->SetFillColor(18);
    pt3->SetTextFont(20);
    pt3->SetTextColor(4);
    pt3-> Draw();
    //Add the legend to the Canvas
    auto leg = new TLegend(0.1,0.8,0.48,0.9);
    //leg->SetHeader("The Legend Title");
    leg->AddEntry(histo,"Histogram filled with the chosen distribution","l");
    leg->AddEntry(reconstructedPdf,"Reconstructed pdf with Alias method","l");
    leg->Draw();
    
    c2->Update();
    c2->SaveAs("RecHisto.pdf");
}

//--------------------------------------------------
Double_t poissonf(Double_t*x,Double_t*par)
{
    return par[0]*TMath::Poisson(x[0],par[1]);
}



//--------------------------------------------------
// This method creates a gaussian histogram and then applying the Alias method to reconstruct the original pdf
void testGaussian(int i){
    
    int histo_bins=100;
    int histo_entries=ENTRIES;
    TH1F* histo;
    
    switch(i){
            
        case 1 :
        {
            //*********** GAUSSIAN ***********//
            std::cout<<"\nAlias test with Gaussian: START.\n";
            //Create an histogram from a Gaussian distribution, with <histo_bins> bins + underflow and overflow, and fill with <histo_entries> entries
            histo=new TH1F("h1", "histo from a gaussian", histo_bins, -3, 3); // size is 100 BIN+2 (underflow ad overflow bins = from 0 to 101)
            histo->FillRandom("gaus", histo_entries);
            break;
        }
            
        case 2 :{
            //*********** POISSON 1 ***********//
            std::cout<<"\nAlias test with POISSON A: START.\n";
            auto *f1 = new TF1("poisson","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", 0, 10); // "xmin" = 0, "xmax" = 10
            f1->SetParameters(4, 4, 4); // you MUST set non-zero initial values for parameters
            //hPois->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
            histo = new TH1F("h1f","Test random numbers",histo_bins,0,10);
            //histo->SetFillColor(45);
            histo->FillRandom("poisson",histo_entries);
            break;
        }
            
        case 3 :
        {
            
            //*********** POISSON 2 ***********//
            std::cout<<"\nAlias test with POISSON B: START.\n";
            TF1 pois("pois",poissonf,10,1,2); // x in [0;10], 2 parameters
            pois.SetParName(0,"Const");
            pois.SetParName(1,"#mu");
            // Create histogram with poisson distribution
            histo = new TH1F("testhi","Poisson distribution",100,0,5);
            //pois.SetParameter(0,3.75654);
            //pois.SetParameter(1,2.95437);
            
            pois.SetParameter(0,1);
            pois.SetParameter(1,1);
            //histo->FillRandom("pois",20000);
            histo->FillRandom("pois",histo_entries);
            break;
        }
            
        case 4 :{
            
            //*********** SQROOT ***********//
            std::cout<<"\nAlias test with SQROOT: START.\n";
            auto form1 = new TFormula("form1","abs(sin(x)/x)");//x*gaus(0) + [3]*abs(sin(x)/x)
            auto sqroot = new TF1("sqroot","x*gaus(0) + [3]*form1",0,10);
            sqroot->SetParameters(10,4,1,20);
            histo = new TH1F("h1f","Test random numbers",histo_bins,0,10);
            //histo->SetFillColor(45);
            histo->FillRandom("sqroot",histo_entries);
            break;
        }
            
        default :{
            std::cout<<"Invalid value. Exiting.\n";
            exit (-1);
        }
            
    }
    
    //double min=histo->GetMinimum();//min content of bins (i.e. 2916)
    //double max=histo->GetMaximum();//max content of bins (i.e. 240263)
    
    double lowEdge=histo->GetXaxis()->GetBinLowEdge(1);         //low edge on the x-axis (i.e. -3)
    double upEdge=histo->GetXaxis()->GetBinUpEdge(histo_bins);  //up edge on the x-axis (i.e. 3)
    
    //std::cout<<"****** min: "<<min<<" max: "<< max<<"\n";
    //std::cout<<"****** lowEdge: "<<lowEdge<<" upEdge: "<< upEdge<<"\n";
    
    int histo_size=histo->GetSize();
    std::cout<<"Histo size: "<<histo_size<<"\n"; //histo_bins+2
    double *pdf = new double[histo_bins]; //pdf[100] --> 0-99
    double sumHisto=0., sumHistoScaled=0.;
    
    
    //Store the bin contents in the pdf vector, but they must be scaled
    for(int i = 0; i < histo_bins ; ++i){
        pdf[i] = histo->GetBinContent(i+1);         //i+1 because bin[0] is the underflow and bin[histo_size]
        //std::cout<<"<<<<<<<gaussian bin content: "<<pdf[i]<<"\n";
        sumHisto+=pdf[i];
    }
    
    Double_t norm = histo->GetEntries();
    //std::cout<<"norm: "<<norm<<"\n"; //histo_entries
    //Scale of the histogram in order to have the pdf
    
    histo->Scale(1/norm);
    //check the integral
    for(int i = 0; i < histo_bins ; ++i){
        pdf[i] = histo->GetBinContent(i+1);
        sumHistoScaled+=pdf[i];
    }
    
    TCanvas *c1 = new TCanvas("c1","histo",200,10,700,500);
    histo->Draw();
    c1->Update();
    
    std::cout<<"sumHisto: "<<sumHisto<<"\n";
    std::cout<<"sumHistoScaled: "<<sumHistoScaled<<"\n";
    std::cout<<"histo->GetSize(): "<<histo->GetSize()<<"\n";
    
    //Create a GUAliasSampler
    Random_t *fRandomState=new Random_t() ;
    int fThreadId=0;
    double fLowEnergyLimit=10;
    double fHighEnergyLimit=1000;
    
    //0 is because we want to test only one at the moment
    GUAliasSampler* fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fLowEnergyLimit, fHighEnergyLimit, 0, histo_bins);
    int z = 0; // Z=0 is convention for atomic independent models
    //Call void GUAliasSampler::BuildAliasTable(int Zelement, const double *pdf)
    fAliasSampler->BuildAliasTable(z, pdf);
    GUAliasTable* aliasTable=fAliasSampler->GetAliasTableManager()->GetAliasTable(0);
    std::cout<<"aliasTable->SizeOfGrid(): "<<aliasTable->SizeOfGrid()<<"\n";
    std::cout<<"aliasTable->SizeOfTable(): "<<aliasTable->SizeOfTable()<<"\n";
    
    std::cout<<"Calling interact alias\n";
    interactAlias(aliasTable, histo, histo_bins, lowEdge, upEdge);
    std::cout<<"Alias test: END.\n";
    c1->SaveAs("histo.pdf");
    
}

//-----------------------------
void fillrandom() {
    // Fill a 1-D histogram from a parametric function
    // To see the output of this macro, click begin_html <a href="gif/fillrandom.gif">here</a>. end_html
    TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,10,700,900);
    c1->SetFillColor(18);
    
    auto pad1 = new TPad("pad1","The pad with the function",0.05,0.50,0.95,0.95,21);
    auto pad2 = new TPad("pad2","The pad with the histogram",0.05,0.05,0.95,0.45,21);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    
    //gBenchmark->Start("fillrandom");
    
    // A function (any dimension) or a formula may reference
    // an already defined formula
    auto form1 = new TFormula("form1","abs(sin(x)/x)");
    auto sqroot = new TF1("sqroot","x*gaus(0) + [3]*form1",0,10);
    sqroot->SetParameters(10,4,1,20);
    
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->GetFrame()->SetFillColor(42);
    pad1->GetFrame()->SetBorderMode(-1);
    pad1->GetFrame()->SetBorderSize(5);
    sqroot->SetLineColor(4);
    sqroot->SetLineWidth(6);
    sqroot->Draw();
    auto lfunction = new TPaveLabel(5,39,9.8,46,"The sqroot function");
    lfunction->SetFillColor(41);
    lfunction->Draw();
    c1->Update();
    
    // Create a one dimensional histogram (one float per bin)
    // and fill it following the distribution in function sqroot.
    pad2->cd();
    pad2->GetFrame()->SetFillColor(42);
    pad2->GetFrame()->SetBorderMode(-1);
    pad2->GetFrame()->SetBorderSize(5);
    auto h1f = new TH1F("h1f","Test random numbers",200,0,10);
    h1f->SetFillColor(45);
    h1f->FillRandom("sqroot",10000);
    h1f->Draw();
    c1->Update();
    
    // Open a ROOT file and save the formula, function and histogram
    TFile myfile("fillrandom.root","RECREATE");
    form1->Write();
    sqroot->Write();
    h1f->Write();
    //gBenchmark->Show("fillrandom");
}


//--------------------------------------------------
int main(int argc,  char *argv[]){
    
    TApplication theApp("App",&argc,argv);
    
    if(argc<2)
    {
        std::cout<<"Usage: ./AliasUnitTest <index> \n <1: gaussian> <2: poissonA> <3: poissonB> <4: sqroot>\n";
        return -1;
    }
    
    struct timeval time;
    gettimeofday(&time,NULL);
    
    // microsecond has 1 000 000
    // Assuming you did not need quite that accuracy
    // Also do not assume the system clock has that accuracy.
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    
    std::cout<<"Alias unit test: START.\n";
    
    /*
     std::cout<<"Alias unit test: START.\n";
     
     //Create a GUAliasSampler
     Random_t *fRandomState;         // necessary to instantiate the GUAliasSampler
     int fThreadId;                  // necessary to instantiate the GUAliasSampler
     double fLowEnergyLimit=10;      // necessary to instantiate the GUAliasSampler
     double fHighEnergyLimit=1000;   // necessary to instantiate the GUAliasSampler
     GUAliasSampler* fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fLowEnergyLimit, fHighEnergyLimit, 0, 6); // 0: only one energy, 6: n. of bins
     //Quantities to check
     bool negativeNAProb=false;
     bool negativeAliasIndex=false;
     bool errorInAliasValues=false;
     bool chi2Valid=false;
     
     
     //Build the pdf table
     size_t sizeOfTable = (fAliasSampler->GetNumEntries()+1) * fAliasSampler->GetSamplesPerEntry();//why the +1?
     
     std::cout<<"Size of table: "<<sizeOfTable<<"\n";
     double *pdf = new double[sizeOfTable]; //the pdf is just generated and stored in a vector
     
     double sum=0.;
     for(int i=0; i<sizeOfTable; i++)
     {
     srand (1);
     pdf[i]=random()%100;
     sum+=pdf[i];
     std::cout<<"pdf["<<i<<"]: "<<pdf[i]<<"\n";
     }
     double sumNormalized=0.;
     for(int i=0; i<sizeOfTable; i++) //normalize the pdf
     {
     pdf[i]/=sum;
     sumNormalized+=pdf[i];
     std::cout<<"pdfNorm["<<i<<"]: "<<pdf[i]<<"\n";
     
     }
     std::cout<<"Normalized sum: "<<sumNormalized<<"\n";
     
     int z = 0; // Z=0 is convention for atomic independent models
     //static_cast<EmModel *>(this)->BuildPdfTable(z, pdf);
     
     
     //Call void GUAliasSampler::BuildAliasTable(int Zelement, const double *pdf)
     fAliasSampler->BuildAliasTable(z, pdf);
     std::cout<<"Alias Table created.\n";
     
     GUAliasTable* aliasTable=fAliasSampler->GetAliasTableManager()->GetAliasTable(0);
     
     std::cout<<"\n\n\n**********\nDebug of the Alias Table for n.Entries: "<<aliasTable->SizeOfGrid()<<"\n";
     for(int i=0; i<aliasTable->SizeOfGrid(); ++i)
     {
     std::cout<<"fpdf["<<i<<"]: "<<aliasTable->fpdf[i]<<"\n";
     std::cout<<"fProbQ["<<i<<"]: "<<aliasTable->fProbQ[i]<<"\n";
     std::cout<<"fAlias["<<i<<"]: "<<aliasTable->fAlias[i]<<"\n";
     if(aliasTable->fProbQ[i]<0) negativeNAProb=true;
     if(aliasTable->fAlias[i]<0) negativeAliasIndex=true;
     if(aliasTable->fProbQ[i]!=1 && aliasTable->fAlias[i]==i) {errorInAliasValues=true;}
     }
     delete[] pdf;*/
    int index = atoi(argv[1]);
    testGaussian(index);
    theApp.Run();
    
}
