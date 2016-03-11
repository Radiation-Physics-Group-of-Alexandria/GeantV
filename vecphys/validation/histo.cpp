/*
 =======
 Author: M. Bandieramonte
 
 NB:
 1. How to compile:
 g++ histo.cpp -o histo `root-config --cflags --glibs`
 
 2. Then execute the file and give the input interactively
 
 Once defined the projectile input energy (i.e. 100 MeV), the program is looking for the following input files:
 geant4_100MeV.root
 scalar_100MeV.root
 vector_100MeV.root
 =======
 */

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <cmath>
#include <unistd.h>
#include <TMarker.h>
#include <TGraph2D.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TAttText.h>
#include <TGraphQQ.h>
#include <TGraphAsymmErrors.h>

#include <stdio.h>
#include <string.h>
#include <algorithm>


#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TH3D.h"
#include "TF1.h"
#include "TImage.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <stdlib.h>
#include "TGaxis.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TLatex.h"

//got
#include <cassert>
#include "Math/GoFTest.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"
//

constexpr double electron_mass_c2 = 0.510998910 * 1.;
constexpr double inv_electron_mass_c2 = 1.0/electron_mass_c2;

using namespace ROOT::Math;

#define DEBUG



//_________
void MSaveBigPDF(double scale=5) {
    TCanvas* old_canv = gPad->GetCanvas();
    
    gROOT->SetBatch(true);
    gROOT->ForceStyle(true);
    
    Int_t orig_msz = gStyle->GetMarkerSize();
    Int_t orig_mst = gStyle->GetMarkerStyle();
    Int_t orig_lt  = gStyle->GetLineWidth();
    
    gStyle->SetMarkerSize(1.0+scale/5);
    gStyle->SetMarkerStyle(20);
    gStyle->SetLineWidth(orig_lt*scale);
    
    TString filename = old_canv->GetName();
    filename += ".png";
    
    
    Int_t old_width  = old_canv->GetWindowWidth();
    Int_t old_height = old_canv->GetWindowHeight();
    
    Int_t new_width = old_width * scale;
    Int_t new_height= old_height* scale;
    
    TCanvas* temp_canvas = new TCanvas("temp", "", new_width, new_height);
    old_canv->DrawClonePad();
    
    temp_canvas->SetLineWidth(orig_lt*40);
    temp_canvas->Draw();
    temp_canvas->SaveAs(filename);
    temp_canvas->Close();
    
    gStyle->SetMarkerSize(orig_msz);
    gStyle->SetMarkerStyle(orig_mst);
    gStyle->SetLineWidth(orig_lt);
    
    gROOT->ForceStyle(false);
    gROOT->SetBatch(false);
    
    std::cout<<"Saving the image as: "<<filename<<"\n";
    
    return;
}


//_________
void drawLatex(double x, double y, const char *s)
{
    //TLatex *t = new TLatex(x,y,Form("#chi^{2}: %s",s));
    TLatex *t = new TLatex(x,y,Form("p-value: %s",s));
    t->SetTextFont(42);
    //t->SetTextAlign(12);
    t->SetNDC();
    t->SetTextColor(6);
    t->SetTextSize(0.048);
    t->Draw();
}

//_________
void drawText(double x, double y, const char *s)
{
    TText *t = new TText(x,y,s);
    t->SetTextFont(42);
    //t->SetTextAlign(12);
    t->SetNDC();
    t->SetTextColor(6);
    t->SetTextSize(0.048);
    t->Draw();
}


//______________________________________________________________________________
void validatePdf(TH1F* eOutScalar, TH1F* eOutVector, TH1F* eOutG4, char* energy)
{
    
    Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    //eOutScalar->GetEntries();
    //std::cout<<"N of histogram entries: "<<eOutScalar->GetEntries()<<"\n";
    
    
    Double_t xScalar[entries], yScalar[entries], zScalar[entries];
    Double_t xVector[entries], yVector[entries], zVector[entries];
    Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
    
    TString pdfFileName=Form("pdf_%sMeV.root",energy);
    TFile *fPfd = new TFile(pdfFileName,"w");
    TGraph *pdfGraph= (TGraph*)fPfd->Get("Pdf2.0"); //Get the first pdf, for fMinX value
    
    
    TString c2Name=Form("Pdf %s MeV",energy);
    
    TCanvas *c2 = new TCanvas("c2",c2Name,200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    //pdfGraph->SetTitle("pdf; mcs; dm");
    
    pdfGraph->Draw();
    c2->Update();
    
    
    //Scale the histograms
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    
    
    TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    norm= eOutVectorScaled->GetEntries();
    eOutVectorScaled->Scale(1/norm);
    
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);
    
    
    for(int j = 0; j < entries ; ++j){
        yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xScalar[j],zScalar[j]); //to map
        
        yVector[j] = eOutVectorScaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xVector[j],zVector[j]); //to map
        
        yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xGeant4[j],zGeant4[j]); //to map
    }
    
    
    //Read the pdf File for 500MeV
    TCanvas *c3 = new TCanvas("c3","Graph comparison",200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    pdfGraph->Draw();
    c3->Update();
    
    
    
    //Create the Graph with x-values taken from the pdfGraph and the y-values taken from the simulation histogram
    TGraph *myCopyVectorGr=new TGraph(entries,xVector,yVector);
    //myCopyVectorGr->SetLineColor(kYellow+10);
    //myCopyVectorGr->SetLineStyle(2);
    //myCopyVectorGr->SetLineWidth(2);
    //myCopyVectorGr->Draw("LP");
    
    //TMarker *mark=new TMarker();
    //mark->SetMarkerStyle(2);
    //mark->SetMarkerSize(1);
    
    TGraph *myCopyScalarGr=new TGraph(entries,xScalar,yScalar);
    myCopyScalarGr->SetLineColor(kYellow+10);
    //myCopyScalarGr->SetLineStyle(4);
    //myCopyScalarGr->SetLineWidth(2);
    
    //myCopyScalarGr->SetMarkerStyle(20);
    myCopyScalarGr->Draw("LP");
    
    //TMarker *mark2=new TMarker();
    //mark2->SetMarkerStyle(4);
    
    
    TGraph *myCopyG4Gr=new TGraph(entries,xGeant4,yGeant4);
    myCopyG4Gr->SetLineColor(kRed);
    //myCopyG4Gr->SetLineStyle(5);
    //myCopyG4Gr->SetLineWidth(2);
    myCopyG4Gr->Draw("LP");
    
    TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
    leg->SetHeader("Legend");
    //leg->AddEntry(my,"Histogram filled with random numbers","f");
    //leg->AddEntry(my,"Function abs(#frac{sin(x)}{x})","l");
    //leg->AddEntry("gr","Graph with error bars","lep");
    TString legName=Form("Pdf for E_in=%s MeV",energy);
    leg->AddEntry(pdfGraph,legName,"l");
    leg->AddEntry(myCopyScalarGr,"E_out Histogram values scaled - Scalar","l");
    //leg->AddEntry(myCopyVectorGr,"E_out Histogram values scaled - Vector","l");
    leg->AddEntry(myCopyG4Gr,"E_out Histogram values scaled - Geant4","l");
    
    leg->Draw();
    c3->Update();
    c3->SaveAs("pdfValidation.pdf");
    
    
    ////Calculate chi-square
    std::cout<<"#E_out Histogram entries: "<<myCopyScalarGr->GetN()<<"#Pdf entries: "<<pdfGraph->GetN()<<"\n";
    double chiSquare_vector=0,chiSquare_scalar=0, chiSquare_geant4=0, expX,expY ;
    for (int i=0; i<pdfGraph->GetN(); i++)
    {
        pdfGraph->GetPoint(i, expX, expY);
        //myCopyScalarGr->GetPoint(i, obsX, obsY);
        //std::cout<<"expX: "<<expX<<" expY: "<<expY<<" obsX: "<<obsX<<" obsY: "<<obsY<<std::endl;
        chiSquare_scalar+=(yScalar[i]-expY)*(yScalar[i]-expY)/expY;
        chiSquare_vector+=(yVector[i]-expY)*(yVector[i]-expY)/expY;
        chiSquare_geant4+=(yGeant4[i]-expY)*(yGeant4[i]-expY)/expY;
    }
    std::cout<<"chiSquare scalar: "<<chiSquare_scalar<<" chiSquare-reduced: "<<chiSquare_scalar/(eOutScalarScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare vector: "<<chiSquare_vector<<" chiSquare-reduced: "<<chiSquare_vector/(eOutVectorScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare geant4: "<<chiSquare_geant4<<" chiSquare-reduced: "<<chiSquare_geant4/(eOutGeant4Scaled->GetSize()-2)<<std::endl;
    
}

//______________________________________________________________________________
void chiSquare_pdf(TH1F* eOutScalar, TH1F* eOutVector, TH1F* eOutG4, int energy)
{
    
    Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    //eOutScalar->GetEntries();
    //std::cout<<"N of histogram entries: "<<eOutScalar->GetEntries()<<"\n";
    
    
    Double_t xScalar[entries], yScalar[entries], zScalar[entries];
    Double_t xVector[entries], yVector[entries], zVector[entries];
    Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
    Double_t pdf[entries], x[entries];
    TGraph *pdfGraph;
    
    double logxmin = std::log(1);
    double dx = (std::log(10000) - logxmin)/99;
    
    //// pdf calculation
    double energy0 = std::exp(logxmin + dx*energy);
    
    double ymin = energy0/(1+2.0*energy0*inv_electron_mass_c2); //MeV
    double dy = (energy0 - ymin)/(1000);
    double yo = ymin + 0.5*dy;
    
    
    double sum = 0.;
    //double integralPdf=0;
    
    for(int j = 0; j < entries ; ++j) {
        //for each output energy bin
        double energy1 = yo + dy*j;
        
        double *grej=new double();
        
        
        
        double E0_m = energy0/0.510998910 ;
        double epsilon = energy1/energy0;
        
        double onecost = (1.- epsilon)/(epsilon*E0_m);
        double sint2   = onecost*(2.-onecost);
        double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
        double xsec = (epsilon + 1./epsilon)*greject;
        
        x[j]=energy1;
        pdf[j] = xsec;
        sum += xsec;
    }
    std::cout<<"xsec sum: "<<sum<<" 1/sum: "<<1/sum<<" \n";
    
    sum=1./sum;
    for(int j = 0; j < entries ; ++j)
        pdf[j]*=sum;
    
    pdfGraph = new TGraph(entries,x,pdf);
    
    TCanvas *c1 = new TCanvas("c1","pdf",200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    pdfGraph->Draw();
    c1->Update();
    
    //TString pdfFileName=Form("pdf_%sMeV.root",energy);
    //TFile *fPfd = new TFile(pdfFileName,"w");
    //TGraph *pdfGraph= (TGraph*)fPfd->Get("Pdf2.0"); //Get the first pdf, for fMinX value
    
    //Scale the histograms
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    
    
    TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    norm= eOutVectorScaled->GetEntries();
    eOutVectorScaled->Scale(1/norm);
    
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);
    
    
    for(int j = 0; j < entries ; ++j){
        yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
        //pdfGraph->GetPoint(j, xScalar[j],zScalar[j]); //to map
        
        yVector[j] = eOutVectorScaled->GetBinContent(j+1);
        //pdfGraph->GetPoint(j, xVector[j],zVector[j]); //to map
        
        yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
        //pdfGraph->GetPoint(j, xGeant4[j],zGeant4[j]); //to map
    }
    
    ////Calculate chi-square
    //std::cout<<"#E_out Histogram entries: "<<myCopyScalarGr->GetN()<<"\n";
    double chiSquare_vector=0,chiSquare_scalar=0, chiSquare_geant4=0, expX,expY ;
    for (int i=0; i<entries; i++)
    {
        chiSquare_scalar+=(yScalar[i]-pdf[i])*(yScalar[i]-pdf[i])/pdf[i];
        chiSquare_vector+=(yVector[i]-pdf[i])*(yVector[i]-pdf[i])/pdf[i];
        chiSquare_geant4+=(yGeant4[i]-pdf[i])*(yGeant4[i]-pdf[i])/pdf[i];
    }
    std::cout<<"chiSquare scalar: "<<chiSquare_scalar<<" chiSquare-reduced: "<<chiSquare_scalar/(eOutScalarScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare vector: "<<chiSquare_vector<<" chiSquare-reduced: "<<chiSquare_vector/(eOutVectorScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare geant4: "<<chiSquare_geant4<<" chiSquare-reduced: "<<chiSquare_geant4/(eOutGeant4Scaled->GetSize()-2)<<std::endl;
    
    //return pdfGraph;
}


//______

void gaxis(){
    // draw TGaxis objects in various formats.
    //To see the output of this macro, click begin_html <a href="gif/gaxis.gif" >here</a> end_html
    
    TCanvas *c1 = new TCanvas("c1","Examples of Gaxis",10,10,700,500);
    
    c1->Range(-10,-1,10,1);
    
    TGaxis *axis1 = new TGaxis(-4.5,-0.2,5.5,-0.2,-6,8,510,"");
    axis1->SetName("axis1");
    axis1->Draw();
    
    TGaxis *axis2 = new TGaxis(-4.5,0.2,5.5,0.2,0.001,10000,510,"G");
    axis2->SetName("axis2");
    axis2->Draw();
    
    TGaxis *axis3 = new TGaxis(-9,-0.8,-9,0.8,-8,8,50510,"");
    axis3->SetName("axis3");
    axis3->SetTitle("axis3");
    axis3->SetTitleOffset(0.5);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(-7,-0.8,-7,0.8,1,10000,50510,"G");
    axis4->SetName("axis4");
    axis4->SetTitle("axis4");
    axis4->Draw();
    
    TGaxis *axis5 = new TGaxis(-4.5,-0.6,5.5,-0.6,1.2,1.32,80506,"-+");
    axis5->SetName("axis5");
    axis5->SetLabelSize(0.03);
    axis5->SetTextFont(72);
    
    axis5->Draw();
    
    TGaxis *axis6 = new TGaxis(-4.5,0.5,5.5,0.5,100,900,50510,"-");
    axis6->SetName("axis6");
    axis6->Draw();
    TGaxis *axis6a = new TGaxis(-5.5,0.85,5.5,0.85,0,4.3e-6,510,"");
    axis6a->SetName("axis6a");
    axis6a->Draw();
    
    TGaxis *axis7 = new TGaxis(8,-0.8,8,0.8,0,9000,50510,"+L");
    axis7->SetName("axis7");
    axis7->Draw();
    
    //one can make axis going top->bottom. However because of a long standing
    //problem, the two x values should not be equal
    TGaxis *axis8 = new TGaxis(6.5,0.8,6.499,-0.8,0,90,50510,"-");
    axis8->SetName("axis8");
    axis8->Draw();
    
    c1->Modified();
    c1->Update();
    
    c1->SaveAs("assi.pdf");
}


//_________________________
double chi2test_Plots(Float_t w, TH1F *hist1, TH1F *hist2, const char *process, const char *variable, const char* energy)
{
    // Note: The parameter w is used to produce the 2 pictures in
    // the TH1::Chi2Test method. The 1st picture is produced with
    // w=0 and the 2nd with w=17 (see TH1::Chi2Test() help).
    
    TH1F *h1 = (TH1F*)(hist1->Clone("Geant4"));
    TH1F *h2 = (TH1F*)(hist2->Clone("GeantV"));
    TH1F *h3 = (TH1F*)(h1->Clone("G4/GeantV"));
    //TH1F *h3 = (TH1F*)(h2->Clone("GeantV/G4"));
    
    
    
    Int_t binIndex;
    Int_t firstx=0,lastx=0;
    Double_t bincontent;
    lastx  = h3->GetXaxis()->GetNbins()+1;
    firstx = h3->GetXaxis()->GetFirst();
    int firstNonZeroBin;
    for (binIndex=firstx;binIndex<=lastx;binIndex++) {
        bincontent = h3->GetBinContent(binIndex);
        if(bincontent!=0)
        {
            firstNonZeroBin=binIndex;
            break;
        }
        
    }
    //h1->GetXaxis()->SetRange(firstNonZeroBin, lastx);
    //h2->GetXaxis()->SetRange(firstNonZeroBin, lastx);
    //h3->GetXaxis()->SetRange(firstNonZeroBin, lastx);
    
    const Int_t n= h1->GetSize()-2; //size: nBins+2
    std::cout<<"Number of entries h1: "<<h1->GetEntries()<<std::endl;
    std::cout<<"Number of entries h2: "<<h2->GetEntries()<<std::endl;
    
    //int firstNonZero=h1->FindFirstBinAbove();
    ///int lastNonZero= h1->FindLastBinAbove();
    
    
    //DEBUG
    
    //h1
    TAxis *xaxis_h1 = h1->GetXaxis();
    int first_h1 = h1->GetXaxis()->GetFirst();//first binIndex
    int last_h1  = h1->GetXaxis()->GetNbins();//last binIndex
    double bblow_h1 = xaxis_h1->GetBinLowEdge(first_h1); //Min value x-axis
    double bbup_h1 = xaxis_h1->GetBinUpEdge(last_h1);
    int firstNonZero_h1=h1->FindFirstBinAbove();
    int lastNonZero_h1= h1->FindLastBinAbove();
    
    //h2
    TAxis *xaxis_h2 = h2->GetXaxis();
    int first_h2 = h2->GetXaxis()->GetFirst();
    int last_h2  = h2->GetXaxis()->GetNbins();
    double bblow_h2 = xaxis_h2->GetBinLowEdge(first_h2);
    double bbup_h2 = xaxis_h2->GetBinUpEdge(last_h2);
    int firstNonZero_h2=h2->FindFirstBinAbove();
    int lastNonZero_h2= h2->FindLastBinAbove();
    
    //h3
    h3->Divide(h2);
    auto minch= new TGraphAsymmErrors();
    minch->Divide(h1,h2, "pois");
    
    TAxis *xaxis_h3 = h3->GetXaxis();
    int first_h3 = h3->GetXaxis()->GetFirst();
    int last_h3  = h3->GetXaxis()->GetNbins();
    double bblow_h3 = xaxis_h3->GetBinLowEdge(first_h3);
    double bbup_h3 = xaxis_h3->GetBinUpEdge(last_h3);
    int firstNonZero_h3=h3->FindFirstBinAbove();
    int lastNonZero_h3=h3->FindLastBinAbove();
    
    
    
    //hres
    auto hres = new TH1D("Nres","Normalized Residuals",n,bblow_h1,bbup_h1); //Problems with empty bins ?? giusto?
    std::vector<double> resLo(h1->GetNbinsX()); //resLo vettore di dimensione uguale al numero di bin di h1
    std::vector<double> res2;
    
    Double_t chi2=h1->Chi2Test(h2,"UU P",resLo.data()); //qui applico il Chi2Test e memorizzo in resLo
    
    for (int i = 0; i < n; ++i) {
        if (h1->GetBinContent(i+1) > 0 || h2->GetBinContent(i+1) > 0) //se uno dei due histo ha un contenuto >0 OR not AND!!!
        {
            hres->SetBinContent(i+1, resLo[i] ); //altrimenti salto un paio di indici
            hres->SetBinError(i+1,1); //commented out
            res2.push_back(resLo[i]); //probabilmente il problema è qui! era resLo[i+1]
        }
    }
    
    
    //NEW LORENZO-try
    /*std::vector<double> resLo(h1->GetNbinsX()+2);
     std::vector<double> res2;
     
     // include underflow-overflow in the test
     h1->GetXaxis()->SetRange(0, h1->GetNbinsX()+1); //important: include underflow and overflow bins
     
     Double_t chi2=h1->Chi2Test(h2,"UU P",resLo.data());
     //if(variable=="AngleOut2")
     //if(  !strcmp(variable,"AngleOut2"))
     //std::cout<<"Occhio numero due: chi2 "<<chi2<<"\n";
     
     auto hres = (TH1D*) h1->Clone();
     //auto hres = new TH1D("Nres","Normalized Residuals",n,bblow_h1,bbup_h1);
     hres->Reset();
     hres->SetName("hres");
     
     for (int i = 0; i < hres->GetNbinsX()+2; ++i) {
     // skip empty bins
     if (h1->GetBinContent(i)==0 && h2->GetBinContent(i)==0) continue;
     hres->SetBinContent(i,resLo[i]);
     hres->SetBinError(i,1);
     res2.push_back(resLo[i]);
     }
     */
    
    
#ifdef DEBUG
    
    std::cout<<"\n\n****** Debug ****\n\n";
    
    if (h1->GetDimension() != h2->GetDimension() )
        std::cout<<"\n******** ALERT! : Histograms of "<< variable<<" have different dimensions. H1: "<<h1->GetDimension()<<", H2: "<<h1->GetDimension()<<"\n";
    else
        std::cout<<"\n******** OK! : Histograms of "<< variable<<" have same dimensions. H1: "<<h1->GetDimension()<<", H2: "<<h1->GetDimension()<<"\n";
    if (h1->GetXaxis()->GetNbins()!= h2->GetXaxis()->GetNbins() )
        std::cout<<"\n******** ALERT! : Histograms of "<< variable<<" have different n.Bins. H1: "<<h1->GetXaxis()->GetNbins()<<", H2: "<<h2->GetXaxis()->GetNbins()<<"\n";
    else
        std::cout<<"\n******** OK! : Histograms of "<< variable<<" have same n.Bins. H1: "<<h1->GetXaxis()->GetNbins()<<", H2: "<<h2->GetXaxis()->GetNbins()<<"\n";
    if (firstNonZero_h1 != firstNonZero_h2)
        std::cout<<"\n******** ALERT! : Histograms of "<< variable<<"  have different first non-zero bin. H1: "<<firstNonZero_h1<<", H2: "<<firstNonZero_h2<<"\n";
    else
        std::cout<<"\n******** OK! : Histograms of "<< variable<<"  have same first non-zero bin. H1: "<<firstNonZero_h1<<", H2: "<<firstNonZero_h2<<"\n";
    if (lastNonZero_h1 != lastNonZero_h2)
        std::cout<<"\n******** ALERT!: Histograms of "<< variable<<" have different last non-zero bin. H1: "<<lastNonZero_h1<<", H2: "<<lastNonZero_h2<<"\n";
    else
        std::cout<<"\n******** OK!: Histograms of "<< variable<<" have same last non-zero bin. H1: "<<lastNonZero_h1<<", H2: "<<lastNonZero_h2<<"\n";
    if (h1->GetXaxis()->GetFirst()!= h2->GetXaxis()->GetFirst())
        std::cout<<"\n******** ALERT!: Histograms of "<< variable<<" have different first bin. H1: "<<h1->GetXaxis()->GetFirst()<<", H2: "<<h2->GetXaxis()->GetFirst()<<"\n";
    else
        std::cout<<"\n******** OK!: Histograms of "<< variable<<" have same first bin. H1: "<<h1->GetXaxis()->GetFirst()<<", H2: "<<h2->GetXaxis()->GetFirst()<<"\n";
    if (h1->GetXaxis()->GetLast()!= h2->GetXaxis()->GetLast())
        std::cout<<"\n******** ALERT!: Histograms of "<< variable<<" have different last bin. H1: "<<h1->GetXaxis()->GetLast()<<", H2: "<<h2->GetXaxis()->GetLast()<<"\n";
    else
        std::cout<<"\n******** OK!: Histograms of "<< variable<<" have same last bin. H1: "<<h1->GetXaxis()->GetLast()<<", H2: "<<h2->GetXaxis()->GetLast()<<"\n";
    
    std::cout<<"\n***** H1 *****\n";
    std::cout<<"h1: GetSize()="<<h1->GetSize()<<"\nh1->GetEntries()="<<h1->GetEntries()<<"\nFirst bin index: "<<first_h1<<" Content FirstBin: "<<h1->GetBinContent(first_h1)<<"\nLast bin index: "<<last_h1<<" Content LastBin: "<<h1->GetBinContent(last_h1)<<"\nMinimum value x-axis: "<<bblow_h1<<"\nMaximum value x-axis: "<<bbup_h1<<"\n";
    std::cout<<"\nFirst Non-zero bin index: "<<firstNonZero_h1<<" Content FirstNZBin: "<<h1->GetBinContent(firstNonZero_h1)<<"\nLast Non-zero bin index: "<<lastNonZero_h1<<" Content LastNZBin: "<<h1->GetBinContent(lastNonZero_h1)<<"\n\n";
    
    std::cout<<"\n***** H2 *****\n";
    std::cout<<"h2: GetSize()="<<h2->GetSize()<<"\nh2->GetEntries()="<<h2->GetEntries()<<"\nFirst bin index: "<<first_h2<<" Content FirstBin: "<<h2->GetBinContent(first_h2)<<"\nLast bin index: "<<last_h2<<" Content LastBin: "<<h2->GetBinContent(last_h2)<<"\nMinimum value x-axis: "<<bblow_h2<<"\nMaximum value x-axis: "<<bbup_h2<<"\n";
    std::cout<<"\nFirst Non-zero bin index: "<<firstNonZero_h2<<" Content FirstNZBin: "<<h2->GetBinContent(firstNonZero_h2)<<"\nLast Non-zero bin index: "<<lastNonZero_h2<<" Content LastNZBin: "<<h2->GetBinContent(lastNonZero_h2)<<"\n\n";
    
    std::cout<<"\n***** H3 *****\n";
    std::cout<<"h3: GetSize()="<<h3->GetSize()<<"\nh3->GetEntries()="<<h3->GetEntries()<<"\nFirst bin index: "<<first_h3<<" Content FirstBin: "<<h3->GetBinContent(first_h3)<<"\nLast bin index: "<<last_h3<<" Content LastBin: "<<h3->GetBinContent(last_h3)<<"\nMinimum value x-axis: "<<bblow_h3<<"\nMaximum value x-axis: "<<bbup_h3<<"\n";
    std::cout<<"\nFirst Non-zero bin index: "<<firstNonZero_h3<<" Content FirstNZBin: "<<h3->GetBinContent(firstNonZero_h3)<<"\nLast Non-zero bin index: "<<lastNonZero_h3<<" Content LastNZBin: "<<h3->GetBinContent(lastNonZero_h3)<<"\n\n";
    
#endif
    
    TAxis *xaxis_hres = hres->GetXaxis();
    int first_hres = hres->GetXaxis()->GetFirst();//primo binIndex
    int last_hres  = hres->GetXaxis()->GetNbins();//ultimo binIndex --> c'era un +1
    double bblow_hres = xaxis_hres->GetBinLowEdge(first_hres); //Edge è il valore nell'asse x
    double bbup_hres = xaxis_hres->GetBinUpEdge(last_hres);
    int firstNonZero_hres=hres->FindFirstBinAbove();
    int lastNonZero_hres= hres->FindLastBinAbove();
    
#ifdef DEBUG
    
    std::cout<<"***** HRES *****\n";
    std::cout<<"hRes: GetSize()="<<hres->GetSize()<<"\nhres->GetEntries()="<<hres->GetEntries()<<"\nFirst bin index: "<<first_hres<<" Content FirstBin: "<<hres->GetBinContent(first_hres)<<"\nLast bin index: "<<last_hres<<" Content LastBin: "<<hres->GetBinContent(last_hres)<<"\nMinimum value x-axis: "<<bblow_hres<<"\nMaximum value x-axis: "<<bbup_hres<<"\n";
    std::cout<<"\nFirst Non-zero bin index: "<<firstNonZero_hres<<" Content FirstNZBin: "<<hres->GetBinContent(firstNonZero_hres)<<"\nLast Non-zero bin index: "<<lastNonZero_hres<<" Content LastNZBin: "<<hres->GetBinContent(lastNonZero_hres)<<"\n\n";
    
    
    int binmin_hres = hres->GetMinimumBin();
    double binpositionMin_hres = hres->GetXaxis()->GetBinCenter(binmin_hres);
    int binmax_hres = hres->GetMaximumBin();
    double binpositionMax_hres = hres->GetXaxis()->GetBinCenter(binmax_hres);
    std::cout<<"Media hist: "<<hres->GetMean()<<"\nMinimum bin: "<<binmin_hres<<" at the x-value: "<<binpositionMin_hres<<" and at the y-value: "<<hres->GetBinContent(binmin_hres)<<"\nMaximum bin: "<<binmax_hres<<" at the x-value: "<<binpositionMax_hres<<" and at the y-value: "<<hres->GetBinContent(binmax_hres)<<"\n\n";
    
    
    if(lastNonZero_hres!=lastNonZero_h1)
        for(int index=1; index<(lastNonZero_h1-lastNonZero_hres)+1; index++)
            std::cout<<"Bin["<<index+lastNonZero_hres<<"]:"<<hres->GetBinContent(index+lastNonZero_hres)<<"\n";
    
    std::cout<<"**** HRES CONTENT ***\n";
    double media=0;
    double content;
    for (int i=1; i<=hres->GetNbinsX(); i++)
    {
        content= hres->GetBinContent(i);
        std::cout<<"hres["<<i<<"]= "<<content<<"\n";
        media+=content;
    }
    
    media=media/hres->GetSize();
    std::cout<<"MEAN: "<<media<<"\n";
    
#endif
    
    //***** Graph for Residuals *****//
    
    hres->GetYaxis()->SetTitle("Normalized Residuals");
    hres->SetMarkerStyle(21);
    hres->SetMarkerColor(6);
    hres->SetMarkerSize(.2);
    //hres->SetStats(kFALSE);
    
    //***** Quantile-Quantile plot: with the Gaussian *****//
    
    TF1 *f = new TF1("f","TMath::Gaus(x,0,1)",-10,10);
    TGraphQQ *qqplot = new TGraphQQ(res2.size(),res2.data(),f);
    //TGraphQQ *qqplot = new TGraphQQ(resLo.size(),resLo.data(),f);
    qqplot->SetMarkerStyle(20);
    qqplot->SetMarkerColor(6);
    qqplot->SetMarkerSize(.4);
    qqplot->SetTitle("Q-Q plot of Normalized Residuals");
    
    
    //Create Canvas
    
    TCanvas *canvasChi2Plots = new TCanvas("Chi-stat Plot","Chi-stat Plot",10,10,700,600);
    canvasChi2Plots->SetFillColor(0);
    canvasChi2Plots->Divide(2,2);
    
    
    //***** Draw Histogramms and Graphs *****//
    
    //cd 1
    canvasChi2Plots->cd(1);
    gPad->SetFrameFillColor(0);
    //gStyle->SetOptStat("kKsSiIrRmMen");
    h1->SetMarkerColor(4);
    h1->SetMarkerSize(.2);
    h1->SetMarkerStyle(20);
    h1->Draw("");
    h1->GetYaxis()->SetTitleOffset(1.5);
    
    TPaveStats *ps1=(TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.5); ps1->SetX2NDC(0.70);
    
    
    h2->Draw("SAME");
    ps1 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.72); ps1->SetX2NDC(0.96);
    gPad->SetFrameFillColor(0);
    drawLatex(.05,0.95,Form("%g",chi2));
    
    
    //cd 2
    canvasChi2Plots->cd(2);
    gPad->SetFrameFillColor(0);
    gPad->SetGrid();
    gPad->SetFrameBorderSize(12);
    
    
    /*h3->SetTitle("Histograms Ratio");
     h3->GetXaxis()->SetTitle("Bins");
     h3->GetYaxis()->SetTitle("G4/GV");
     h3->SetMarkerColor(6);
     h3->SetMarkerSize(.2);
     h3->SetMarkerStyle(20);
     //h3->Draw("");*/
    
    minch->SetTitle("Histograms Ratio");
    minch->GetXaxis()->SetTitle("Bins");
    minch->GetYaxis()->SetTitle("G4/GV");
    minch->GetYaxis()->SetTitleOffset(1.4);
    minch->SetMarkerColor(4);
    minch->SetMarkerSize(.6);
    minch->SetMarkerStyle(21);
    
    
    minch->Draw("ALP");
    ps1 = (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.70); ps1->SetX2NDC(0.96);
    
    //cd 3
    canvasChi2Plots->cd(3);
    gPad->SetFrameFillColor(0);
    gPad->SetGrid();
    hres->SetStats(0);
    //resgr->Draw("APL"); //invece di plottare questo:
    //hres->Draw("E");
    hres->SetMarkerStyle(20);
    hres->Fit("pol0");
    hres->Draw();
    
    
    //cd 4
    canvasChi2Plots->cd(4);
    gPad->SetFrameFillColor(0);
    qqplot->Draw("AP");
    canvasChi2Plots->Update();
    
    //si potrebbe evitare
    TString chiTestName=Form("%s-%s-%sMeV-ChiTest.pdf", process,variable,energy);
    canvasChi2Plots->SaveAs(chiTestName);
    
    TString nameFileMerge=Form("%s-Validation.pdf", process);
    canvasChi2Plots->Print(nameFileMerge);
    canvasChi2Plots->Close();
    return chi2;
    
}

// need to use Functor1D
double landau(double x) { return ROOT::Math::landau_pdf(x); }

void goftest() {
    
    // ------------------------------------------------------------------------
    // C a s e  1 :  C r e a t e   l o g N o r m a l  r a n d o m  s a m p l e
    // ------------------------------------------------------------------------
    
    UInt_t nEvents1 = 1000;
    
    //ROOT::Math::Random<ROOT::Math::GSLRngMT> r;
    TF1 * f1 = new TF1("logNormal","ROOT::Math::lognormal_pdf(x,[0],[1])",0,500);
    // set the lognormal parameters (m and s)
    f1->SetParameters(4.0,1.0);
    f1->SetNpx(1000);
    
    
    Double_t* sample1 = new Double_t[nEvents1];
    
    auto h1smp = new TH1D("h1smp", "LogNormal distribution histogram", 100, 0, 500);
    h1smp->SetStats(kFALSE);
    
    for (UInt_t i = 0; i < nEvents1; ++i) {
        //Double_t data = f1->GetRandom();
        Double_t data = gRandom->Gaus(4,1);
        data = TMath::Exp(data);
        sample1[i] = data;
        h1smp->Fill(data); //istogramma riempito con dati da Gaussiana
    }
    // normalize correctly the histogram using the entries inside
    h1smp->Scale( ROOT::Math::lognormal_cdf(500.,4.,1) / nEvents1, "width");
    
    TCanvas* c = new TCanvas("c","1-Sample and 2-Samples GoF Tests");
    c->Divide(1, 2);
    TPad * pad = (TPad *)c->cd(1);
    h1smp->Draw();
    h1smp->SetLineColor(kBlue);
    pad->SetLogy();
    f1->SetNpx(100); // use same points as histo for drawing
    f1->SetLineColor(kRed);
    f1->Draw("SAME");
    
    // -----------------------------------------
    // C r e a t e   G o F T e s t  o b j e c t
    // -----------------------------------------
    
    ROOT::Math::GoFTest* goftest_1 = new ROOT::Math::GoFTest(nEvents1, sample1, ROOT::Math::GoFTest::kLogNormal);
    
    /* Possible calls for the Anderson - DarlingTest test */
    /*----------------------------------------------------*/
    
    /* a) Returning the Anderson-Darling standardized test statistic */
    Double_t A2_1 = goftest_1-> AndersonDarlingTest("t");
    Double_t A2_2 = (*goftest_1)(ROOT::Math::GoFTest::kAD, "t");
    assert(A2_1 == A2_2);
    
    /* b) Returning the p-value for the Anderson-Darling test statistic */
    Double_t pvalueAD_1 = goftest_1-> AndersonDarlingTest(); // p-value is the default choice
    Double_t pvalueAD_2 = (*goftest_1)(); // p-value and Anderson - Darling Test are the default choices
    assert(pvalueAD_1 == pvalueAD_2);
    
    /* Rebuild the test using the default 1-sample construtor */
    delete goftest_1;
    goftest_1 = new ROOT::Math::GoFTest(nEvents1, sample1 ); // User must then input a distribution type option
    goftest_1->SetDistribution(ROOT::Math::GoFTest::kLogNormal);
    
    
    /* Possible calls for the Kolmogorov - Smirnov test */
    /*--------------------------------------------------*/
    
    /* a) Returning the Kolmogorov-Smirnov standardized test statistic */
    Double_t Dn_1 = goftest_1-> KolmogorovSmirnovTest("t");
    Double_t Dn_2 = (*goftest_1)(ROOT::Math::GoFTest::kKS, "t");
    assert(Dn_1 == Dn_2);
    
    /* b) Returning the p-value for the Kolmogorov-Smirnov test statistic */
    Double_t pvalueKS_1 = goftest_1-> KolmogorovSmirnovTest();
    Double_t pvalueKS_2 = (*goftest_1)(ROOT::Math::GoFTest::kKS);
    assert(pvalueKS_1 == pvalueKS_2);
    
    /* Valid but incorrect calls for the 2-samples methods of the 1-samples constucted goftest_1 */
#ifdef TEST_ERROR_MESSAGE
    Double_t A2 = (*goftest_1)(ROOT::Math::GoFTest::kAD2s, "t");     // Issues error message
    Double_t pvalueKS = (*goftest_1)(ROOT::Math::GoFTest::kKS2s);    // Issues error message
    assert(A2 == pvalueKS);
#endif
    
    TPaveText* pt1 = new TPaveText(0.58, 0.6, 0.88, 0.80, "brNDC");
    Char_t str1[50];
    sprintf(str1, "p-value for A-D 1-smp test: %f", pvalueAD_1); //AdamsonSmith
    pt1->AddText(str1);
    pt1->SetFillColor(18);
    pt1->SetTextFont(20);
    pt1->SetTextColor(4);
    Char_t str2[50];
    sprintf(str2, "p-value for K-S 1-smp test: %f", pvalueKS_1); //Kolmogorov-Smirnov
    pt1->AddText(str2);
    pt1->Draw();
    
    // ------------------------------------------------------------------------
    // C a s e  2 :  C r e a t e   G a u s s i a n  r a n d o m  s a m p l e s
    // ------------------------------------------------------------------------
    
    UInt_t nEvents2 = 2000;
    
    Double_t* sample2 = new Double_t[nEvents2];
    
    auto h2smps_1 = new TH1D("h2smps_1", "Gaussian distribution histograms", 100, 0, 500);
    h2smps_1->SetStats(kFALSE);
    
    auto h2smps_2 = new TH1D("h2smps_2", "Gaussian distribution histograms", 100, 0, 500);
    h2smps_2->SetStats(kFALSE);
    
    TRandom3 r;
    for (UInt_t i = 0; i < nEvents1; ++i) {
        Double_t data = r.Gaus(300, 50);
        sample1[i] = data;
        h2smps_1->Fill(data);
    }
    h2smps_1->Scale(1. / nEvents1, "width");
    c->cd(2);
    h2smps_1->Draw();
    h2smps_1->SetLineColor(kBlue);
    
    for (UInt_t i = 0; i < nEvents2; ++i) {
        Double_t data = r.Gaus(300, 50);
        sample2[i] = data;
        h2smps_2->Fill(data);
    }
    h2smps_2->Scale(1. / nEvents2, "width");
    h2smps_2->Draw("SAME");
    h2smps_2->SetLineColor(kRed);
    
    // -----------------------------------------
    // C r e a t e   G o F T e s t  o b j e c t
    // -----------------------------------------
    
    ROOT::Math::GoFTest* goftest_2 = new ROOT::Math::GoFTest(nEvents1, sample1, nEvents2, sample2);
    
    /* Possible calls for the Anderson - DarlingTest test */
    /*----------------------------------------------------*/
    
    /* a) Returning the Anderson-Darling standardized test statistic */
    A2_1 = goftest_2->AndersonDarling2SamplesTest("t");
    A2_2 = (*goftest_2)(ROOT::Math::GoFTest::kAD2s, "t");
    assert(A2_1 == A2_2);
    
    /* b) Returning the p-value for the Anderson-Darling test statistic */
    pvalueAD_1 = goftest_2-> AndersonDarling2SamplesTest(); // p-value is the default choice
    pvalueAD_2 = (*goftest_2)(ROOT::Math::GoFTest::kAD2s);  // p-value is the default choices
    assert(pvalueAD_1 == pvalueAD_2);
    
    /* Possible calls for the Kolmogorov - Smirnov test */
    /*--------------------------------------------------*/
    
    /* a) Returning the Kolmogorov-Smirnov standardized test statistic */
    Dn_1 = goftest_2-> KolmogorovSmirnov2SamplesTest("t");
    Dn_2 = (*goftest_2)(ROOT::Math::GoFTest::kKS2s, "t");
    assert(Dn_1 == Dn_2);
    
    /* b) Returning the p-value for the Kolmogorov-Smirnov test statistic */
    pvalueKS_1 = goftest_2-> KolmogorovSmirnov2SamplesTest();
    pvalueKS_2 = (*goftest_2)(ROOT::Math::GoFTest::kKS2s);
    assert(pvalueKS_1 == pvalueKS_2);
    
#ifdef TEST_ERROR_MESSAGE
    /* Valid but incorrect calls for the 1-sample methods of the 2-samples constucted goftest_2 */
    A2 = (*goftest_2)(ROOT::Math::GoFTest::kAD, "t");     // Issues error message
    pvalueKS = (*goftest_2)(ROOT::Math::GoFTest::kKS);    // Issues error message
    assert(A2 == pvalueKS);
#endif
    
    TPaveText* pt2 = new TPaveText(0.13, 0.6, 0.43, 0.8, "brNDC");
    sprintf(str1, "p-value for A-D 2-smps test: %f", pvalueAD_1);
    pt2->AddText(str1);
    pt2->SetFillColor(18);
    pt2->SetTextFont(20);
    pt2->SetTextColor(4);
    sprintf(str2, "p-value for K-S 2-smps test: %f", pvalueKS_1);
    pt2-> AddText(str2);
    pt2-> Draw();
    
    
    // ------------------------------------------------------------------------
    // C a s e  3 :  C r e a t e   L a n d a u  r a n d o m  s a m p l e
    // ------------------------------------------------------------------------
    
    UInt_t nEvents3 = 1000;
    
    Double_t* sample3 = new Double_t[nEvents3];
    for (UInt_t i = 0; i < nEvents3; ++i) {
        Double_t data = r.Landau();
        sample3[i] = data;
    }
    
    // ------------------------------------------
    // C r e a t e   G o F T e s t  o b j e c t s
    // ------------------------------------------
    
    /* Possible constructors for the user input distribution */
    /*-------------------------------------------------------*/
    
    /* a) User input PDF */
    ROOT::Math::Functor1D f(&landau);
    double min = 3*TMath::MinElement(nEvents3, sample3);
    double max = 3*TMath::MaxElement(nEvents3, sample3);
    ROOT::Math::GoFTest* goftest_3a = new ROOT::Math::GoFTest(nEvents3, sample3, f,  ROOT::Math::GoFTest::kPDF, min,max);  // need to specify am interval
    /* b) User input CDF */
    ROOT::Math::Functor1D fI(&TMath::LandauI);
    ROOT::Math::GoFTest* goftest_3b = new ROOT::Math::GoFTest(nEvents3, sample3, fI, ROOT::Math::GoFTest::kCDF,min,max);
    
    
    /* Returning the p-value for the Anderson-Darling test statistic */
    pvalueAD_1 = goftest_3a-> AndersonDarlingTest(); // p-value is the default choice
    
    pvalueAD_2 = (*goftest_3b)(); // p-value and Anderson - Darling Test are the default choices
    
    /* Checking consistency between both tests */
    std::cout << " \n\nTEST with LANDAU distribution:\t";
    if (TMath::Abs(pvalueAD_1 - pvalueAD_2) > 1.E-1 * pvalueAD_2) {
        std::cout << "FAILED " << std::endl;
        Error("goftest","Error in comparing testing using Landau and Landau CDF");
        std::cerr << " pvalues are " << pvalueAD_1 << "  " << pvalueAD_2 << std::endl;
    }
    else
        std::cout << "OK ( pvalues = " << pvalueAD_2 << "  )" << std::endl;
    
    c->SaveAs("GoF_test.pdf");
    //std::cout<<"Ho finito, ma ho disegnato qualcosa?\n";
}

void goftestSingleMod(const char *phisicsModel, const char*  variable, const char* energy){
    
    TString histPath= Form("%s/%s",phisicsModel, variable);
    TString histName= Form("%s/%s/%sMeV",phisicsModel, variable, energy);
    TString histoFileNamePdf= Form("%s-%s-%sMeV.pdf",phisicsModel, variable, energy);
    TString g4HistoName=Form("%s/geant4",phisicsModel);
    TString gvHistoScalarName=Form("%s/geantVscalar",phisicsModel);
    
    TString geant4RootFileName= Form("geant4_%sMeV.root",energy);
    TString scalarRootFileName= Form("scalar_%sMeV.root",energy);
    
    //Geant 4
    TFile *fG4 = new TFile(geant4RootFileName,"w");
    TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
    if(eOutG4) eOutG4->SetName("geant4");
    else std::cout<<"eOutG4 is null\n";
    
    //GeantV scalar
    TFile *fScalar = new TFile(scalarRootFileName,"w");
    TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
    if(eOutScalar) eOutScalar->SetName("geantVscalar");
    else std::cout<<"eOutScalar is null\n";
    
    //GeantV vector
    //TFile *fVector = new TFile(vectorRootFileName,"w");
    //TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
    //if(eOutVector) eOutVector->SetName("geantVvector");
    //else std::cout<<"eOutVector is null\n";
    
    TCanvas* c = new TCanvas("c","GoF Tests");
    Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    
    Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
    Double_t xScalar[entries], yScalar[entries], zScalar[entries];
    //Double_t xVector[entries], yVector[entries], zVector[entries];
    
    //Scale the histograms
    double norm=1.;
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    //eOutGeant4Scaled->Scale(1/norm);
    
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    //eOutScalarScaled->Scale(1/norm);
    
    //TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    //norm= eOutVectorScaled->GetEntries();
    //eOutVectorScaled->Scale(1/norm);
    
    for(int j = 0; j < entries ; ++j){
        yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
        yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
        //yVector[j] = eOutVectorScaled->GetBinContent(j+1);
    }
    
    UInt_t nEvents2 = entries;
    
    Double_t* sample2 = new Double_t[nEvents2];
    
    //TH1D* h2smps_1 = new TH1D("h2smps_1", "Gaussian distribution histograms", 100, 0, 500);
    //h2smps_1->SetStats(kFALSE);
    eOutGeant4Scaled->SetStats(kFALSE);
    eOutGeant4Scaled->Draw();
    eOutGeant4Scaled->SetLineColor(kBlue);
    
    //
    eOutScalarScaled->SetStats(kFALSE);
    eOutScalarScaled->Draw("SAME");
    eOutScalarScaled->SetLineColor(kRed);
    
    // -----------------------------------------
    // C r e a t e   G o F T e s t  o b j e c t
    // -----------------------------------------
    
    ROOT::Math::GoFTest* goftest_2 = new ROOT::Math::GoFTest(nEvents2, yGeant4, nEvents2, yScalar);
    
    /* Possible calls for the Anderson - DarlingTest test */
    /*----------------------------------------------------*/
    
    
    /* a) Returning the Anderson-Darling standardized test statistic */
    Double_t A2_1 = goftest_2->AndersonDarling2SamplesTest("t");
    Double_t A2_2 = (*goftest_2)(ROOT::Math::GoFTest::kAD2s, "t");
    assert(A2_1 == A2_2);
    
    /* b) Returning the p-value for the Anderson-Darling test statistic */
    Double_t pvalueAD_1 = goftest_2-> AndersonDarling2SamplesTest(); // p-value is the default choice
    Double_t pvalueAD_2 = (*goftest_2)(ROOT::Math::GoFTest::kAD2s);  // p-value is the default choices
    assert(pvalueAD_1 == pvalueAD_2);
    
    /* Possible calls for the Kolmogorov - Smirnov test */
    /*--------------------------------------------------*/
    
    /* a) Returning the Kolmogorov-Smirnov standardized test statistic */
    Double_t Dn_1 = goftest_2-> KolmogorovSmirnov2SamplesTest("t");
    Double_t Dn_2 = (*goftest_2)(ROOT::Math::GoFTest::kKS2s, "t");
    assert(Dn_1 == Dn_2);
    
    /* b) Returning the p-value for the Kolmogorov-Smirnov test statistic */
    Double_t pvalueKS_1 = goftest_2-> KolmogorovSmirnov2SamplesTest();
    Double_t pvalueKS_2 = (*goftest_2)(ROOT::Math::GoFTest::kKS2s);
    assert(pvalueKS_1 == pvalueKS_2);
    
#ifdef TEST_ERROR_MESSAGE
    /* Valid but incorrect calls for the 1-sample methods of the 2-samples constucted goftest_2 */
    Double_t A2 = (*goftest_2)(ROOT::Math::GoFTest::kAD, "t");     // Issues error message
    Double_t pvalueKS = (*goftest_2)(ROOT::Math::GoFTest::kKS);    // Issues error message
    assert(A2 == pvalueKS);
#endif
    
    TPaveText* pt2 = new TPaveText(0.13, 0.6, 0.43, 0.8, "brNDC");
    Char_t str1[50];
    sprintf(str1, "p-value for A-D 2-smps test: %f", pvalueAD_1);
    pt2->AddText(str1);
    pt2->SetFillColor(18);
    pt2->SetTextFont(20);
    pt2->SetTextColor(4);
    Char_t str2[50];
    sprintf(str2, "p-value for K-S 2-smps test: %f", pvalueKS_1);
    pt2-> AddText(str2);
    pt2-> Draw();
    
    c->Update();
    c->SaveAs("nientepopodimenoche.pdf");
    
}

//void goftestAllMod(int entries1, Double_t* sample1, int entries2, Double_t* sample2, Double_t pvalueAD_1, Double_t pvalueKS_1){
void goftestAllMod(TH1F *histo1, TH1F *histo2, Double_t & pvalueAD_1, Double_t & pvalueKS_1){
    
    UInt_t nEvents1 = histo1->GetSize()-2;
    UInt_t nEvents2 = histo2->GetSize()-2;
    std::cout<<"nEvents1: "<<nEvents1<<" nEvents2: "<<nEvents2<<"\n";
    
    // -----------------------------------------
    // Create the samples vectors
    // -----------------------------------------
    Double_t sample1[nEvents1], sample2[nEvents2];
    
    
    //If necessary scale the histograms
    //Double_t norm=histo1->GetEntries();
    //TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    //norm= eOutGeant4Scaled->GetEntries();
    //eOutGeant4Scaled->Scale(1/norm);
    
    //TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    //norm= eOutScalarScaled->GetEntries();
    //eOutScalarScaled->Scale(1/norm);
    
    //TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    //norm= eOutVectorScaled->GetEntries();
    //eOutVectorScaled->Scale(1/norm);
    
    for(int i = 0; i < nEvents1 ; ++i){
        sample1[i] = histo1->GetBinContent(i+1);
        std::cout<<sample1[i]<<" ";
    }
    std::cout<<"\n";
    for(int i = 0; i < nEvents2 ; ++i)
        sample2[i] = histo2->GetBinContent(i+1);
    
    
    // -----------------------------------------
    // C r e a t e   G o F T e s t  o b j e c t
    // -----------------------------------------
    
    ROOT::Math::GoFTest* goftest_2 = new ROOT::Math::GoFTest(nEvents1, sample1, nEvents2, sample2);
    
    /* Possible calls for the Anderson - DarlingTest test */
    /*----------------------------------------------------*/
    
    
    /* a) Returning the Anderson-Darling standardized test statistic */
    Double_t A2_1 = goftest_2->AndersonDarling2SamplesTest("t");
    Double_t A2_2 = (*goftest_2)(ROOT::Math::GoFTest::kAD2s, "t");
    assert(A2_1 == A2_2);
    
    /* b) Returning the p-value for the Anderson-Darling test statistic */
    pvalueAD_1 = goftest_2-> AndersonDarling2SamplesTest(); // p-value is the default choice
    Double_t pvalueAD_2 = (*goftest_2)(ROOT::Math::GoFTest::kAD2s);  // p-value is the default choices
    assert(pvalueAD_1 == pvalueAD_2);
    
    /* Possible calls for the Kolmogorov - Smirnov test */
    /*--------------------------------------------------*/
    
    /* a) Returning the Kolmogorov-Smirnov standardized test statistic */
    Double_t Dn_1 = goftest_2-> KolmogorovSmirnov2SamplesTest("t");
    Double_t Dn_2 = (*goftest_2)(ROOT::Math::GoFTest::kKS2s, "t");
    assert(Dn_1 == Dn_2);
    
    /* b) Returning the p-value for the Kolmogorov-Smirnov test statistic */
    pvalueKS_1 = goftest_2-> KolmogorovSmirnov2SamplesTest();
    Double_t pvalueKS_2 = (*goftest_2)(ROOT::Math::GoFTest::kKS2s);
    assert(pvalueKS_1 == pvalueKS_2);
    
#ifdef TEST_ERROR_MESSAGE
    /* Valid but incorrect calls for the 1-sample methods of the 2-samples constucted goftest_2 */
    Double_t A2 = (*goftest_2)(ROOT::Math::GoFTest::kAD, "t");     // Issues error message
    Double_t pvalueKS = (*goftest_2)(ROOT::Math::GoFTest::kKS);    // Issues error message
    assert(A2 == pvalueKS);
#endif
}


//______________________________________________________________________________
double chiSquare(TH1F* eOutG4, TH1F* eOutScalar)
{
    
    int entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    double chiSquare =0. ,chiSquare_scaled =0. , exp=0., obs=0. ;
    
    //calculate the chi-square of the scaled histo
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);
    
    
    for (int i=0; i<entries; i++)
    {
        exp = eOutG4->GetBinContent(i); //to verify
        obs = eOutScalar->GetBinContent(i);
        
        //std::cout<<"exp:  "<<exp<<" obs:  "<<obs<<std::endl;
        if(exp!=0)
            chiSquare+=((obs-exp)*(obs-exp)/exp);
        
        exp = eOutGeant4Scaled->GetBinContent(i); //to verify
        obs = eOutScalarScaled->GetBinContent(i);
        if(exp!=0)
            chiSquare_scaled+=((obs-exp)*(obs-exp)/exp);
        
    }
    std::cout<<"chiSquare "<<chiSquare<<std::endl;
    std::cout<<"chiSquare_scaled "<<chiSquare_scaled<<std::endl;
    
    return chiSquare_scaled;
}

//______________________________________________________________________________
void genAllHisto(const char *process, const char* energy, double *chi2PV, double *adPV, double *ksPV)
{
    
    TString processName= process;
    //TString histogramName= variable;
    TString geant4RootFileName= Form("geant4.root");
    TString scalarRootFileName= Form("scalar.root");
    TString vectorRootFileName= Form("vector.root");
    
    TString g4HistoName=Form("%s/geant4",process);
    TString gvHistoScalarName=Form("%s/geantVscalar",process);
    TString gvHistoVectorName=Form("%s/geantVvector",process);
    
    TString histPath= Form("%s/%s",process, "EnergyOut1");
    TString histName= Form("%s/%s/%sMeV",process, "EnergyOut1", energy);
    //TString histoFileNamePdf= Form("%s-%s-%sMeV.pdf",process, variable, energy);
    TString histoFileNameEps= Form("%s-%sMeV.eps",process, energy);
    TString histoFileNamePdf= Form("%s-%sMeV.pdf",process, energy);
    
    //Geant 4
    TFile *fG4 = new TFile(geant4RootFileName,"w");
    TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
    if(eOutG4) eOutG4->SetName(" geant4 ");
    else std::cout<<"eOutG4 is null\n";
    
    //GeantV scalar
    TFile *fScalar = new TFile(scalarRootFileName,"w");
    TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
    if(eOutScalar) eOutScalar->SetName("geantVscalar");
    else std::cout<<"eOutScalar is null\n";
    
    //GeantV vector
    TFile *fVector = new TFile(vectorRootFileName,"w");
    TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
    if(eOutVector) eOutVector->SetName("geantVvector");
    else std::cout<<"eOutVector is null\n";
    
    TCanvas *firstPage=new TCanvas(process,process,800,200,1200,1000);
    drawText(.4,0.45,Form("Model: %s, energyIn: %s MeV ",process, energy));
    
    TString nameFileMerge=Form("%s-Validation.pdf", process);
    firstPage->Print(nameFileMerge);
    
    
    TCanvas *MyC1 = new TCanvas(process,process,800,200,1200,1000);
    MyC1->Divide(2,2);
    
    
    //1
    MyC1->cd(1);
    gPad->SetLogy();
    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Scattered Photon Energy [MeV]");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("EnergyOut1");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    
    gStyle->SetEndErrorSize(3);
    gStyle->SetErrorX(1.);
    eOutG4->Draw("E");
    
    
    
    MyC1->Update();
    
    TPaveStats *ps1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.55);
    MyC1->Modified();
    MyC1->Update();
    eOutScalar->SetLineColor(kYellow+10);
    eOutScalar->Draw("sames");
    
    
    
    MyC1->Update();
    
    TPaveStats *ps2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.6); ps2->SetX2NDC(0.75);
    ps2->SetTextColor(kYellow+10);
    ps2->SetOptStat(1111);
    
    eOutVector->SetLineColor(kMagenta);
    //eOutVector->Draw("][sames");
    eOutVector->Draw("sames");
    
    
    
    MyC1->Update();
    TPaveStats *ps3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.8); ps3->SetX2NDC(0.95);
    ps3->SetTextColor(kMagenta);
    TPaveText *t = new TPaveText(0.0, 0.95, 0.3, 1.0, "brNDC"); // left-up
    //t->AddText(histName);
    
    //double res[eOutG4->GetSize()];
    double res[eOutG4->GetNbinsX()+2];
    
    //Added for uoFlows_Analysis
    // include underflow-overflow in the test
    //eOutG4->GetXaxis()->SetRange(0, eOutG4->GetNbinsX()+1);
    double pValueP=eOutG4->Chi2Test(eOutScalar,"UU P",res);
    //std::cout<<"Occhio qui E1: pValueP "<<pValueP<<"\n";
    //drawText(.1,0.95,Form("chiQuadro: %f",chi2 ));
    
    //drawLatex(.1,0.95,Form("%g",chi2));
    
    
    /*AGGIUNTO PER AD e KS test*/
    Double_t pvalueAD_1, pvalueKS_1;
    goftestAllMod(eOutG4, eOutScalar, pvalueAD_1, pvalueKS_1);
    
    TPaveText* pt1 = new TPaveText(0.13, 0.55, 0.43, 0.75, "brNDC");
    Char_t str[50];
    sprintf(str, "p-value for Pearson chi2 test: %g", pValueP);
    Char_t str1[50];
    sprintf(str1, "p-value for A-D 2-smps test: %f", pvalueAD_1);
    pt1->AddText(str);
    pt1->AddText(str1);
    pt1->SetFillColor(18);
    pt1->SetTextFont(20);
    pt1->SetTextColor(4);
    Char_t str2[50];
    sprintf(str2, "p-value for K-S 2-smps test: %f", pvalueKS_1);
    pt1-> AddText(str2);
    pt1-> Draw();
    
    chi2PV[0]=pValueP;
    adPV[0]=pvalueAD_1;
    ksPV[0]=pvalueKS_1;
    std::cout<<"EnergyOut1 -  \nPearson p-Value: "<<pValueP<<" \nA-D p-Value: "<<pvalueAD_1<<" \nK-S p-Value: "<<pvalueKS_1<<"\n ";
    
    MyC1->Modified();
    MyC1->Update();
    
    chi2test_Plots(0, eOutG4, eOutScalar, process, "EnergyOut1", energy);
    
    
    //2
    MyC1->cd(2);
    gPad->SetLogy();
    histPath= Form("%s/%s",process, "EnergyOut2");
    histName= Form("%s/%s/%sMeV",process, "EnergyOut2", energy);
    eOutG4 = (TH1F*)fG4->Get(histPath);
    eOutScalar = (TH1F*)fScalar->Get(histPath);
    eOutVector = (TH1F*)fVector->Get(histPath);
    eOutG4->SetName(" geant4 ");
    eOutScalar->SetName("geantVscalar");
    eOutVector->SetName("geantVvector");
    eOutScalar->SetLineColor(kYellow+10);
    eOutVector->SetLineColor(kMagenta);
    
    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Electron Energy [MeV]");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("EnergyOut2");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    
    
    eOutG4->Draw("E");
    eOutScalar->Draw("sames");
    eOutVector->Draw("sames");
    MyC1->Update();
    TPaveStats *psEnOut1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    psEnOut1->SetX1NDC(0.4); psEnOut1->SetX2NDC(0.55);
    //psEnOut1->SetTextColor(kMagenta);
    
    TPaveStats *psEnOut2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    psEnOut2->SetX1NDC(0.6); psEnOut2->SetX2NDC(0.75);
    psEnOut2->SetTextColor(kYellow+10);
    psEnOut2->SetOptStat(1111);
    
    //chiS=chiSquare(eOutG4,eOutScalar);
    //chiSquareString= Form("%f",chiS);
    //drawtext(.1, 0.92, chiSquareString);
    
    TPaveStats *psEnOut3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    psEnOut3->SetX1NDC(0.8); psEnOut3->SetX2NDC(0.95);
    psEnOut3->SetTextColor(kMagenta);
    
    //double mio=chi2test(0, eOutG4, eOutScalar, process, "EnergyOut2", energy);
    
    TText *t22 = new TText(0.05,0.8,"This is pad22");
    t22->SetTextSize(10);
    t22->Draw("sames");
    MyC1->Modified();
    MyC1->Update();
    
    
    double resE2[eOutG4->GetNbinsX()+2];
    //Added for uoFlows_Analysis
    // include underflow-overflow in the test
    //eOutG4->GetXaxis()->SetRange(0, eOutG4->GetNbinsX()+1);
    pValueP=eOutG4->Chi2Test(eOutScalar,"UU P",resE2);
    //std::cout<<"Occhio qui E2: pValueP "<<pValueP<<"\n";
    
    
    /*AGGIUNTO PER AD e KS test*/
    //Double_t pvalueAD_1, pvalueKS_1;
    goftestAllMod(eOutG4, eOutScalar, pvalueAD_1, pvalueKS_1);
    
    TPaveText* pt2 = new TPaveText(0.13, 0.55, 0.43, 0.75, "brNDC");
    //Char_t str[50];
    sprintf(str, "p-value for Pearson chi2 test: %g", pValueP);
    //Char_t str1[50];
    sprintf(str1, "p-value for A-D 2-smps test: %f", pvalueAD_1);
    pt2->AddText(str);
    pt2->AddText(str1);
    pt2->SetFillColor(18);
    pt2->SetTextFont(20);
    pt2->SetTextColor(4);
    //Char_t str2[50];
    sprintf(str2, "p-value for K-S 2-smps test: %f", pvalueKS_1);
    pt2-> AddText(str2);
    pt2-> Draw();
    std::cout<<"EnergyOut2 -  \nPearson p-Value: "<<pValueP<<" \nA-D p-Value: "<<pvalueAD_1<<" \nK-S p-Value: "<<pvalueKS_1<<"\n ";
    
    
    chi2PV[1]=pValueP;
    adPV[1]=pvalueAD_1;
    ksPV[1]=pvalueKS_1;
    /**/
    
    
    //drawText(.1,0.95,Form("chiQuadro: %f",chi2 ));
    //drawLatex(.1,0.95,Form("%g",chi2));
    
    chi2test_Plots(0, eOutG4, eOutScalar, process, "EnergyOut2", energy);
    
    
    //3
    MyC1->cd(3);
    gPad->SetLogy();
    histPath= Form("%s/%s",process, "AngleOut1");
    histName= Form("%s/%s/%sMeV",process, "AngleOut1", energy);
    eOutG4 = (TH1F*)fG4->Get(histPath);
    eOutScalar = (TH1F*)fScalar->Get(histPath);
    eOutVector = (TH1F*)fVector->Get(histPath);
    eOutG4->SetName(" geant4 ");
    eOutScalar->SetName("geantVscalar");
    eOutVector->SetName("geantVvector");
    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Scattered Photon Angle");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("AngleOut1");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    
    //eOutG4->GetXaxis()->SetTitle("AngleOut1");
    
    eOutG4->Draw("E");
    eOutScalar->SetLineColor(kYellow+10);
    eOutVector->SetLineColor(kMagenta);
    eOutScalar->Draw("sames");
    eOutVector->Draw("sames");
    MyC1->Update();
    
    //chiS=chiSquare(eOutG4,eOutScalar);
    //chiSquareString= Form("%f",chiS);
    //drawLatex(.1, 0.92, chiSquareString);
    
    double resA1[eOutG4->GetNbinsX()+2];
    //Added for uoFlows_Analysis
    // include underflow-overflow in the test
    //eOutG4->GetXaxis()->SetRange(0, eOutG4->GetNbinsX()+1);
    pValueP=eOutG4->Chi2Test(eOutScalar,"UU P",resA1);
    //std::cout<<"Occhio qui A1: pValueP "<<pValueP<<"\n";
    
    /*AGGIUNTO PER AD e KS test*/
    //Double_t pvalueAD_1, pvalueKS_1;
    goftestAllMod(eOutG4, eOutScalar, pvalueAD_1, pvalueKS_1);
    
    TPaveText* pt3 = new TPaveText(0.13, 0.55, 0.43, 0.75, "brNDC");
    //Char_t str[50];
    sprintf(str, "p-value for Pearson chi2 test: %g", pValueP);
    //Char_t str1[50];
    sprintf(str1, "p-value for A-D 2-smps test: %f", pvalueAD_1);
    pt3->AddText(str);
    pt3->AddText(str1);
    pt3->SetFillColor(18);
    pt3->SetTextFont(20);
    pt3->SetTextColor(4);
    //Char_t str2[50];
    sprintf(str2, "p-value for K-S 2-smps test: %f", pvalueKS_1);
    pt3-> AddText(str2);
    pt3-> Draw();
    
    std::cout<<"AngleOut1 -  \nPearson p-Value: "<<pValueP<<" \nA-D p-Value: "<<pvalueAD_1<<" \nK-S p-Value: "<<pvalueKS_1<<"\n ";
    
    chi2PV[2]=pValueP;
    adPV[2]=pvalueAD_1;
    ksPV[2]=pvalueKS_1;
    /**/
    
    //drawText(.1,0.95,Form("chiQuadro: %f",chi2 ));
    //drawLatex(.1,0.95,Form("%g",chi2));
    
    TPaveStats *psAOut1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    psAOut1->SetX1NDC(0.4); psAOut1->SetX2NDC(0.55);
    
    TPaveStats *psAOut2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    psAOut2->SetX1NDC(0.6); psAOut2->SetX2NDC(0.75);
    psAOut2->SetTextColor(kYellow+10);
    psAOut2->SetOptStat(1111);
    
    TPaveStats *psAOut3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    psAOut3->SetX1NDC(0.8); psAOut3->SetX2NDC(0.95);
    psAOut3->SetTextColor(kMagenta);
    chi2test_Plots(0, eOutG4, eOutScalar, process,"AngleOut1", energy);
    
    //4
    MyC1->cd(4);
    gPad->SetLogy();
    histPath= Form("%s/%s",process, "AngleOut2");
    histName= Form("%s/%s/%sMeV",process, "AngleOut2", energy);
    eOutG4 = (TH1F*)fG4->Get(histPath);
    eOutScalar = (TH1F*)fScalar->Get(histPath);
    eOutVector = (TH1F*)fVector->Get(histPath);
    eOutG4->SetName("geant4");
    eOutScalar->SetName("geantVscalar");
    eOutVector->SetName("geantVvector");
    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Electron Angle");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("AngleOut2");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    //eOutG4->GetXaxis()->SetTitle("AngleOut2");
    eOutG4->Draw("E");
    eOutScalar->SetLineColor(kYellow+10);
    eOutVector->SetLineColor(kMagenta);
    eOutScalar->Draw("sames");
    eOutVector->Draw("sames");
    MyC1->Update();
    
    //chiS=chiSquare(eOutG4,eOutScalar);
    //chiSquareString= Form("%f",chiS);
    //drawLatex(.1, 0.92, chiSquareString);
    
    
    double resA2[eOutG4->GetNbinsX()+2];
    //Added for uoFlows_Analysis
    // include underflow-overflow in the test
    //eOutG4->GetXaxis()->SetRange(0, eOutG4->GetNbinsX()+1);
    pValueP=eOutG4->Chi2Test(eOutScalar,"UU P",resA2);
    //std::cout<<"Occhio qui A2: pValueP "<<pValueP<<"\n";
    
    /*AGGIUNTO PER AD e KS test*/
    //Double_t pvalueAD_1, pvalueKS_1;
    goftestAllMod(eOutG4, eOutScalar, pvalueAD_1, pvalueKS_1);
    
    TPaveText* pt4 = new TPaveText(0.13, 0.55, 0.43, 0.75, "brNDC");
    //Char_t str[50];
    sprintf(str, "p-value for Pearson chi2 test: %g", pValueP);
    //Char_t str1[50];
    sprintf(str1, "p-value for A-D 2-smps test: %f", pvalueAD_1);
    pt4->AddText(str);
    pt4->AddText(str1);
    pt4->SetFillColor(18);
    pt4->SetTextFont(20);
    pt4->SetTextColor(4);
    //Char_t str2[50];
    sprintf(str2, "p-value for K-S 2-smps test: %f", pvalueKS_1);
    pt4-> AddText(str2);
    pt4-> Draw();
    std::cout<<"AngleOut2 -  \nPearson p-Value: "<<pValueP<<" \nA-D p-Value: "<<pvalueAD_1<<" \nK-S p-Value: "<<pvalueKS_1<<"\n ";
    
    
    chi2PV[3]=pValueP;
    adPV[3]=pvalueAD_1;
    ksPV[3]=pvalueKS_1;
    /**/
    
    
    //drawText(.1,0.95,Form("chiQuadro: %f",chi2 ));
    //drawLatex(.1,0.95,Form("%g",chi2));
    
    
    
    TPaveStats *psAngleOut1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    psAngleOut1->SetX1NDC(0.4); psAngleOut1->SetX2NDC(0.55);
    
    TPaveStats *psAngleOut2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    psAngleOut2->SetX1NDC(0.6); psAngleOut2->SetX2NDC(0.75);
    psAngleOut2->SetTextColor(kYellow+10);
    psAngleOut2->SetOptStat(1111);
    
    TPaveStats *psAngleOut3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    psAngleOut3->SetX1NDC(0.8); psAngleOut3->SetX2NDC(0.95);
    psAngleOut3->SetTextColor(kMagenta);
    
    chi2test_Plots(0, eOutG4, eOutScalar, process, "AngleOut2", energy);
    MyC1->SaveAs(histoFileNameEps);
    MyC1->Print(nameFileMerge);
}



//______________________________________________________________________________
void genSingleHisto(const char *process, const char *variable, const char* energy)
{
    TString processName= process;
    TString histogramName= variable;
    
    TString histPath= Form("%s/%s",process, variable);
    TString histName= Form("%s/%s/%sMeV",process, variable, energy);
    TString histoFileNamePdf= Form("%s-%s-%sMeV.pdf",process, variable, energy);
    TString histoFileNameEps= Form("%s-%s-%sMeV.eps",process, variable, energy);
    
    TString g4HistoName=Form("%s/geant4",process);
    TString gvHistoScalarName=Form("%s/geantVscalar",process);
    TString gvHistoVectorName=Form("%s/geantVvector",process);
    
    
    TString geant4RootFileName= Form("geant4_%sMeV.root",energy);
    TString scalarRootFileName= Form("scalar_%sMeV.root",energy);
    TString vectorRootFileName= Form("vector_%sMeV.root",energy);
    //TString pdfFileName=Form("pdf_%sMeV.root",energy);
    
    //Geant 4
    TFile *fG4 = new TFile(geant4RootFileName,"w");
    //fG4->ls();
    //std::cout<<"Path to the histogram:"<< histPath<<std::endl;
    TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
    if(eOutG4) eOutG4->SetName("geant4");
    else std::cout<<"eOutG4 is null\n";
    
    //GeantV scalar
    TFile *fScalar = new TFile(scalarRootFileName,"w");
    TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
    if(eOutScalar) eOutScalar->SetName("geantVscalar");
    else std::cout<<"eOutScalar is null\n";
    
    //GeantV vector
    TFile *fVector = new TFile(vectorRootFileName,"w");
    TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
    if(eOutVector) eOutVector->SetName("geantVvector");
    else std::cout<<"eOutVector is null\n";
    
    //int intenergy= atoi(energy);
    //TGraph my=chiSquare_pdf(eOutG4, eOutScalar, eOutVector, intenergy);
    double chiS= chiSquare(eOutG4, eOutScalar);
    
    std::cout<<"chiS: "<<chiS<<"\n";
    
    TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
    c1->SetLogy();
    
    
    if(!strcmp(process,"KleinNishina"))
    {
        std::cout<<"esatto\n";
        if(!strcmp(variable,"EnergyIn"))
        {
            eOutG4->GetXaxis()->SetTitle("EnergyIn [MeV]");
            eOutG4->GetYaxis()->SetTitle("dN/dE");
            
        }else
            if(!strcmp(variable,"EnergyOut1"))
            {
                eOutG4->GetXaxis()->SetTitle("Scattered Photon Energy [MeV]");
                eOutG4->GetYaxis()->SetTitle("dN/dE");
                
            }else
                if(!strcmp(variable,"EnergyOut2"))
                {
                    eOutG4->GetXaxis()->SetTitle("Electron Energy [MeV]");
                    eOutG4->GetYaxis()->SetTitle("dN/dE");
                    
                }else
                    //angle of the scatterred photon
                    if(!strcmp(variable,"AngleOut1"))
                    {
                        eOutG4->GetXaxis()->SetTitle("Scattered Photon Angle");
                        //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                        
                    }else
                        if(!strcmp(variable,"AngleOut2"))
                        {
                            eOutG4->GetXaxis()->SetTitle("Electron Angle");
                            //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                        }
    } else
    {
        std::cout<<"sbagliato\n";
        if(!strcmp(variable,"EnergyIn"))
        {
            eOutG4->GetXaxis()->SetTitle("EnergyIn [MeV]");
            eOutG4->GetYaxis()->SetTitle("dN/dE");
            
        }else
            if(!strcmp(variable,"EnergyOut1"))
            {
                eOutG4->GetXaxis()->SetTitle("EnergyOut1 [MeV]");
                eOutG4->GetYaxis()->SetTitle("dN/dE");
                
            }else
                if(!strcmp(variable,"EnergyOut2"))
                {
                    eOutG4->GetXaxis()->SetTitle("EnergyOut2 [MeV]");
                    eOutG4->GetYaxis()->SetTitle("dN/dE");
                    
                }else
                    
                    if(!strcmp(variable,"AngleOut1"))
                    {
                        eOutG4->GetXaxis()->SetTitle("AngleOut1");
                        //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                        
                    }else
                        if(!strcmp(variable,"AngleOut2"))
                        {
                            eOutG4->GetXaxis()->SetTitle("AngleOut2");
                            //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                        }
        
        
    }
    
    eOutG4->Draw("E");
    //eOutG4->SetLineColor(kGreen);
    c1->Update();
    
    
    TPaveStats *ps1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    //ps1->SetTextSize(0.005);
    ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.55);
    //ps1->SetX1NDC(0.6); ps1->SetX2NDC(0.75);
    
    c1->Modified();
    c1->Update();
    
    eOutScalar->SetLineColor(kYellow+10);
    eOutScalar->Draw("sames");
    
    c1->Update();
    
    TPaveStats *ps2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.6); ps2->SetX2NDC(0.75);
    //ps2->SetX1NDC(0.8); ps2->SetX2NDC(0.95);
    ps2->SetTextColor(kYellow+10);
    ps2->SetOptStat(1111);
    
    
    eOutVector->SetLineColor(kMagenta);
    eOutVector->Draw("sames");
    
    c1->Update();
    
    
    TPaveStats *ps3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.8); ps3->SetX2NDC(0.95);
    ps3->SetTextColor(kMagenta);
    
    TPaveText *t = new TPaveText(0.0, 0.95, 0.3, 1.0, "brNDC"); // left-up
    t->AddText(histName);
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    t->SetBorderSize(0);
    t->SetFillColor(gStyle->GetTitleFillColor());
    t->Draw();
    
    
    //TString chiSquareString= Form("%f",chiS);
    //drawLatex(.1, 0.92, chiSquareString);
    
    double res[eOutG4->GetSize()];
    double chi2=eOutG4->Chi2Test(eOutScalar,"UU P",res);
    drawLatex(.1,0.95,Form("%g",chi2));
    
    c1->Modified();
    c1->Update();
    
    //c1->SaveAs(histoFileNameEps);
    c1->SaveAs(histoFileNamePdf);
    
    
    chi2test_Plots(0, eOutG4, eOutScalar, process, variable, energy);
    
    
    //validatePdf(eOutScalar, eOutVector, eOutG4, energy);
    /*c1->Destructor();
     if(c1!=NULL)
     {
     c1->Clear();
     c1->Closed();
     c1 = NULL;
     }*/
}

//______________________________________________________________________________
void scaleHisto(TH1F* eOutScalar, TH1F* eOutVector, TH1F* eOutG4)
{
    
    //Scale the histograms and plot them
    TCanvas *c0 = new TCanvas("c0","Scaled histograms",200,10,700,500);
    
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    eOutScalarScaled->SetLineColor(kMagenta);
    eOutScalarScaled->Draw("][sames");
    c0->Update();
    
    
    TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    norm= eOutVectorScaled->GetEntries();
    eOutVectorScaled->Scale(1/norm);
    eOutVectorScaled->SetLineColor(kBlue);
    eOutVectorScaled->Draw("][sames");
    c0->Update();
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);
    eOutGeant4Scaled->SetLineColor(kGreen+2);
    eOutGeant4Scaled->Draw("][sames");
    c0->Update();
    
}



void table(Float_t x1, Float_t x2, Float_t yrange, TText *t, char **symbol, Bool_t octal);

void pstable(char **symbol, const char *process)
{
    
    TString tableName=Form("%s_pValuesTable.pdf", process);
    TString nameFileMerge=Form("%s-Validation.pdf", process);
    
    Float_t xrange = 80;
    Float_t yrange = 65;
    //Int_t w = 1150;
    Int_t w = 800;
    Int_t h = w*yrange/xrange;
    
    TCanvas *c3 = new TCanvas("c3","c3",240,20,w,h);
    //TCanvas *c3 = new TCanvas("c3","c3",240,20,600,h);
    c3->Range(0,0,xrange,yrange);
    
    TText *t = new TText(0,0,"a");
    t->SetTextSize(0.02); //Table content
    t->SetTextFont(62);
    t->SetTextAlign(22);
    
    table(0.5,xrange-0.5,yrange,t,symbol,0);
    
    //table(0.5,0.5*xrange-0.5,yrange,t,symbol,0);
    //table(0.5*xrange+0.3,xrange+0.5,yrange,t,symbol,0);
    TText *tlabel = new TText(0,0,"a");
    tlabel->SetTextSize(0.04);
    tlabel->SetTextColor(kAzure+5);
    tlabel->DrawText(0.5*xrange-10,yrange-4,"P-value Table");
    
    /*TText *tlabel1 = new TText(0,0,"a");
     tlabel1->SetTextSize(0.04);
     tlabel1->SetTextColor(kAzure+4);
     tlabel1->DrawText(0.5*xrange-10,yrange-8,"P-value Table");
     
     TText *tlabel2 = new TText(0,0,"a");
     tlabel2->SetTextSize(0.04);
     tlabel2->SetTextColor(kAzure+9);
     tlabel2->DrawText(0.5*xrange-10,yrange-10,"P-value Table");*/
    
    c3->Modified();
    c3->Update();
    c3->Print(tableName);
    c3->Print(nameFileMerge);
    
}


void table(Float_t x1, Float_t x2, Float_t yrange, TText *t, char **symbol, Bool_t octal)
{
    Int_t i;
    Int_t n = 0;
    for (i=0;i<1000;i++) {
        std::cout<<"symbol["<<i<<"]:"<<symbol[i]<<"\n";
        if (!strcmp(symbol[i],"END")) break;
        n++;
    }
    
    Float_t y1  = 1.5;
    Float_t y2  = yrange - 5.5;
    Float_t dx  = (x2-x1)/5;
    Float_t dy  = (y2 - 1 -y1)/(n/5+1);
    Float_t y   = y2 - 1 - 0.7*dy;
    Float_t xc0 = x1  + 0.5*dx;
    Float_t xc1 = xc0 + dx;
    Float_t xc2 = xc1 + dx;
    Float_t xc3 = xc2 + dx;
    Float_t xc4 = xc3 + dx;
    TLine *line = new TLine();
    line->DrawLine(x1,y1,x1,y2);//linea verticale sx
    line->DrawLine(x1,y1,x2,y1);//linea orizzontale bassa
    line->DrawLine(x1,y2,x2,y2);//linea orizzontale alta
    line->DrawLine(x2,y1,x2,y2);//linea verticale dx
    line->DrawLine(x1,y2-1.5,x2,y2-1.5); // linea orizzontale prima riga
    line->DrawLine(x1+  dx,y1,x1+  dx,y2);
    line->DrawLine(x1+2*dx,y1,x1+2*dx,y2);
    line->DrawLine(x1+3*dx,y1,x1+3*dx,y2);
    line->DrawLine(x1+4*dx,y1,x1+4*dx,y2);
    TText *tit = new TText(0,0,"a");
    tit->SetTextSize(0.02);//First raw table
    tit->SetTextFont(72);
    tit->SetTextAlign(22);
    tit->DrawText(xc0,y2-0.6,"EnergyIn");
    tit->DrawText(xc1,y2-0.6,"ValidationQ");
    tit->DrawText(xc2,y2-0.6,"Chi-2 test p-value");
    tit->DrawText(xc3,y2-0.6,"A-D test p-value");
    tit->DrawText(xc4,y2-0.6,"K-S test p-value");
    char text[15]; //BUG: il testo era troppo lungo!
    
    std::cout<<"N. entries: "<<n<<"\n";
    for (i=0;i<n;i+=5) {
        
        sprintf(text,"%s",symbol[i]);
        t->DrawText(xc0,y,text);
        sprintf(text,"%s",symbol[i+1]);
        t->DrawText(xc1,y,text);
        sprintf(text,"%s",symbol[i+2]);
        t->DrawText(xc2,y,text);
        sprintf(text,"%s",symbol[i+3]);
        t->DrawText(xc3,y,text);
        sprintf(text,"%s",symbol[i+4]);
        t->DrawText(xc4,y,text);
        y -= dy;
        
    }
}

void createPvTable(){
    
    std::cout<<"createPvTable, ciao.\n";
}



//______________________________________________________________________________
int main( int argc,  char *argv[])
{
    TApplication theApp("App",&argc,argv);
    
    /*
     TH1F *h1 = new TH1F("h1", "Gaussian plots", 100, -5, 5);
     TH1F *h2 = new TH1F("h2", "h2", 100, -5, 5);
     h1->FillRandom("gaus");
     h2->FillRandom("gaus");
     
     // Define the Canvas
     TCanvas *canvasProva = new TCanvas("canvasProva", "canvasProva", 800, 800);
     h1->Draw();
     canvasProva->SaveAs("primo.pdf");
     TCanvas *canvasProva2 = new TCanvas("canvasProva2", "canvasProva2", 800, 800);
     h2->Draw();
     canvasProva2->SaveAs("secondo.pdf");
     
     
     TCanvas *canvasMerge=new TCanvas("canvasMerge", "canvasMerge", 1600, 800);
     canvasMerge->Divide(2, 1);
     
     
     canvasMerge->cd(1);
     canvasProva->DrawClonePad();
     
     canvasMerge->cd(2);
     canvasProva2->DrawClonePad();
     canvasProva->Print("Merge.pdf(");*/
    
    TString nameFileMerge;
    char** symbolData;
    int i=10;
    
    const char* GUPhysicsModelName[5] = {
        "KleinNishina",  //Compton - gamma
        "BetheHeitler",  //Conversion - pair production - gamma
        "SauterGavrila", // Photo-Electric Effect - gamma
        "MollerBhabha",  // Ionization -electron
        "SeltzerBerger"  // Bremsstrahlung
    };
    
    const char* GUPhysicsModelVariable[5] = {
        "EnergyIn",
        "EnergyOut1",
        "EnergyOut2",
        "AngleOut1",
        "AngleOut2"
    };
    
    const char *pvalues[] =
    {"10 KeV","En_out1","0.650629","1","1","","En_out2","0.650629","1","1",
        "","Angle_out1","0.650629","1","1", "","Angle_out1","0.650629","1","1",
        "","En_out1","0.650629","1","1", "10 KeV","En_out1","0.650629","1","1",
        "10 KeV","En_out1","0.650629","1","1","END"};
    
    
    
    //gaxis();
    //TString modelName;
    //TString variableName;
    TString  in1, in2, in3[10], *in4;
    const char *phisicsModel, *variable, *energy[10], *nEnergies;
    int nEn;
    
    TString myInput;
    std::cout<<"Starting.\n";
    
    
    TCanvas *firstPageValidation=new TCanvas("Validation","Validation",800,200,1200,1000);
    drawText(.6,.45,Form("G4 vs GeantV validation"));
    //firstPageValidation->Print("Merge.pdf(");
    
    
    do{
        
        std::cout<<"Press 'D' for the default execution, press 'C' for the customized execution, 'T' to run statistical Tests, 'CT' to create the p-Value Table.\n";
        
        std::cin>>myInput;
        //All models, all validation quantities, one energy
        if(myInput.EqualTo("D")|| myInput.EqualTo("d"))
        {
            std::cout<<"Insert the energy of the projectile [Mev]\n";
            std::cin>>in3[i];
            energy[i]=in3[i].Data();
            std::cout << "Starting default mode.\n";
            
            //for(int i=0; i<5 ; i++)
            //genAllHisto(GUPhysicsModelName[i], energy);
            
        }
        //One model, all validation quantities, different energies
        else if(myInput.EqualTo("C")|| myInput.EqualTo("c"))
        {
            std::cout<<"PhysicsModel: KleinNishina, BetheHeitler, SauterGavrila, MollerBhabha, SeltzerBerger.\n";
            std::cin>>in1;
            phisicsModel=in1.Data();
            //std::cout<<"Variable: EnergyIn, EnergyOut1, EnergyOut2, AngleOut1, AngleOut2, or all.\n";
            //std::cin>>in2;
            //variable=in2.Data();
            variable="all";
            
            nameFileMerge=Form("%s-Validation.pdf",phisicsModel);
            TString cazzu=Form("%s-Validation.pdf(",phisicsModel);
            firstPageValidation->Print(cazzu);
            
            std::cout<<"How many energies do you want to test?\n";
            std::cin>>in2;
            nEnergies=in2.Data();
            nEn=atoi(nEnergies);
            
            TString* energy_string;
            TString* chi2_string;
            TString* AD_string;
            TString* KS_string;
            
            
            energy_string=new TString[nEn];
            chi2_string=new TString[nEn*4];
            AD_string=new TString[nEn*4];
            KS_string=new TString[nEn*4];
            
            
            double chi2_pval[4], AD_pval[4], KS_pval[4];
            symbolData=new char*[nEn*20+1];
            
            for (int i=0; i<nEn; i++) //energie
            {
                std::cout<<"Energy of the projectile [Mev]\n";
                std::cin>>in3[i];
                energy[i]=in3[i].Data();
                
                energy_string[i]=Form("%s MeV",energy[i]);
                
                
                genAllHisto(phisicsModel,energy[i], chi2_pval, AD_pval, KS_pval);
                
                
                for (int k=0; k<5; k++) //colonne
                    for (int j=0; j<4; j++) //righe
                    {
                        if(k==0)
                            if(j==0)
                            {
                                
                                symbolData[i*20]=(char*)energy_string[i].Data();
                            }
                            else
                            {
                                //sprintf( symbolData[i*20+j*5], "");
                                symbolData[i*20+j*5]=(char*)"";
                            }
                            else
                                if(k==1)
                                {
                                    symbolData[i*20+j*5+k]=(char*)GUPhysicsModelVariable[j+1];
                                    
                                }
                                else
                                    if(k==2)
                                    {
                                        chi2_string[4*i+j]= Form("%g",chi2_pval[j]);
                                        symbolData[i*20+j*5+k]=(char *)chi2_string[4*i+j].Data();
                                        std::cout<<"chi2_pval symbolData["<<i*20+j*5+k<<"]: "<<symbolData[i*20+j*5+k]<<"\n";
                                    }
                                    else
                                        if(k==3)
                                        {
                                            AD_string[4*i+j]= Form("%g",AD_pval[j]);
                                            symbolData[i*20+j*5+k]=(char *)AD_string[4*i+j].Data();
                                            std::cout<<"AD_pval symbolData["<<i*20+j*5+k<<"]: "<<symbolData[i*20+j*5+k]<<"\n";
                                        }
                                        else
                                            if(k==4)
                                            {
                                                KS_string[4*i+j]= Form("%g",KS_pval[j]);
                                                symbolData[i*20+j*5+k]=(char *)KS_string[4*i+j].Data();
                                                std::cout<<"KS_pval symbolData["<<i*20+j*5+k<<"]: "<<symbolData[i*20+j*5+k]<<"\n";
                                            }
                        
                        
                    }
                
            }
            symbolData[nEn*20]=(char*)"END";
            
            /*std::cout<<"Stampiamo sta caxxo de tabella\n";
             std::cout<<"Energie "<<energy[0]<<" e " << energy[1]<<"\n";
             for (int stamico=0; stamico<nEn*20; stamico++)
             std::cout<<"symbolData["<<stamico<<"]: "<<symbolData[stamico]<<"\n";*/
            
            pstable(symbolData, phisicsModel);
            
        }
        else if(myInput.EqualTo("T")|| myInput.EqualTo("t"))
        {
            
            goftest();
            std::cout<<"PhysicsModel: KleinNishina, BetheHeitler, SauterGavrila, MollerBhabha, SeltzerBerger.\n";
            std::cin>>in1;
            phisicsModel=in1.Data();
            std::cout<<"Variable: EnergyIn, EnergyOut1, EnergyOut2, AngleOut1, AngleOut2.\n";
            std::cin>>in2;
            variable=in2.Data();
            
            std::cout<<"Energy of the projectile [Mev]\n";
            std::cin>>in3[i];
            energy[i]=in3[i].Data();
            //if(in2.EqualTo("all")||in2.EqualTo("All"))
            //  goftestAllMod(phisicsModel,energy);
            //else
            goftestSingleMod(phisicsModel, variable, energy[i]);
            
            
        }
        
        else if(myInput.EqualTo("CT")|| myInput.EqualTo("ct")|| myInput.EqualTo("Ct")|| myInput.EqualTo("cT"))
        {
            createPvTable();
        }
        
        std::cout<<"Press 'c' to continue, 'e' to exit\n";
        std::cin>>myInput;
        
    } while(myInput.EqualTo("c")|| myInput.EqualTo("C"));
    
    
    /*
     if (argc<2)
     {
     std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n \n<ENERGY> [MeV] \n";
     std::cout << "Default mode usage: "<< argv[0] <<" <ENERGY> [MeV] \n";
     return 0;
     }
     
     int argIndex = 1;
     bool default=true;
     
     while (argIndex < argc){
     
     if (!strcmp(argv[argIndex],"-default")) {
     ++argIndex;
     
     ++argIndex;
     }
     
     else if (!strcmp(argv[argIndex],"-model")) {
     ++argIndex;
     modelName= theApp.Argv(argIndex);
     ++argIndex;
     }
     
     else if (!strcmp(argv[argIndex],"-histo")) {
     ++argIndex;
     variableName= theApp.Argv(argIndex);
     ++argIndex;
     }
     else if (!strcmp(argv[argIndex],"-energy")) {
     ++argIndex;
     energy= theApp.Argv(argIndex);
     ++argIndex;
     }
     
     }
     
     
     
     if(argc==2)
     {
     std::cout << "Starting default mode.\n";
     for(int i=0; i<5 ; i++)
     for(int j=0; j<5; j++)
     {
     genSingleHisto(GUPhysicsModelName[i], GUPhysicsModelVariable[j], theApp.Argv(1), i+j);
     std::cout << "Faccio altro, magari .\n";
     }
     std::cout << "Run completed.\n";
     theApp.Run();
     return 0;
     }
     
     
     else if (argc<4)
     {
     std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n \n<ENERGY> [MeV] \n";
     return 0;
     }
     
     
     else
     if( !(strcmp(argv[1], "KleinNishina")==0) && !(strcmp(argv[1], "BetheHeitler")==0) && !(strcmp(argv[1], "SauterGravila")==0) && !(strcmp(argv[1], "MollerBhabha")==0) && !(strcmp(argv[1], "SeltzerBerger")==0))
     {
     std::cout<<"Wrong PHYSICSMODEL_NAME.\n";
     std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n\n <ENERGY> [MeV] \n";
     return 0;
     }
     
     else if( !(strcmp(argv[2], "EnergyIn")==0) && !(strcmp(argv[2], "EnergyOut1")==0) && !(strcmp(argv[2], "EnergyOut2")==0) && !(strcmp(argv[2], "AngleOut1")==0) && !(strcmp(argv[2], "AngleOut2")==0))
     {
     std::cout<<"Wrong HISTOGRAM_NAME.\n";
     std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n \n <ENERGY> [MeV] \n";
     return 0;
     }
     
     std::cout << "Starting.\n";
     genSingleHisto(theApp.Argv(1), theApp.Argv(2), theApp.Argv(3),1);
     
     TCanvas * newCanvas=new TCanvas("c1","c",200,10,700,500);
     
     double y = 0.95;
     for (int f = 12; f<=152; f+=10) {
     if (f!=142) drawtext(0.02,y, f,"ABCDEFGH abcdefgh 0123456789 @#$");
     else drawtext(0.02,y, f,"ABCD efgh 01234 @#$");
     y -= 0.065;
     }
     genSingleHisto("KleinNishina", "AngleOut2", "500",6);*/
    
    
    /*
     TString processName= theApp.Argv(1);
     TString histogramName= theApp.Argv(2);
     
     TString histPath= Form("%s/%s",theApp.Argv(1), theApp.Argv(2));
     TString histName= Form("%s/%s/%sMeV",theApp.Argv(1), theApp.Argv(2), theApp.Argv(3));
     TString histoFileName= Form("%s-%s-%sMeV.pdf",theApp.Argv(1), theApp.Argv(2), theApp.Argv(3));
     
     TString g4HistoName=Form("%s/geant4",theApp.Argv(1));
     TString gvHistoScalarName=Form("%s/geantVscalar",theApp.Argv(1));
     TString gvHistoVectorName=Form("%s/geantVvector",theApp.Argv(1));
     
     TString geant4RootFileName= Form("geant4_%sMeV.root",theApp.Argv(3));
     TString scalarRootFileName= Form("scalar_%sMeV.root",theApp.Argv(3));
     TString vectorRootFileName= Form("vector_%sMeV.root",theApp.Argv(3));
     TString pdfFileName=Form("pdf_%sMeV.root",theApp.Argv(3));
     
     //Geant 4
     
     TFile *fG4 = new TFile(geant4RootFileName,"w");
     fG4->ls();
     std::cout<<"Path to the histogram:"<< histPath<<std::endl;
     //TH1F *eOutG4 = (TH1F*)fG4->Get("KleinNishina /EnergyOut1");
     //eOutG4->SetName("KleinNishina/geant4");
     
     TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
     if(eOutG4) eOutG4->SetName("geant4");
     else std::cout<<"eOutG4 is null\n";
     
     
     //GeantV scalar
     TFile *fScalar = new TFile(scalarRootFileName,"w");
     //TH1F *eOutScalar = (TH1F*)fScalar->Get("KleinNishina /EnergyOut1");
     //eOutScalar->SetName("KleinNishina/scalar");
     
     TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
     
     if(eOutScalar) eOutScalar->SetName("geantVscalar");
     else std::cout<<"eOutScalar is null\n";
     
     
     //GeantV vector
     TFile *fVector = new TFile(vectorRootFileName,"w");
     //TH1F *eOutScalar = (TH1F*)fScalar->Get("KleinNishina /EnergyOut1");
     //eOutScalar->SetName("KleinNishina/scalar");
     
     TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
     if(eOutVector) eOutVector->SetName("geantVvector");
     else std::cout<<"eOutVector is null\n";
     
     */
    
    
    
    
    
    /*
     //Plot the not scaled histograms
     TCanvas *c1 = new TCanvas("c1","Physics validation",200,10,700,500);
     eOutG4->Draw();
     c1->Update();
     
     TPaveStats *ps1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
     ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.55);
     //ps1->SetX1NDC(0.5); ps1->SetX2NDC(0.7);
     
     c1->Modified();
     c1->Update();
     
     eOutScalar->SetLineColor(kYellow+10);
     eOutScalar->Draw("][sames");
     c1->Update();
     
     TPaveStats *ps2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
     ps2->SetX1NDC(0.6); ps2->SetX2NDC(0.75);
     //ps2->SetX1NDC(0.75); ps2->SetX2NDC(0.95);
     ps2->SetTextColor(kYellow+10);
     
     eOutVector->SetLineColor(kMagenta);
     eOutVector->Draw("][sames");
     c1->Update();
     
     TPaveStats *ps3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
     ps3->SetX1NDC(0.8); ps3->SetX2NDC(0.95);
     ps3->SetTextColor(kMagenta);
     
     TPaveText *t = new TPaveText(0.0, 0.95, 0.3, 1.0, "brNDC"); // left-up
     t->AddText(histName);
     
     gROOT->SetStyle("Plain");
     //gStyle->SetOptStat(0000);
     gStyle->SetOptFit(1111);
     gStyle->SetOptTitle(0);
     t->SetBorderSize(0);
     t->SetFillColor(gStyle->GetTitleFillColor());
     t->Draw();
     
     c1->Modified();
     c1->Update();
     */
    
    
    
    /*
     
     //PDF validation
     Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
     //eOutScalar->GetEntries();
     //std::cout<<"N of histogram entries: "<<eOutScalar->GetEntries()<<"\n";
     
     
     Double_t xScalar[entries], yScalar[entries], zScalar[entries];
     Double_t xVector[entries], yVector[entries], zVector[entries];
     Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
     
     TFile *fPfd = new TFile(pdfFileName,"w");
     TGraph *pdfGraph= (TGraph*)fPfd->Get("Pdf2.0"); //Get the first pdf, for fMinX value
     
     
     TString c2Name=Form("Pdf %s MeV",theApp.Argv(3));
     
     TCanvas *c2 = new TCanvas("c2",c2Name,200,10,700,500);
     pdfGraph->SetLineColor(kBlue);
     //pdfGraph->SetTitle("pdf; mcs; dm");
     
     pdfGraph->Draw();
     c2->Update();
     
     
     
     for(int j = 0; j < entries ; ++j){
     yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
     pdfGraph->GetPoint(j, xScalar[j],zScalar[j]); //to map
     
     yVector[j] = eOutVectorScaled->GetBinContent(j+1);
     pdfGraph->GetPoint(j, xVector[j],zVector[j]); //to map
     
     yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
     pdfGraph->GetPoint(j, xGeant4[j],zGeant4[j]); //to map
     }
     
     
     //Read the pdf File for 500MeV
     TCanvas *c3 = new TCanvas("c3","Graph comparison",200,10,700,500);
     pdfGraph->SetLineColor(kBlue);
     pdfGraph->Draw();
     c3->Update();
     
     
     
     //Create the Graph with x-values taken from the pdfGraph and the y-values taken from the simulation histogram
     TGraph *myCopyVectorGr=new TGraph(entries,xVector,yVector);
     //myCopyVectorGr->SetLineColor(kYellow+10);
     //myCopyVectorGr->SetLineStyle(2);
     //myCopyVectorGr->SetLineWidth(2);
     //myCopyVectorGr->Draw("LP");
     
     //TMarker *mark=new TMarker();
     //mark->SetMarkerStyle(2);
     //mark->SetMarkerSize(1);
     
     TGraph *myCopyScalarGr=new TGraph(entries,xScalar,yScalar);
     myCopyScalarGr->SetLineColor(kYellow+10);
     //myCopyScalarGr->SetLineStyle(4);
     //myCopyScalarGr->SetLineWidth(2);
     
     //myCopyScalarGr->SetMarkerStyle(20);
     myCopyScalarGr->Draw("LP");
     
     //TMarker *mark2=new TMarker();
     //mark2->SetMarkerStyle(4);
     
     
     TGraph *myCopyG4Gr=new TGraph(entries,xGeant4,yGeant4);
     myCopyG4Gr->SetLineColor(kRed);
     //myCopyG4Gr->SetLineStyle(5);
     //myCopyG4Gr->SetLineWidth(2);
     myCopyG4Gr->Draw("LP");
     
     TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
     leg->SetHeader("Legend");
     //leg->AddEntry(my,"Histogram filled with random numbers","f");
     //leg->AddEntry(my,"Function abs(#frac{sin(x)}{x})","l");
     //leg->AddEntry("gr","Graph with error bars","lep");
     TString legName=Form("Pdf for E_in=%s MeV",theApp.Argv(3));
     leg->AddEntry(pdfGraph,legName,"l");
     leg->AddEntry(myCopyScalarGr,"E_out Histogram values scaled - Scalar","l");
     //leg->AddEntry(myCopyVectorGr,"E_out Histogram values scaled - Vector","l");
     leg->AddEntry(myCopyG4Gr,"E_out Histogram values scaled - Geant4","l");
     
     leg->Draw();
     c3->Update();
     std::cout<<"5\n";
     
     ////Calculate chi-square
     std::cout<<"#E_out Histogram entries: "<<myCopyScalarGr->GetN()<<"#Pdf entries: "<<pdfGraph->GetN()<<"\n";
     double chiSquare_vector=0,chiSquare_scalar=0, chiSquare_geant4=0, expX,expY ;
     for (int i=0; i<pdfGraph->GetN(); i++)
     {
     pdfGraph->GetPoint(i, expX, expY);
     //myCopyScalarGr->GetPoint(i, obsX, obsY);
     //std::cout<<"expX: "<<expX<<" expY: "<<expY<<" obsX: "<<obsX<<" obsY: "<<obsY<<std::endl;
     chiSquare_scalar+=(yScalar[i]-expY)*(yScalar[i]-expY)/expY;
     chiSquare_vector+=(yVector[i]-expY)*(yVector[i]-expY)/expY;
     chiSquare_geant4+=(yGeant4[i]-expY)*(yGeant4[i]-expY)/expY;
     }
     std::cout<<"chiSquare scalar: "<<chiSquare_scalar<<" chiSquare-reduced: "<<chiSquare_scalar/(eOutScalarScaled->GetSize()-2)<<std::endl;
     std::cout<<"chiSquare vector: "<<chiSquare_vector<<" chiSquare-reduced: "<<chiSquare_vector/(eOutVectorScaled->GetSize()-2)<<std::endl;
     std::cout<<"chiSquare geant4: "<<chiSquare_geant4<<" chiSquare-reduced: "<<chiSquare_geant4/(eOutGeant4Scaled->GetSize()-2)<<std::endl;
     
     
     //gr[i]->SetTitle(graphName);
     //gr[i]->GetXaxis()->SetTitle(" EnergyOut [MeV]");
     //gr[i]->GetYaxis()->SetTitle("Cross section");
     //myCopyGr->Write();
     c1->SaveAs(histoFileName);
     */
    
    //canvasProva2->Print("Merge.pdf");
    //canvasMerge->Print("Merge.pdf)");
    std::cout << "\nRun completed.\n";
    
    //theApp.Run();
    TCanvas *lastPageValidation=new TCanvas("Validation","Validation",800,200,1200,1000);
    //drawText(.6,0.45,"Thanks");
    //lastPageValidation->Print("Merge.pdf)");
    
    TString cazzucazzu=Form("%s-Validation.pdf)",phisicsModel);
    lastPageValidation->Print(cazzucazzu);
    return 0;
}
