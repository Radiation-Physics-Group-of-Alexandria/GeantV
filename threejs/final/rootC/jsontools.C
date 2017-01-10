#if defined(G__WIN32) && defined(__CINT__) && !defined(__MAKECINT__)
{
   Info("jsontools.C", "Has to be run in compiled mode ... doing this for you.");
   gSystem->CompileMacro("jsontools.C");
   jsontools();
}
#else

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoBoolNode.h"
#include "TGeoMatrix.h"
#include "TGeoCompositeShape.h"
#include "TBuffer3D.h"
#include "TString.h"
#include "TColor.h"
#include "Riostream.h"
#include <algorithm>
#include <queue>
#include <math.h>
#ifndef VAR_STORE
#define  VAR_STORE
#include "varStore.h"
#endif

#include <geomtools.C>

using namespace std;


int getMatrixVR(varStore *prm,TGeoMatrix *matrix){
  int nv=0;
  Double_t v[16];

  if (!matrix->IsIdentity() ) {

    Double_t v[16];

    matrix->GetHomogenousMatrix(v);

    Double_t t;
    t=v[3];
    v[3]=v[12];
    v[12]=t;

    t=v[7];
    v[7]=v[13];
    v[13]=t;

    t=v[11];
    v[11]=v[14];
    v[14]=t;

    v[15]=1;
    nv=16;
    for (int i=0;i<nv;i++)prm->Mtrx[i]=v[i];
  } 
  prm->nM=nv;
  return nv;
}



int exportVolGeometry(varStore *prm,TGeoShape *shape,Int_t geomId){
  
  TString clsname  = shape->ClassName();

  if (clsname == "TGeoShapeAssembly" ) {

    ExportGeomAssembly(prm);
  } else if (clsname == "TGeoBBox") {
    
    Double_t  BDX=((TGeoBBox *)shape)->GetDX();
    Double_t  BDY=((TGeoBBox *)shape)->GetDY();
    Double_t  BDZ=((TGeoBBox *)shape)->GetDZ();

    ExportGeomBBox(prm,geomId, BDX,BDY,BDZ );

  } else if (clsname == "TGeoParaboloid") {

    Double_t  Rlo=((TGeoParaboloid *)shape)->GetRlo(); 
    Double_t  Rhi=((TGeoParaboloid *)shape)->GetRhi();
    Double_t  Dz=((TGeoParaboloid *)shape)->GetDz();

    ExportGeomParabloid(prm,geomId, Rlo,Rhi,Dz);

  } else if (clsname == "TGeoSphere") {
    Double_t  Rmin= ((TGeoSphere *)shape)->GetRmin();
    Double_t  Rmax= ((TGeoSphere *)shape)->GetRmax();
    Double_t  thetaStart= ((TGeoSphere *)shape)->GetTheta1()*TMath::DegToRad();
    Double_t  thetaLength= ((TGeoSphere *)shape)->GetTheta2()*TMath::DegToRad() - thetaStart; 
    Double_t  phiStart= ((TGeoSphere *)shape)->GetPhi1()*TMath::DegToRad();
    Double_t  phiLength= ((TGeoSphere *)shape)->GetPhi2()*TMath::DegToRad() - phiStart;

    ExportGeomSphere(prm,geomId, Rmin, Rmax, thetaStart,thetaLength,phiStart,phiLength);

  } else if (clsname == "TGeoConeSeg") {
    Double_t Dz=((TGeoConeSeg *)shape)->GetDz();
    Double_t Rmin1=((TGeoConeSeg *)shape)->GetRmin1();
    Double_t Rmin2=((TGeoConeSeg *)shape)->GetRmin2();
    Double_t Rmax1=((TGeoConeSeg *)shape)->GetRmax1();
    Double_t Rmax2=((TGeoConeSeg *)shape)->GetRmax2();
    Double_t phiStart=((TGeoConeSeg *)shape)->GetPhi1()*TMath::DegToRad();
    Double_t phiLength=((TGeoConeSeg *)shape)->GetPhi2()*TMath::DegToRad()-phiStart;

    ExportGeomConeSeg(prm,geomId, Dz, Rmin1,Rmin2,Rmax1,Rmax2,phiStart,phiLength );

  } else if (clsname == "TGeoCone") {

    Double_t Dz=((TGeoCone *)shape)->GetDz();
    Double_t Rmin1=((TGeoCone *)shape)->GetRmin1();
    Double_t Rmin2=((TGeoCone *)shape)->GetRmin2();
    Double_t Rmax1=((TGeoCone *)shape)->GetRmax1();
    Double_t Rmax2=((TGeoCone *)shape)->GetRmax2();

    ExportGeomConeSeg(prm,geomId, Dz, Rmin1,Rmin2,Rmax1,Rmax2);

  } else if (clsname == "TGeoTubeSeg") {

    Double_t Dz=((TGeoTubeSeg *)shape)->GetDz();
    Double_t Rmin=((TGeoTubeSeg *)shape)->GetRmin();
    Double_t Rmax=((TGeoTubeSeg *)shape)->GetRmax();
    Double_t phiStart=((TGeoTubeSeg *)shape)->GetPhi1()*TMath::DegToRad();
    Double_t phiLength=((TGeoTubeSeg *)shape)->GetPhi2()*TMath::DegToRad()- phiStart;

    ExportGeomConeSeg(prm ,geomId, Dz, Rmin, Rmin, Rmax, Rmax, phiStart, phiLength);
  } else if (clsname == "TGeoTube") {

    Double_t Dz=((TGeoTube *)shape)->GetDz();
    Double_t Rmin=((TGeoTube *)shape)->GetRmin();
    Double_t Rmax=((TGeoTube *)shape)->GetRmax();

    ExportGeomTub(prm,geomId, Rmin,Rmax,Dz);

  } else if (clsname == "TGeoPcon") {
    Double_t PerimR[100],PerimZ[100];
    Double_t phiStart=((TGeoPcon *)shape)->GetPhi1()*TMath::DegToRad();
    Double_t phiLength=((TGeoPcon *)shape)->GetDphi()*TMath::DegToRad();
    
    int nZ       =((TGeoPcon *)shape)->GetNz();
    for (int i=0;i<nZ;i++){
      Double_t pZ = ((TGeoPcon *)shape)->GetZ(i);
      Double_t pRmin = ((TGeoPcon *)shape)->GetRmin(i);
      Double_t pRmax = ((TGeoPcon *)shape)->GetRmax(i);

      PerimR[i]=pRmax ;
      PerimZ[i]=pZ;
      PerimR[2*nZ-1-i]=pRmin ;
      PerimZ[2*nZ-1-i]=pZ;
    }
    PerimR[2*nZ]=PerimR[0];
    PerimZ[2*nZ]=PerimZ[0];

   ExportGeomConeSeg1(prm,geomId, 2*nZ,  PerimR, PerimZ, phiStart, phiLength );

  } else if (clsname == "TGeoTorus") {

    Double_t R   =((TGeoTorus *)shape)->GetR();
    Double_t Rmin=((TGeoTorus *)shape)->GetRmin();
    Double_t Rmax=((TGeoTorus *)shape)->GetRmax();
    Double_t phiStart=((TGeoTorus *)shape)->GetPhi1()*TMath::DegToRad();
    Double_t phiLength=((TGeoTorus *)shape)->GetDphi()*TMath::DegToRad();

    ExportGeomTorus(prm,geomId, R,  Rmin, Rmax, phiStart, phiLength );
 
  } else if (clsname == "TGeoPgon") {


    Double_t PerimR[100],PerimZ[100];

    int nedg=((TGeoPgon *)shape)->GetNedges();
    Double_t phiStart=((TGeoPgon *)shape)->GetPhi1()*TMath::DegToRad();
    Double_t phiLength=((TGeoPgon *)shape)->GetDphi()*TMath::DegToRad();

    int nZ       =((TGeoPgon *)shape)->GetNz();
    for (int i=0;i<nZ;i++){

      Double_t pZ = ((TGeoPgon *)shape)->GetZ(i);
      Double_t pRmin = ((TGeoPgon *)shape)->GetRmin(i);
      Double_t pRmax = ((TGeoPgon *)shape)->GetRmax(i);
      
      PerimR[i]=pRmax ;
      PerimZ[i]=pZ;
      PerimR[2*nZ-1-i]=pRmin ;
      PerimZ[2*nZ-1-i]=pZ;
    }
    PerimR[2*nZ]=PerimR[0];
    PerimZ[2*nZ]=PerimZ[0];
    for (int i=0;i<nZ;i++){

      for (int j=0;j<=nZ;j++){
        if ((PerimR[j]==PerimR[2*nZ-1-i])&& (PerimZ[j]==PerimZ[2*nZ-1-i])){
          PerimR[j]+=0.0001;
        }
      }
    }
 
    ExportGeomConeSeg2(prm,geomId, nedg, 2*nZ,  PerimR, PerimZ, phiStart, phiLength );
    

  } else if (clsname == "TGeoHype") {


    Double_t  Rmin=((TGeoHype *)shape)->GetRmin(); 
    Double_t  Rmax=((TGeoHype *)shape)->GetRmax();
    Double_t  StIn=((TGeoHype *)shape)->GetStIn(); 
    Double_t  StOut=((TGeoHype *)shape)->GetStOut();
    Double_t  Dz=((TGeoHype *)shape)->GetDz();

    ExportGeomHype(prm,geomId,Rmin,Rmax,StIn,StOut,Dz);

  } else if (clsname == "TGeoScaledShape") {
  } else if (clsname == "TGeoArb8") {

    Double_t x[8],y[8];
    Double_t dz=((TGeoArb8 *)shape)->GetDz();   
    Double_t z[]={-dz,-dz,-dz,-dz,dz,dz,dz,dz};

    for (int i=0;i<8;i++){
      x[i]=((TGeoArb8 *)shape)->GetVertices()[2*i];
      y[i]=((TGeoArb8 *)shape)->GetVertices()[2*i+1];
    }

    ExportGeomQuad(prm,geomId, 8, x, y, z, 1);

  } else if (clsname == "TGeoPara") {

    /*     Parallelepiped class.

        It has 6 parameters :

        dx, dy, dz - half lengths in X, Y, Z
        alpha - angle w.r.t the Y axis from center of low Y edge to center of high Y edge [deg]
        theta, phi - polar and azimuthal angles of the segment between low and high Z surfaces [deg]

    */

    Double_t dx=((TGeoPara *)shape)->GetX();
    Double_t dy=((TGeoPara *)shape)->GetY();
    Double_t dz=((TGeoPara *)shape)->GetZ();
    Double_t alpha=((TGeoPara *)shape)->GetAlpha()*TMath::DegToRad();
    Double_t theta=((TGeoPara *)shape)->GetTheta()*TMath::DegToRad();
    Double_t phi=((TGeoPara *)shape)->GetPhi()*TMath::DegToRad();

    ExportGeomPara(prm,geomId,dx,dy,dz,alpha,theta,phi);


  } else if (clsname == "TGeoTrap") {
    /*
      TRAP is a general trapezoid, i.e.

      one for which the faces perpendicular to z are trapezia and their centres are not the same x, y. 
      It has 11 parameters: 
      the half length in z, 
      the polar angles from the centre of the face at low z to that at high z, 
      H1 the half length in y at low z, 
      LB1 the half length in x at low z and y low edge, 
      LB2 the half length in x at low z and y high edge, 
      TH1 the angle w.r.t. the y axis from the centre of low y edge to the centre of the high y edge, and 
      H2, LB2, LH2, TH2, the corresponding quantities at high z.

    */
    Double_t dz=    ((TGeoTrap *)shape)->GetDz();
    Double_t theta= ((TGeoTrap *)shape)->GetTheta()*TMath::DegToRad();
    Double_t phi=   ((TGeoTrap *)shape)->GetPhi()*TMath::DegToRad();
    Double_t bl1=   ((TGeoTrap *)shape)->GetBl1();
    Double_t tl1=   ((TGeoTrap *)shape)->GetTl1();
    Double_t bl2=   ((TGeoTrap *)shape)->GetBl2();
    Double_t tl2=   ((TGeoTrap *)shape)->GetTl2();
    Double_t h1=    ((TGeoTrap *)shape)->GetH1();
    Double_t h2=    ((TGeoTrap *)shape)->GetH2();

    Double_t alpha1=((TGeoTrap *)shape)->GetAlpha1()*TMath::DegToRad();
    Double_t alpha2=((TGeoTrap *)shape)->GetAlpha2()*TMath::DegToRad();

    ExportGeomTrap(prm,geomId,    dz,theta,phi,bl1,tl1,bl2,tl2,h1,h2,alpha1,alpha2 );


  } else if (clsname == "TGeoGtra") {


    /* Gtra is a twisted trapezoid.

       i.e. one for which the faces perpendicular to z are trapezia and their centres are not the same x, y. It has 12 parameters: 
       the half length in z, 
       the polar angles from the centre of the face at low z to that at high z, 
       twist, 
       H1 the half length in y at low z, 
       LB1 the half length in x at low z and y low edge, 
       LB2 the half length in x at low z and y high edge, 
       TH1 the angle w.r.t. the y axis from the centre of low y edge to the centre of the high y edge, and 
       H2, LB2, LH2, TH2, the corresponding quantities at high z.

    */

    Double_t dz=    ((TGeoGtra *)shape)->GetDz();
    Double_t theta= ((TGeoGtra *)shape)->GetTheta()*TMath::DegToRad();
    Double_t phi=   ((TGeoGtra *)shape)->GetPhi()*TMath::DegToRad();
    Double_t bl1=   ((TGeoGtra *)shape)->GetBl1();
    Double_t tl1=   ((TGeoGtra *)shape)->GetTl1();
    Double_t bl2=   ((TGeoGtra *)shape)->GetBl2();
    Double_t tl2=   ((TGeoGtra *)shape)->GetTl2();
    Double_t h1=    ((TGeoGtra *)shape)->GetH1();
    Double_t h2=    ((TGeoGtra *)shape)->GetH2();
    Double_t alpha1=((TGeoGtra *)shape)->GetAlpha1()*TMath::DegToRad();
    Double_t alpha2=((TGeoGtra *)shape)->GetAlpha2()*TMath::DegToRad();

    Double_t twist= ((TGeoGtra *)shape)->GetTwistAngle()*TMath::DegToRad();

    ExportGeomGTra(prm,geomId,    dz,theta,phi,bl1,tl1,bl2,tl2,h1,h2,alpha1,alpha2,twist );

  } else if (clsname == "TGeoTrd1") {

    Double_t dx1=  ((TGeoTrd1 *)shape)->GetDx1();
    Double_t dx2=  ((TGeoTrd1 *)shape)->GetDx2();
    Double_t dy1=  ((TGeoTrd1 *)shape)->GetDy();
    Double_t dy2=  dy1;
    Double_t dz=   ((TGeoTrd1 *)shape)->GetDz();

    ExportGeomTrd(prm,geomId,  dx1,dx2,dy1,dy2,dz);

  } else if (clsname == "TGeoTrd2") {

    Double_t dx1=  ((TGeoTrd2 *)shape)->GetDx1();
    Double_t dx2=  ((TGeoTrd2 *)shape)->GetDx2();
    Double_t dy1=  ((TGeoTrd2 *)shape)->GetDy1();
    Double_t dy2=  ((TGeoTrd2 *)shape)->GetDy2();
    Double_t dz=   ((TGeoTrd2 *)shape)->GetDz();

    ExportGeomTrd(prm,geomId,  dx1,dx2,dy1,dy2,dz);

  } else if (clsname == "TGeoCtub") {

    Double_t  Rmin=((TGeoCtub *)shape)->GetRmin(); 
    Double_t  Rmax=((TGeoCtub *)shape)->GetRmax();
    Double_t  Dz=((TGeoCtub *)shape)->GetDz();

    ExportGeomTub(prm,geomId,Rmin,Rmax,Dz);

  } else if (clsname == "TGeoEltu") {
    Double_t A=((TGeoEltu *)shape)->GetA(); //dx
    Double_t B=((TGeoEltu *)shape)->GetB(); //dy
    Double_t Dz=((TGeoEltu *)shape)->GetDz(); //dz

    ExportGeomEltub(prm,geomId, A,B,Dz);
         
  } else if (clsname == "TGeoXtru") {

    Double_t x[1000],y[1000];
    Double_t z[100],px[100],py[100],f[100];

    int np=((TGeoXtru *)shape)->GetNvert();
    int nz=((TGeoXtru *)shape)->GetNz();
    
    for (int i = 0; i < np; i++) {
      x[i]=((TGeoXtru *)shape)->GetX(i);
      y[i]=((TGeoXtru *)shape)->GetY(i);
    }

    for (int i = 0; i < nz; i++) {
      z[i] =((TGeoXtru *)shape)->GetZ(i);
      px[i]=((TGeoXtru *)shape)->GetXOffset(i);
      py[i]=((TGeoXtru *)shape)->GetYOffset(i);
      f[i] =((TGeoXtru *)shape)->GetScale(i);
    }

    ExportGeomExtrude(prm,geomId,  np,nz,x,y,z,px,py,f); 
  }
  return prm->kVolume;
}


int defJSONcolors(ofstream &PrimFile){
  int ncol=TColor::GetNumberOfColors();

  PrimFile << "\n\"Colors\":[";

  for (int i=0;i<ncol;i++){

    TColor *color = gROOT->GetColor(i);
    PrimFile << "\t{\"cn\":\"";
    PrimFile <<color->GetName()<<"\",\"cl\":";

    TString ColorHex=color->AsHexString();
    ColorHex.Remove(0,1);
    PrimFile <<"\"#"<<ColorHex<<"\"}";
    if (i<ncol-1) PrimFile <<",\n";
  }
  PrimFile <<"\n]";
  return ncol;
}


bool checkVol(varStore *prm,Int_t geomId){
   bool check=prm->MapVol[geomId];
   prm->MapVol[geomId]=false;
   return check;
}

TGeoCompositeShape * Conv(TGeoCompositeShape *shape){
  return shape;
}

int ExportCurrentVolume(varStore *prm,TGeoVolume *vol, TGeoShape *shape,TGeoMatrix *matrix,TString Comp,TString CompL,TString CompR){

  TString guid, guidVol,CompC,CompName,CompNameL,CompNameR;
  TGeoCompositeShape * Cshape;

  if (shape->IsComposite()) {
    CompC=Comp;
    Cshape=Conv((TGeoCompositeShape *)shape);

    CompL=Comp+"-L";
    CompC=CompL;

    ExportCurrentVolume(prm,vol,Cshape->GetBoolNode()->GetLeftShape(),Cshape->GetBoolNode()->GetLeftMatrix(),CompC,CompL,CompR);

    CompR=Comp+"-R";
    CompC=CompR;
    ExportCurrentVolume(prm,vol,Cshape->GetBoolNode()->GetRightShape(),Cshape->GetBoolNode()->GetRightMatrix(),CompC,CompL,CompR);

    CompName =prm->current->GetVolume()->GetName()+Comp;
    CompNameL=prm->current->GetVolume()->GetName()+CompL;
    CompNameR=prm->current->GetVolume()->GetName()+CompR;
    TString op=" + ";
    TGeoBoolNode::EGeoBoolType boolType = ((TGeoCompositeShape *)shape)->GetBoolNode()->GetBooleanOperator();
    TString operation="";
    switch (boolType) {
      case TGeoBoolNode::kGeoUnion:
        operation="Union";
        break;
      case TGeoBoolNode::kGeoSubtraction:
        operation="Subtraction";
        op=" - ";
        break;
      case TGeoBoolNode::kGeoIntersection:
        operation="Intersection";
        op=" ^ ";
        break;
    }
    // 
    getMatrixVR(prm,matrix);
    ExportJsonComposite(prm,CompNameL,CompNameR,operation,CompName);
    //cout << "\n  COMP!!\t"<<prm->current->GetNumber()<<"\t"<<vol->GetNumber()<<"\t"<<CompName<<" = "<<" ["<<CompNameL<<op<<CompNameR<<"] ";

  } else {

    getMatrixVR(prm,matrix);
    if (Comp==""){
        TString clsname=shape->ClassName();
        exportVolGeometry(prm,shape,vol->GetNumber());
      	if (clsname != "TGeoShapeAssembly" )writeMesh(prm,vol->GetNumber(),gGeoManager->CountNodes(vol,1000,1));
        //cout << "\n  norm!!\t"<<prm->current->GetNumber()<<"\t"<<vol->GetNumber()<<"\t"<<vol->GetName()<<Comp;
        //cout << "\n  NORM!!\t"<<prm->current->GetNumber()<<"\t"<<vol->GetNumber()<<"\t"<<vol->GetName()<<Comp;
    } else {
        prm->nMaxVol+=1;
        //cout << "\n  COMP!!\t"<<prm->current->GetNumber()<<"\t"<<prm->nMaxVol<<"\t"<<vol->GetName()<<Comp;
        exportVolGeometry(prm,shape,prm->nMaxVol);
      	writeCompMesh(prm,prm->nMaxVol,Comp);
    }

  }
  return prm->kVolume;
}


void initClonesFile(varStore *prm){
  prm->ClonesFile.open(prm->ClonesFilename);
  prm->ClonesFile<< "{\n\"metadata\": { \"version\": 1.0,\"type\": \"Object\", \"generator\": \"CogevitoExporter\",\"File\":\""<<prm->PrimFilename<<"\"},";

}
void closeClonesFile(varStore *prm){
  prm->ClonesFile << "\n]}";
  prm->ClonesFile.close();
}

void initPrimFile(varStore *prm){
  prm->PrimFile.open(prm->PrimFilename);
  prm->PrimFile<< "{\n\"metadata\": { \"version\": 1.0,\"type\": \"Object\", \"generator\": \"CogevitoExporter\",\"File\":\""<<prm->PrimFilename<<"\"},";
}

void closePrimFile(varStore *prm){
  prm->PrimFile << "\n]}";
  prm->PrimFile.close();
}


int ExportPhyVols(varStore *prm){

  TGeoNode *current;
  TGeoVolume *vol;
  TGeoShape *shape;
  prm->level=0;

  initPrimFile(prm);  
  initClonesFile(prm);

  TString clsname,volname;

  prm->nColors=defJSONcolors(prm->PrimFile);

  ExportStart(prm);

  TGeoIterator iter(prm->top);
  int nvs=1;

  while ((current = iter.Next())) {

    prm->level=iter.GetLevel();
    prm->current=current;
    vol=current->GetVolume();
    shape=vol->GetShape();
    clsname  = shape->ClassName();

    if (prm->level <= prm->MaxVisiLevel) {
      if ( checkVol(prm,vol->GetNumber()) ) {
        //if (current->GetNumber()>=2) cout <<"\nWARNING copy number: "<<vol->GetNumber();
        ExportCurrentVolume(prm,vol,shape,current->GetMatrix(),"","","");
        //cout <<"\t\t   --->"<<nvs++;
        volname=vol->GetName();
        //if (volname=="SC5D") break;
      } else {
        //if (current->GetNumber()<2) cout <<"\nWARNING copy number: "<<vol->GetNumber();
        getMatrixVR(prm,current->GetMatrix());
        cloneVolume1(prm,gGeoManager->CountNodes(vol,1000,1));
        //cout << "\n  Clnd!!\t"<<current->GetNumber()<<"\t"<<vol->GetNumber()<<"\t"<<vol->GetName()<<"\t"<<clsname;
        iter.Skip();
      }
      prm->previousLevel=prm->level;
    } else {
      iter.Skip();
    }
  }
  cout <<"\n"<<gGeoManager->CountNodes(prm->top,1000,1)<<" Nodes have been defined in the file: "<< prm->PrimFilename<<"\n";

  ExportEnd(prm);
  closeClonesFile(prm);
  closePrimFile(prm);

  return(prm->kVolume);
}


int jsontools(){
  return(0);
}

#endif
