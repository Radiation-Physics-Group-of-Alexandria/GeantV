
TString getColor(varStore *prm){
  //Int_t kC=prm->current->GetColour();
  Int_t kC=prm->current->GetVolume()->GetLineColor( ) ;
  TColor *color = gROOT->GetColor(kC);
  TString ColorHex=color->AsHexString();
  ColorHex.Remove(0,1);

  TString ColorString="\"#"+ColorHex+"\",\"cn\":\""+color->GetName()+"\"";
  return ColorString;
}


TString getMother(varStore *prm){
  if (prm->level<=1) {
    return "root";
  } else {
    return prm->current->GetMotherVolume()->GetName();
  }
}

TString getVolName(varStore *prm){
  //return prm->current->GetVolume()->GetName()+TString::Format("-%d",prm->current->GetNumber());
  return prm->current->GetVolume()->GetName();
}

TString getCloneVolName(varStore *prm){
  return prm->current->GetVolume()->GetName()+TString::Format("-%d",prm->current->GetNumber());
  //return prm->current->GetVolume()->GetName());
}

int getVisible(varStore *prm){
  if (prm->level==1){
    return false;
  }
  return prm->current->IsVisible();

}

void ExportEnd(varStore *prm){
  prm->PrimFile << ",\n\t{\"fn\":\"END\"}";
  prm->ClonesFile << ",\n\t{\"fn\":\"END\"}";
}

void ExportStart(varStore *prm){
  prm->PrimFile<<",\n\"Geometry\":[";
  prm->PrimFile << "\n\t{\"fn\":\"START\"}";
  prm->ClonesFile<<"\n\"Clones\":[";
  prm->ClonesFile << "\n\t{\"fn\":\"START\"}";
}


void writeMesh(varStore *prm,Int_t geomId,Int_t nv){
  prm->PrimFile << ",\n\t{\"fn\":\"mesh\",\"gid\":"<<geomId<<",\"nv\":"<<nv<<",\"lv\":"<<prm->level <<",\"m\":\""<<getMother(prm) <<"\",\"cl\":" << getColor(prm)  <<",\"vs\":"<<getVisible(prm) <<",\"n\":\"" <<getVolName(prm)<<"\",\"v\":[";
  for (int i=0;i<prm->nM;i++) {
    prm->PrimFile <<prm->Mtrx[i];
    if (i<(prm->nM-1))prm->PrimFile <<",";
  }; 
  prm->PrimFile << "]}";
}

void writeCompMesh(varStore *prm,Int_t geomId,TString Comp){
  prm->PrimFile << ",\n\t{\"fn\":\"mesh\",\"gid\":"<<geomId<<",\"nv\":0,\"lv\":"<<prm->level <<",\"m\":\"root\",\"cl\":" << getColor(prm)  <<",\"vs\":"<<getVisible(prm) <<",\"n\":\"" <<getVolName(prm) <<Comp<<"\",\"v\":[";
  for (int i=0;i<prm->nM;i++) {
    prm->PrimFile <<prm->Mtrx[i];
    if (i<(prm->nM-1))prm->PrimFile <<",";
  }; 
  prm->PrimFile << "]}";
}


void writeGeomObject(varStore *prm,Int_t geomId,TString Vtype){
  prm->PrimFile << ",\n\t{\"fn\":\""<<Vtype<<"\",\"gid\":"<<geomId<<",";
}

void ExportGeomAssembly(varStore *prm){
  int nv = gGeoManager->CountNodes(prm->current->GetVolume(),1000,1);
  prm->PrimFile << ",\n\t{\"fn\":\"grp\",\"nv\":"<<nv<<",\"lv\":"<<prm->level <<",\"m\":\""<<getMother(prm)<<"\",\"vs\":"<<getVisible(prm)  <<",\"n\":\"" <<getVolName(prm)<<"\",\"v\":[";
  for (int i=0;i<prm->nM;i++) {
    prm->PrimFile <<prm->Mtrx[i];
    if (i<(prm->nM-1))prm->PrimFile <<",";
  }; 
  prm->PrimFile << "]}";

}


void ExportGeomBBox(varStore *prm,Int_t geomId,   Double_t  DX, Double_t  DY, Double_t  DZ){       
  // create a cube geometry in Three-JS + material + mesh
  writeGeomObject(prm,geomId,"box");
  prm->PrimFile <<"\"x\":"<<2*DX<<",\"y\":"<<2*DY<<",\"z\":"<<2*DZ<<"}";
}


void ExportGeomConeSeg2(varStore *prm,Int_t geomId,int nedg, int nZ,  Double_t PerimR[], Double_t PerimZ[], Double_t  phiStart=0, Double_t  phiLength=2*TMath::Pi()){
  // create a lathe geometry in Three-JS + material + mesh
  writeGeomObject(prm,geomId,"pgn");

  prm->PrimFile <<"\"r\":[";
  for (int i=0;i<nZ;i++) prm->PrimFile <<PerimR[i]<<",";
  prm->PrimFile << PerimR[nZ]<<"],\"z\":[";
  for (int i=0;i<nZ;i++) prm->PrimFile <<PerimZ[i]<<",";
  prm->PrimFile << PerimZ[nZ]<<"],";
  prm->PrimFile <<"\"nedg\":"<<nedg<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";
}

void ExportGeomConeSeg1(varStore *prm,Int_t geomId, int nZ,  Double_t PerimR[], Double_t PerimZ[], Double_t  phiStart=0, Double_t  phiLength=2*TMath::Pi()){
  // create a lathe geometry in Three-JS + material + mesh
  writeGeomObject(prm,geomId,"lth");

  prm->PrimFile <<"\"r\":[";
  for (int i=0;i<nZ;i++) prm->PrimFile <<PerimR[i]<<",";
  prm->PrimFile << PerimR[nZ]<<"],\"z\":[";
  for (int i=0;i<nZ;i++) prm->PrimFile <<PerimZ[i]<<",";
  prm->PrimFile << PerimZ[nZ]<<"],";
  prm->PrimFile <<"\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";
}

void ExportGeomConeSeg(varStore *prm,Int_t geomId, Double_t  Dz, Double_t  Rmin1, Double_t  Rmin2, Double_t  Rmax1, Double_t  Rmax2 , Double_t phiStart=0.,Double_t phiLength=2*TMath::Pi()){

  writeGeomObject(prm,geomId,"csg");
  prm->PrimFile <<"\"dz\":"<<Dz<<",\"rmin1\":"<<Rmin1<<",\"rmin2\":"<<Rmin2<<",\"rmax1\":"<<Rmax1<<",\"rmax2\":"<<Rmax2<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";
}

void ExportGeomTorus(varStore *prm,Int_t geomId,  Double_t R, Double_t  Rmin, Double_t Rmax, Double_t phiStart=0.,Double_t phiLength=2*TMath::Pi() ){

  writeGeomObject(prm,geomId,"tor");
  prm->PrimFile <<"\"r\":"<<R<<",\"rmin\":"<<Rmin<<",\"rmax\":"<<Rmax<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";
}

void ExportGeomParabloid(varStore *prm,Int_t geomId, Double_t  Rlo, Double_t  Rhi, Double_t  Dz){

  writeGeomObject(prm,geomId,"prb");
  prm->PrimFile <<"\"rlo\":"<<Rlo<<",\"rhi\":"<<Rhi<<",\"dz\":"<<Dz<<"}";
}

void ExportGeomHype(varStore *prm,Int_t geomId, Double_t  Rmin, Double_t  Rmax,Double_t  StIn, Double_t  StOut, Double_t  Dz){

  StIn=TMath::DegToRad()*StIn;
  StOut=TMath::DegToRad()*StOut;

  writeGeomObject(prm,geomId,"hpb");
  prm->PrimFile <<"\"rmin\":"<<Rmin<<",\"rmax\":"<<Rmax<<",\"stin\":"<<StIn<<",\"stout\":"<<StOut<<",\"dz\":"<<Dz<<"}";
}


void ExportGeomExtrude(varStore *prm,Int_t geomId,  int np,int nz,Double_t x[],Double_t y[],Double_t z[],Double_t px[],Double_t py[],Double_t f[]){

  writeGeomObject(prm,geomId,"exr");

  prm->PrimFile <<"\"x\":[";
  for (int i=0;i<np;i++) prm->PrimFile<<x[i]<<",";
  prm->PrimFile<<x[0]<<"],";

  prm->PrimFile <<"\"y\":[";
  for (int i=0;i<np;i++) prm->PrimFile<<y[i]<<",";
  prm->PrimFile<<y[0]<<"],";

  prm->PrimFile <<"\"z\":[";
  for (int i=0;i<nz-1;i++) prm->PrimFile<<z[i]<<",";
  prm->PrimFile<<z[nz-1]<<"],";

  prm->PrimFile <<"\"px\":[";
  for (int i=0;i<nz-1;i++) prm->PrimFile<<px[i]<<",";
  prm->PrimFile<<px[nz-1]<<"],";

  prm->PrimFile <<"\"py\":[";
  for (int i=0;i<nz-1;i++) prm->PrimFile<<py[i]<<",";
  prm->PrimFile<<py[nz-1]<<"],";

  prm->PrimFile <<"\"f\":[";
  for (int i=0;i<nz-1;i++) prm->PrimFile<<f[i]<<",";
  prm->PrimFile<<f[nz-1]<<"]}";

}


void ExportGeomSphere(varStore *prm,Int_t geomId, Double_t  Rmin, Double_t  Rmax, Double_t thetaStart=0., Double_t thetaLength=TMath::Pi(), Double_t phiStart=0., Double_t phiLength=2*TMath::Pi()){

  writeGeomObject(prm,geomId,"sph");
  prm->PrimFile <<"\"rmin\":"<<Rmin<<",\"rmax\":"<<Rmax<<",\"ths\":"<<thetaStart<<",\"thl\":"<<thetaLength<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";
}


void ExportGeomQuad(varStore *prm,Int_t geomId,  int np,  Double_t x[], Double_t y[], Double_t z[], int nz=1){

  writeGeomObject(prm,geomId,"tra");

  prm->PrimFile <<"\"x\":[";
  for (int i=0;i<np-1;i++) prm->PrimFile<<x[i]<<",";
  prm->PrimFile<<x[np-1]<<"],";  
 
  prm->PrimFile <<"\"y\":[";
  for (int i=0;i<np-1;i++) prm->PrimFile<<y[i]<<",";
  prm->PrimFile<<y[np-1]<<"],";  

  prm->PrimFile <<"\"z\":[";
  for (int i=0;i<np-1;i++) prm->PrimFile<<z[i]<<",";
  prm->PrimFile<<z[np-1]<<"],";  

  prm->PrimFile <<"\"nz\":"<<nz<<"}";

}

void ExportGeomPara(varStore *prm,Int_t geomId,  Double_t dx,Double_t dy,Double_t dz,Double_t alpha,Double_t theta,Double_t phi){

  writeGeomObject(prm,geomId,"par");
  prm->PrimFile <<"\"dx\":"<<dx<<",\"dy\":"<<dy<<",\"dz\":"<<dz<<",\"al\":"<<alpha<<",\"th\":"<<theta<<",\"ph\":"<<phi<<"}";
}


void ExportGeomTrap(varStore *prm,Int_t geomId,    Double_t dz,Double_t theta,Double_t phi,Double_t bl1,Double_t tl1,Double_t bl2,Double_t tl2,Double_t h1,Double_t h2,Double_t alpha1,Double_t alpha2 ){

  writeGeomObject(prm,geomId,"trp");
  prm->PrimFile <<"\"dz\":"<<dz<<",\"th\":"<<theta<<",\"ph\":"<<phi<<",\"bl1\":"<<bl1<<",\"tl1\":"<<tl1<<",\"bl2\":"<<bl2<<",\"tl2\":"<<tl2<<",\"h1\":"<<h1<<",\"h2\":"<<h2<<",\"al1\":"<<alpha1<<",\"al2\":"<<alpha2<<"}";
}


void ExportGeomGTra(varStore *prm,Int_t geomId,    Double_t dz,Double_t theta,Double_t phi,Double_t bl1,Double_t tl1,Double_t bl2,Double_t tl2,Double_t h1,Double_t h2,Double_t alpha1,Double_t alpha2,Double_t twist ){

  writeGeomObject(prm,geomId,"gtr");
  prm->PrimFile <<"\"dz\":"<<dz<<",\"th\":"<<theta<<",\"ph\":"<<phi<<",\"bl1\":"<<bl1<<",\"tl1\":"<<tl1<<",\"bl2\":"<<bl2<<",\"tl2\":"<<tl2<<",\"h1\":"<<h1<<",\"h2\":"<<h2<<",\"al1\":"<<alpha1<<",\"al2\":"<<alpha2<<",\"twist\":"<<twist<<"}";
}


void ExportGeomTrd(varStore *prm,Int_t geomId,    Double_t dx1,  Double_t dx2,  Double_t dy1,  Double_t dy2,  Double_t dz){

  writeGeomObject(prm,geomId,"trd");
  prm->PrimFile <<"\"dx1\":"<<dx1<<",\"dx2\":"<<dx2<<",\"dy1\":"<<dy1<<",\"dy2\":"<<dy2<<",\"dz\":"<<dz<<"}";
}


void ExportGeomTub(varStore *prm,Int_t geomId,  Double_t Rmin,Double_t Rmax,Double_t Dz, Double_t phiStart=0., Double_t phiLength=2*TMath::Pi()){

//  if (Rmin==0) {

//    writeGeomObject(prm,geomId,"cyl");
//    prm->PrimFile <<"\"ncs\":"<<prm->ncseg<<",\"nzs\":"<<prm->nzseg<<",\"dz\":"<<2*Dz<<",\"rtop\":"<<Rmax<<",\"rbot\":"<<Rmax<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";

//  } else {

    writeGeomObject(prm,geomId,"csg");
    prm->PrimFile <<"\"dz\":"<<Dz<<",\"rmin1\":"<<Rmin<<",\"rmin2\":"<<Rmin<<",\"rmax1\":"<<Rmax<<",\"rmax2\":"<<Rmax<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";
//  }

}


void ExportGeomEltub(varStore *prm,Int_t geomId,  Double_t A,Double_t B,Double_t Dz, Double_t phiStart=0., Double_t phiLength=2*TMath::Pi()){

    writeGeomObject(prm,geomId,"elc");
    prm->PrimFile <<"\"dz\":"<<2*Dz<<",\"a\":"<<A<<",\"b\":"<<B<<",\"phs\":"<<phiStart<<",\"phl\":"<<phiLength<<"}";

}


void ExportJsonComposite(varStore *prm,TString CompNameL,TString CompNameR,TString operation,TString CompName){

  prm->PrimFile << ",\n\t{\"fn\":\"cmp\",\"nv\":"<<gGeoManager->CountNodes(prm->current->GetVolume(),1000,1)<<",\"lv\":"<<prm->level <<",\"m\":\""<<getMother(prm) <<"\",\"cl\":" << getColor(prm)  <<",\"vs\":"<<getVisible(prm) <<",\"v\":[";
  for (int i=0;i<prm->nM;i++) {
    prm->PrimFile <<prm->Mtrx[i];
    if (i<(prm->nM-1))prm->PrimFile <<",";
  }; 
  prm->PrimFile << "],\"lpn\":\""<<CompNameL<<"\",\"rpn\":\""<<CompNameR<<"\",\"opr\":\""<<operation<<"\",\"cpn\":\""<<CompName<<"\"}";
}


void cloneVolume1(varStore *prm,Int_t nv){
 prm->ClonesFile << ",\n\t{\"fn\":\"clone\",\"nv\":"<<nv<<",\"lv\":"<<prm->level  <<",\"m\":\""<<getMother(prm) <<"\",\"ovn\":\"" << getVolName(prm) <<"\",\"cn\":"<<prm->current->GetNumber() <<",\"v\":[";
  for (int i=0;i<prm->nM;i++) {
    prm->ClonesFile <<prm->Mtrx[i];
    if (i<(prm->nM-1))prm->ClonesFile <<",";
  }; 
  prm->ClonesFile << "]}";
}

void cloneVolume(varStore *prm,Int_t nv){
 prm->PrimFile   << ",\n\t{\"fn\":\"clone\",\"nv\":"<<nv<<",\"lv\":"<<prm->level  <<",\"m\":\""<<getMother(prm) <<"\",\"ovn\":\"" << getVolName(prm) <<"\",\"cn\":"<<prm->current->GetNumber() <<",\"v\":[";
  for (int i=0;i<prm->nM;i++) {
    prm->PrimFile <<prm->Mtrx[i];
    if (i<(prm->nM-1))prm->PrimFile <<",";
  }; 
  prm->PrimFile << "]}";
}



int geomtools(){
  return 0;
}

