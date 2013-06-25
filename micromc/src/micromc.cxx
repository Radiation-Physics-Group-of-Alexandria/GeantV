#include <TFile.h>
#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TList.h>
#include <TRandom.h>
#include <GeantTrack.h>

#include <TPartIndex.h>
#include <TMXsec.h>

void GenerateEvent(Double_t avemult, Double_t energy, Double_t fVertex[3]);
Double_t SampleMaxwell(Double_t emean);

static TList particleStack;

Int_t main (int argc, char *argv[]) {

   for(Int_t i=0; i<argc; ++i) {
      printf("argv[%d] = %s\n",i,argv[i]);
   }

   Int_t nevent=1;
   if(argc>1) sscanf(argv[1],"%d",&nevent);

   Double_t avemult = 10.;
   if(argc>2) sscanf(argv[2],"%lf",&avemult);
   
   Double_t energy = 10.;
   if(argc>3) sscanf(argv[3],"%lf",&energy);
		 
   printf("Generating %d events with ave multiplicity %f and energy %f\n",nevent,avemult, energy);

   const Char_t *geofile="http://root.cern.ch/files/cms.root";
   TGeoManager *geom = TGeoManager::Import(geofile);
   
   // loop materials

   TFile *f = new TFile("xsec.root");
   TList *matlist = (TList*) geom->GetListOfMaterials();
   TIter next(matlist);
   TGeoMaterial *mat=0;
   TGeoMixture *mix=0;
   Int_t nmater = matlist->GetEntries();
   while((mat = (TGeoMaterial*) next())) {
      if(!mat->IsUsed()) continue;
      Int_t nelem = mat->GetNelements();
      Int_t *z = new Int_t[nelem];
      Int_t *a = new Int_t[nelem];
      Float_t *w = new Float_t[nelem];
      for(Int_t iel=0; iel<nelem; ++iel) {
	 Double_t ad;
	 Double_t zd;
	 Double_t wd;
	 mat->GetElementProp(ad,zd,wd,iel);
	 a[iel]=ad;
	 z[iel]=zd;
	 w[iel]=wd;
	 printf("Mixture %s element %s z %d a %d\n",
		mat->GetName(), mat->GetElement(iel)->GetName(),
		z[iel],a[iel]);
      }
      mat->SetFWExtension(
	new TGeoRCExtension(
	   new TMXsec(mat->GetName(),mat->GetTitle(),
	   z,a,w,nelem,mat->GetDensity(),kTRUE))); 
      //      myObject = mat->GetExtension()->GetUserObject();
      delete [] a;
      delete [] z;
      delete [] w;
   }

   for(Int_t iev=0; iev<nevent; ++iev) {
      // should define a vertex, origin for the moment
      Double_t vertex[3]={0,0,0};
      GenerateEvent(avemult, energy, vertex);
   }
   /*
   Double_t dir[3];
   Double_t pos[3];
   TGeoNode *current=0;
   TGeoNode *nexnode=0;
   TIter next(particleStack);
   for(Int_t iev=0; iev<nevent; ++iev) {
      GenerateEvent();
      next.Reset();
      GeantTrack *tr=0;
      while((tr=(GeantTrack*)next())) {
	 Int_t G5index = TPartIndex::I()->PartIndex(tr->pdg);
	 tr->Direction(dir);
	 x[0]=tr->xpos;
	 x[1]=tr->ypos;
	 x[2]=tr->zpos;
	 // where am I
	 current = geom->InitTrack(pos,dir);
	 Double_t ken = 
	 while(!geom->IsOutside()) {
	    mat = current->GetVolume()->GetMaterial();
	    Double_t xlen = mat->GetUserField()->Xlength(G5index,;
	    nexnode = geom->FindNextBoundaryAndStep(xlen);
	    Double_t snext = geom->GetStep();
	    if(snext>xlen) {
	       //phys wins
	    } else {
	       // geom wins
	       dirnew = something;
	       geom->SetCurrentDirection(dirnew);
		  
	    }
	 }
	 
	    
      }
   */
   return 0;
}

#define NPART 11

void GenerateEvent(Double_t avemult, Double_t energy, Double_t fVertex[3]) {
   static Bool_t first=kTRUE;
   static const Int_t kMaxPart=NPART;
   static const Char_t* G5name[NPART] = {"pi+","pi-","proton","antiproton","neutron","antineutron","e-","e+","gamma"
					  "mu+","mu-"};
   static const Species_t G5species[NPART] = {kHadron, kHadron, kHadron, kHadron, kHadron, kHadron, 
					      kLepton, kLepton, kLepton, kLepton, kLepton};
   static Int_t G5part[NPART];
   static Float_t G5prob[NPART] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

   const Double_t etamin = -3, etamax = 3;

   // Initialise simple generator
   if(first) {
      particleStack.SetOwner();
      Double_t sumprob=0;
      for(Int_t ip=0; ip<kMaxPart; ++ip) {
	 G5part[ip] = TPartIndex::I()->PartIndex(G5name[ip]);
	 printf("part %s code %d\n",G5name[ip],G5part[ip]);
	 sumprob += G5prob[ip];
      }
      for(Int_t ip=0; ip<kMaxPart; ++ip) {
	 G5prob[ip]/=sumprob;
	 if(ip) G5prob[ip]+=G5prob[ip-1];
      }
      first=kFALSE;
   }
   
   particleStack.Clear();
   Int_t ntracks = gRandom->Poisson(avemult)+0.5;
   for (Int_t i=0; i<ntracks; i++) {
      GeantTrack *track=new GeantTrack();
      Double_t prob = gRandom->Uniform();
      for(Int_t j=0; j<kMaxPart; ++j) {
	 if(prob <= G5prob[j]) {
	    track->fG5code = G5part[j];
	    track->pdg = TPartIndex::I()->PDG(G5part[j]);
	    track->species = G5species[j];
	    printf("Generating a %s",TDatabasePDG::Instance()->GetParticle(track->pdg)->GetName());
	    //	    pdgCount[j]++;
	    break;
	 }
      }   
      if(!track->pdg) Fatal("ImportTracks","No particle generated!");
      TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track->pdg);
      track->charge = part->Charge()/3.;
      track->mass   = part->Mass();
      track->xpos = fVertex[0];
      track->ypos = fVertex[1];
      track->zpos = fVertex[2];
      Double_t ekin = SampleMaxwell(energy/(avemult*1.5));
      track->e = ekin+track->mass;
      Double_t p = TMath::Sqrt(ekin*(2*ekin+track->mass));
      Double_t eta = gRandom->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
      //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
      Double_t phi = TMath::TwoPi()*gRandom->Rndm();
      track->px = p*TMath::Sin(theta)*TMath::Cos(phi);
      track->py = p*TMath::Sin(theta)*TMath::Sin(phi);
      track->pz = p*TMath::Cos(theta);
      track->frombdr = kFALSE;
      Int_t itrack = track->particle;
      
      particleStack.Add(track);
   }
//      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
}

Double_t SampleMaxwell(Double_t emean) 
{
   Double_t th = gRandom->Uniform()*TMath::TwoPi();
   Double_t rho = TMath::Sqrt(-TMath::Log(gRandom->Uniform()));
   Double_t mx = rho*TMath::Sin(th);
   return emean*(-TMath::Log(gRandom->Uniform())+mx*mx);
}

