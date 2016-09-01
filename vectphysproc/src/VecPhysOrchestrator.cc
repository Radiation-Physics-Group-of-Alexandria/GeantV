/*
 Author: M. Bandieramonte
 */
#include "VecPhysOrchestrator.h"

#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "globals.h"
#include "GeantTrack.h"
#include "Geant/Error.h"

//mb:start
#include "base/Global.h"
#include "Geant/Typedefs.h"
#include "GeantFwd.h"
#include "Geant/Config.h"
#include "TPartIndex.h"
//mb:end

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#else
#include "TGeoBranchArray.h"
#endif

#ifdef USE_ROOT
#include "Rtypes.h"
ClassImp( VecPhysOrchestrator )
#endif


//------------------------------------------------------------------------------
VecPhysOrchestrator::VecPhysOrchestrator()
: fProcessId(12), fEnergyLimit(0.), fVComptonModel(0), fTargetElements(0), fParentTrackIndices(0),
fPrimaryTracks(0), fSecondaryTracks(0) {fComptonTotTracks=0.;}

//------------------------------------------------------------------------------
VecPhysOrchestrator::VecPhysOrchestrator(int processId, double energyLimit, int numtracks)
: fProcessId(processId), fEnergyLimit(energyLimit), fPrimaryTracks(0), fSecondaryTracks(0) {
    fComptonTotTracks=0.;
    fPrimaryTracks = new GUTrack_v();
    fSecondaryTracks = new GUTrack_v();
    Allocator(numtracks);                         // allocate numtracks GUTrack_v for primaries
    // and secondaries and so on; will be increased if needed
    //fVComptonModel = new vecphys::ComptonKleinNishina(0,-1); // create the vector physics model
    fVComptonModel = new vecphys::cxx::ComptonKleinNishina(0,-1); // create the vector physics model
    //SamplingMethod sampleType=kAlias
    fVComptonModel->SetSamplingMethod(vecphys::kAlias);
    
    //static vecphys::cxx::ComptonKleinNishina model(0,-1);
    
    std::cout << "------------ processId is: = " << fProcessId << std::endl;
}

//------------------------------------------------------------------------------
VecPhysOrchestrator::~VecPhysOrchestrator() {
    Deallocator();
    delete fPrimaryTracks;
    delete fSecondaryTracks;
    delete fVComptonModel;
}

void VecPhysOrchestrator::CheckEnergyConservation(GeantTrack_v &gTrackV, int numtracks, GUTrack_v &primaries, GUTrack_v &secondaries){
    
    //std::cout<<"***CheckEnergyConservation: START.\n";
    int primariesCompton=primaries.numTracks;
    //int secondariesCompton=secondaries.numTracks;
    //std::cout<<"CheckEnergyConservation\n";
    double check;
    for (int i = 0; i < numtracks; ++i)//all the tracks
    {
        for (int j=0; j<primariesCompton; ++j)
            if(primaries.parentId[j]==i) // i: tot track, j itera solo sui Compton -> j-esima track di fPrimaries corrisponde a i-esima track originaria
            {
                //std::cout<<"CheckEnergyConservation: "<<j<<"-th Compton track\n";
                check=(gTrackV.fEV[i]-gTrackV.fMassV[i])*1000;
                check-=primaries.E[j];
                //check=check-secondaries.E[secondaries.parentId[j]];
                check=check-secondaries.E[j]; //nel caso del Compton è una corrispondenza 1-1 (only valid for COMPTON)

                //std::cout<<"CheckEnergyConservation: gTrackV.fEV["<<i<<"]-gTrackV.fMassV["<<i<<"]: "<< (gTrackV.fEV[i]-gTrackV.fMassV[i])*1000<< " \nprimaries.E["<<j<<"]: "<<primaries.E[j]<<"\nsecondaries.E["<<j<<"]: "<<secondaries.E[j]<<"\n";
                if(check!=0)
                {
                    std::cout<<"CheckEnergyConservation, failed. Exiting.\n";
                    exit(0);
                }
                //else
                //    std::cout<<"CheckEnergyConservation: OK!\n";
            }
    }
    //std::cout<<"***CheckEnergyConservation: END.\n";
}

void VecPhysOrchestrator::CheckDirectionUnitVector(GeantTrack_v &gTrackV, GUTrack_v &primaries, GUTrack_v &secondaries){
    
    int numtracks=gTrackV.fNtracks;
    int primariesCompton=primaries.numTracks;
    double checkModule, checkX, checkY, checkZ ;
    double epsilon=1.e-10;
    for (int i = 0; i < numtracks; ++i)//all the GeantV tracks
    {
        for (int j=0; j<primariesCompton; ++j)
            if(primaries.parentId[j]==i) // i: tot track, j itera solo sui Compton -> j-esima track di fPrimaries corrisponde a i-esima track originaria
            {
                checkX=gTrackV.fXdirV[i]*gTrackV.fXdirV[i];
                checkY=gTrackV.fYdirV[i]*gTrackV.fYdirV[i];
                checkZ=gTrackV.fZdirV[i]*gTrackV.fZdirV[i];
                checkModule=sqrt(checkX+checkY+checkZ);
                checkModule= abs(checkModule-1);
                if(checkModule>epsilon)
                {
                    std::cout<<"CheckDirectionUnitVector, failed. checkModule: "<<checkModule<<". Exiting.\n";
                    exit(0);
                }
            }
    }
    //std::cout<<"***CheckDirectionUnitVector: OK.\n";
}


//_____________________________________________________________________________
// (oldxdir, oldydir, oldzdir) are the direction vector of parent track in lab.
// frame; direction vector of the current track, measured from local Z is in
// (newxdir, newydir, newzdir); here we rotate it to lab. frame -> taken from TabPhysics ONE TRACK
void VecPhysOrchestrator::RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack_v &track, int index) {
    const double one = 1.0;
    const double zero = 0.0;
    const double amin = 1.0e-10;
    const double one5 = 1.5;
    const double half = 0.5;
    
    double cosTheta0 = oldZdir;
    double sinTheta0 = sqrt(oldXdir * oldXdir + oldYdir * oldYdir);
    double cosPhi0;
    double sinPhi0;
    
    if (sinTheta0 > amin) {
        cosPhi0 = oldXdir / sinTheta0;
        sinPhi0 = oldYdir / sinTheta0;
    } else {
        cosPhi0 = one;
        sinPhi0 = zero;
    }
    
    double h0 = track.fXdirV[index];
    double h1 = sinTheta0 * track.fZdirV[index] + cosTheta0 * h0;
    double h2 = track.fYdirV[index];

    track.fXdirV[index] = h1 * cosPhi0 - h2 * sinPhi0;
    track.fYdirV[index] = h1 * sinPhi0 + h2 * cosPhi0;
    track.fZdirV[index] = track.fZdirV[index] * cosTheta0 - h0 * sinTheta0;
    
    // renormalization: -use 1-th order Taylor aprx. of 1/sqrt(x) around 1.0
    // that should be almost exact since the vector almost normalized!
    double delta = one5 - half * (track.fXdirV[index] * track.fXdirV[index] + track.fYdirV[index] * track.fYdirV[index] + track.fZdirV[index] * track.fZdirV[index]);
    track.fXdirV[index] *= delta;
    track.fYdirV[index] *= delta;
    track.fZdirV[index] *= delta;
}



void VecPhysOrchestrator::DebugTracksEnergies(GeantTrack_v &gTrackV, int numtracks, GUTrack_v &primaries, bool checkIdentity)
{
    std::cout<<"\n***DebugTracksEnergies for "<<numtracks<<" tracks. START\n";
    int primariesCompton=primaries.numTracks;
    for (int i = 0; i < numtracks; ++i)
    {
        for (int j=0; j<primariesCompton; ++j)
            if(primaries.parentId[j]==i)
            {
                if(!checkIdentity)
                {
                    std::cout<<"DebugTracksEnergies after conversion:\nKin En. copied on the primaries track: fPrimaries.E["<<j<<"]:  "<<primaries.E[j];
                    std::cout<<"\nTot En. original track, gTrackV.fEV["<<i<<"]: "<<gTrackV.fEV[i];
                    std::cout<<"\nRest mass original track, gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<"\n";
                }
                    
                if((gTrackV.fEV[i]-gTrackV.fMassV[i])!=primaries.E[j] && checkIdentity)
                {
                    std::cout<<"DebugTracksEnergies:\nKin En. copied on the primaries track: fPrimaries.E["<<j<<"]:  "<<primaries.E[j];
                    std::cout<<"\nTot En. original track, gTrackV.fEV["<<i<<"]: "<<gTrackV.fEV[i];
                    std::cout<<"\nRest mass original track, gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<"\n";
                    std::cout<<"DegugTracksEnergies ERROR, copy failed. Exiting. bye bye.\n";
                    exit(0);
                }
                if(primaries.E[j]<0)
                {
                    std::cout<<"DegugTracksEnergies ERROR, negative primary energy. Exiting, bye bye.\n";
                    exit(0);
                }
                
            }
    
    }
    std::cout<<"***DegugTracksEnergies for "<<numtracks<<" tracks successfully done. END\n";
}

//------------------------------------------------------------------------------
int VecPhysOrchestrator::ApplyPostStepProcess(GeantTrack_v &gTrackV, int numtracks, int tid) {
    
    
    //std::cout<<"\n\n*** ApplyPostStepProcess. START ***\n";
    FilterPrimaryTracks(gTrackV, numtracks); //update the fPrimaryTracks copying data of the tracks undergoing Compton and store the index that the tracks had in the original gTrackV in fParentTrackIndices
    CheckDirectionUnitVector(gTrackV, *fPrimaryTracks, *fSecondaryTracks);
    if (fPrimaryTracks->numTracks == 0) // if there is no track with Compton -> return
    {
        //std::cout<<"No tracks with compton, returning.\n";
        return 0;
    }
    //std::cout<<fPrimaryTracks->numTracks<<" tracks with compton.\n";
    //DebugTracksEnergies(gTrackV, numtracks, *fPrimaryTracks, 1); //[GeV] --> no negative values
     //std::cout<<"ConvertEnergiesToVecPhys, START.\n";
    ConvertEnergiesToVecPhys(); //[GeV]->[MeV]
    //std::cout<<"ConvertEnergiesToVecPhys, END.\n";
    //DebugTracksEnergies(gTrackV, numtracks, *fPrimaryTracks, 0); // ---> still not negative after conversion
    PerformInteraction(); //call KleinNishina Compton
    CheckEnergyConservation(gTrackV, numtracks, *fPrimaryTracks, *fSecondaryTracks);
    ConvertEnergiesFromVecPhys(); //[MeV]->[GeV]
    return WriteBackTracks(gTrackV, tid);
    
    
}


//------------------------------------------------------------------------------
// This method should scan all the tracks and check the "compton" ones
// count number of tracks in GeantTrack_v with selected process = Compton

void VecPhysOrchestrator::FilterPrimaryTracks(GeantTrack_v &gTrackV, int numtracks)
{
    //std::cout<<"***FilterPrimaryTracks: START.\n";
    int numComptonTracks = 0;   //n. of tracks undergoing Compton
    int totNumTracks = 0;       //total n. of tracks - undergoing physics, it should be equal to numtracks. To check
    for (int i = 0; i < numtracks; ++i)
    {
        ++totNumTracks;
        //std::cout<<"DEBUG ALL: Mass - gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<" and  Process - gTrackV.fProcessV["<<i<<"]: "<<gTrackV.fProcessV[i]<<"\n";
        if (gTrackV.fProcessV[i] == fProcessId)
            ++numComptonTracks;
    }
    fComptonTotTracks+=numComptonTracks;
    if (numComptonTracks < 1) {
        fPrimaryTracks->numTracks = 0;
        return;
    }
    
    // check if we have enough space: if not make sure that we have
    if (numComptonTracks > fPrimaryTracks->capacity)
    {
        std::cout<<"Need to allocate more space.\n";
        Deallocator();
        Allocator(numComptonTracks);
    }
     // form the input GUTrack_v with the corresponding primary track members
     fPrimaryTracks->numTracks = numComptonTracks;
    
    //double positiveEnergies[numComptonTracks];
    int j = 0;
    double momentum;
    for (int i = 0; i < numtracks && j < numComptonTracks; ++i) { //for all the tracks undergoing physics and only a n. times equal to the numComptonTracks
        if (gTrackV.fProcessV[i] == fProcessId) {
            fParentTrackIndices[j] = i;  // store index of this track in GeantTrack_v to store it back after the interaction
            momentum = gTrackV.fPV[i];   // total momentum
            fPrimaryTracks->px[j] = momentum * gTrackV.fXdirV[i];   // 3-momentum (px, py, pz)
            fPrimaryTracks->py[j] = momentum * gTrackV.fYdirV[i];
            fPrimaryTracks->pz[j] = momentum * gTrackV.fZdirV[i];
            fPrimaryTracks->E[j] = gTrackV.fEV[i] - gTrackV.fMassV[i] ; // Kinetic energy (NB: RestMass for gamma=0)
            if(gTrackV.fMassV[i]!=0) {std::cout<<"Error! Rest mass for gamma not equals to zero. Exiting.\n"; exit(0);}
            
            //std::cout<<"DEBUG ONLY COMPTON: gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<" and  gTrackV.fProcessV["<<i<<"]: "<<gTrackV.fProcessV[i]<<", track: "<<gTrackV.fParticleV[i]<<" with fPrimaryTracks->E["<<j<<"]: "<<fPrimaryTracks->E[j]<< " and gTrackV.fEV["<<i<<"]: "<<gTrackV.fEV[i]<<" \n";
            
            //positiveEnergies[j]=fPrimaryTracks->E[j];
            fTargetElements[j] = gTrackV.fEindexV[i]; // Z of the target atom
            fPrimaryTracks->id[j] = j; //fPrimaryTracks id is j<=i - need to initialize the id on the PrimaryTracks because it will be used in the interact method when storing the secondaries. Quite useless.
            fPrimaryTracks->parentId[j] = i; //fPrimaryTracks parentId is the same as fParentTrackIndices
            ++j;
        }
    }
    /*std::cout<<"fPrimaryTracks Ready and:\n";
    
    for(int k=0; k<numRealTracks; k++)
    {
        std::cout<<"fPrimaryTracks->E["<<k<<"]: "<<fPrimaryTracks->E[k]<<"\n";
        std::cout<<"energiePositive["<<k<<"]: "<<energiePositive[k]<<"\n";
    }
   
    std::cout<<"***FilterPrimaryTracks: END.\n";*/
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::Allocator(int size) {
    if (size < 1)
        size = 1;
    fTargetElements = new int[size];
    fParentTrackIndices = new int[size];
    GUTrackAllocator(*fPrimaryTracks, size);
    GUTrackAllocator(*fSecondaryTracks, size);
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::Deallocator() {
    if (fTargetElements)
        delete fTargetElements;
    fTargetElements = 0;
    if (fParentTrackIndices)
        delete fParentTrackIndices;
    fParentTrackIndices = 0;
    if (fPrimaryTracks)
        GUTrackDeallocator(*fPrimaryTracks);
    if (fSecondaryTracks)
        GUTrackDeallocator(*fSecondaryTracks);
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::GUTrackAllocator(GUTrack_v &gutrack_v, int size) {
    if (gutrack_v.capacity > 0)
        GUTrackDeallocator(gutrack_v);
    
    // keep only those members that are really necessary
    gutrack_v.capacity = size; // maximum number of tracks that can be stored
    gutrack_v.numTracks = 0;   // real number of tracks stored
    //   gutrack_v.status        = new int[size];  // status of the tarcks stored: (we don't need this now) ???
    //   gutrack_v.particleType  = new int[size];  // type of the particles stored: what is this exactly ?
    gutrack_v.id       = new int[size];  // ??
    gutrack_v.parentId = new int[size]; // corresponding parent index in another GUTrack_v
    //   gutrack_v.proc          = new int[size];  // process index (we don't need this now) ???
    //   gutrack_v.x             = new double[size];   // (x,y,z) position (we don't need this now) ???
    //   gutrack_v.y             = new double[size];
    //   gutrack_v.z             = new double[size];
    gutrack_v.px = new double[size]; // momentum (px,py,pz)
    gutrack_v.py = new double[size];
    gutrack_v.pz = new double[size];
    gutrack_v.E =  new double[size]; // total energy (kinetic E would be good as well)
    //   gutrack_v.q             = new double[size];   // charge (do we need this now ???)
    //   gutrack_v.s             = new double[size];   // current step length (do we need this now ???)
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::GUTrackDeallocator(GUTrack_v &gutrack_v) {
    // keep only those members that are really necessary
    //   delete gutrack_v.status;
    //   delete gutrack_v.particleType;
    delete gutrack_v.id;
    delete gutrack_v.parentId;
    //   delete gutrack_v.proc;
    //   delete gutrack_v.x;
    //   delete gutrack_v.y;
    //   delete gutrack_v.z;
    delete gutrack_v.px;
    delete gutrack_v.py;
    delete gutrack_v.pz;
    delete gutrack_v.E;
    //   delete gutrack_v.q;
    //   delete gutrack_v.s;
    
    gutrack_v.capacity = 0;
}

void VecPhysOrchestrator::DebugPrimaryAndSecondaryTracks(){
    
    //std::cout<<"***DebugPrimaryAndSecondaryTracks: START.\n";
    int numPrimaryTracks=fPrimaryTracks->numTracks;
    for(int i=0; i<numPrimaryTracks; ++i)
        std::cout<<"fPrimaryTracks->E["<<i<< "]: "<<fPrimaryTracks->E[i]<<"\n";
    
    int numSecondaryTracks=fSecondaryTracks->numTracks;
    for(int i=0; i<numSecondaryTracks; ++i)
    std::cout<<"fSecondaryTracks->E["<<i<< "]: "<<fSecondaryTracks->E[i]<<"\n";

    std::cout<<"***DebugPrimaryAndSecondaryTracks: END.\n";

}

void VecPhysOrchestrator::PerformInteraction() {
    //DebugPrimaryAndSecondaryTracks();
    //mb: primafVComptonModel->Interact<vecCore::backend::VcVector>(*fPrimaryTracks, fTargetElements, *fSecondaryTracks);
    fVComptonModel->Interact<vecCore::backend::VcVector>(*fPrimaryTracks, fTargetElements, *fSecondaryTracks);
    int nTracks=fPrimaryTracks->numTracks;
    for(int i=0; i<nTracks; ++i )
    {
        if(fPrimaryTracks->E[i]<0 ||fSecondaryTracks->E[i]<0 )
        {
            std::cout<<"fPrimaryTracks->E["<<i<<"]: "<<fPrimaryTracks->E[i]<<" and fSecondaryTracks->E["<<i<<"]: "<<fSecondaryTracks->E[i]<<"\n";
            exit(0);
        }
    }
    //model.Interact<backend::VcVector>(itrack_soa, targetElements, otrack_soa);
    //DebugPrimaryAndSecondaryTracks();
}




//------------------------------------------------------------------------------
void VecPhysOrchestrator::FilterTracksForTabPhys(GeantTrack_v &gTrackV, GeantTrack_v &gTabulatedPhysicsTracks, int numtracksIn, int &numtracksOut, int  *parentTrackIndices)
{

    int numRealTracks = 0, totNumTracks = 0;
    for (int i = 0; i < numtracksIn; ++i)
    {
        ++totNumTracks;
        if (gTrackV.fProcessV[i] != fProcessId)
        {
            ++numRealTracks;
        }
        
    }
    numtracksOut = numRealTracks;
    int j = 0;
    for (int i = 0; i < numtracksIn && j < numRealTracks; ++i) {
        if (gTrackV.fProcessV[i] != fProcessId) {
            parentTrackIndices[j] = i;                  // store index of this track in GeantTrack_v, need to update back
            gTabulatedPhysicsTracks.fPV[j]= gTrackV.fPV[i];          // total momentum
            gTabulatedPhysicsTracks.fXdirV[j] = gTrackV.fXdirV[i];   // 3-momentum
            gTabulatedPhysicsTracks.fYdirV[j] = gTrackV.fYdirV[i];
            gTabulatedPhysicsTracks.fZdirV[j] = gTrackV.fZdirV[i];
            gTabulatedPhysicsTracks.fEV[j] = gTrackV.fEV[i];                  // total energy
            gTabulatedPhysicsTracks.fEindexV[j] = gTrackV.fEindexV[i];        // Z of the target atom
            ++j;
        }
    }
}

void VecPhysOrchestrator::ConvertEnergiesToVecPhys(){

    /*
     GeantV unit is [GeV] -> [1 MeV] = 10^-3
     VecPhys unit is [MeV] -> [1 MeV] = 1
     To call VecPhys we need to convert the energies from GeV to MeV
     */
    int index=fPrimaryTracks->numTracks; //contains the number of tracks undergoing Compton
    for (int i=0; i<index; ++i)
    {
        fPrimaryTracks->E[i]*=1000;
        fPrimaryTracks->px[i]*=1000;
        fPrimaryTracks->py[i]*=1000;
        fPrimaryTracks->pz[i]*=1000;
    }
}
void VecPhysOrchestrator::ConvertEnergiesFromVecPhys(){
    /*
     GeantV unit is [GeV] -> [1 MeV] = 10^-3
     VecPhys unit is [MeV] -> [1 MeV] = 1
     To go back from VecPhys to GeantV we need to convert the energies from MeV to GeV
     */
    double unitConv=1./1000;
    int indexPrimary=fPrimaryTracks->numTracks;
    for (int i=0; i<indexPrimary; ++i)
    {
        fPrimaryTracks->E[i]*=unitConv;
        fPrimaryTracks->px[i]*=unitConv;
        fPrimaryTracks->py[i]*=unitConv;
        fPrimaryTracks->pz[i]*=unitConv;
    }
    int indexSecondary=fSecondaryTracks->numTracks;
    for (int i=0; i<indexSecondary; ++i)
    {
        fSecondaryTracks->E[i]*=unitConv;
        fSecondaryTracks->px[i]*=unitConv;
        fSecondaryTracks->py[i]*=unitConv;
        fSecondaryTracks->pz[i]*=unitConv;
    }
}

int VecPhysOrchestrator::WriteBackTracks(GeantTrack_v &gTrackV, int tid) {
    //At this point I expect to have fPrimaries and fSecondaries energies in [GeV]
    
    
    // 1. update primary tracks in GeantTrack_v to their post-interaction state
    
    int numPrimaryTracks = fPrimaryTracks->numTracks; // number of primary tracks has been used
    //std::cout<<"WriteBackTracks: start, fPrimaryTracks: "<< numPrimaryTracks<<"\n";
    for (int ip = 0; ip < numPrimaryTracks; ++ip) {
        int indxP = fParentTrackIndices[ip]; // primary track indices in GeantTrack_v
        
        // set the selected process value for each Compton tracks in GeantTrack_v to -1
        // important to not process again the same track
        gTrackV.fProcessV[indxP] = -1;
        
        // check the post-interaction kinetic energy of the primary
        //double kinE = (fPrimaryTracks->E[ip]) - gTrackV.fMassV[indxP]; // should be [GeV]
        double kinE = fPrimaryTracks->E[ip]; //This is the kinEn
        
        if(kinE<0)
        {
           std::cout<<"fPrimaryTracks->E[ip]: "<<fPrimaryTracks->E[ip]<<" e gTrackV.fMassV[indxP]: "<<gTrackV.fMassV[indxP]<<"\n";
            std::cout<<"kinE for the primary is negative, ERROR. Exiting, bye bye.\n";
            exit(0);
        }
        
        
        if (kinE > fEnergyLimit) {
            // above tracking cut -> survived the phyisics interaction
            gTrackV.fEV[indxP]= fPrimaryTracks->E[ip]+gTrackV.fMassV[indxP];// update total E [GeV]
            //double totalP = std::sqrt(kinE * (kinE + 2.0 * gTrackV.fMassV[indxP])); // total P in [GeV]
            double totalP =sqrt(fPrimaryTracks->px[ip]*fPrimaryTracks->px[ip]+fPrimaryTracks->py[ip]*fPrimaryTracks->py[ip]+fPrimaryTracks->pz[ip]*fPrimaryTracks->pz[ip]);
            //std::cout<<"**** WriteBackTracks, primary mass: "<<gTrackV.fMassV[indxP]<<" and KinE: fPrimaryTracks->E["<<ip<<"]: "<<fPrimaryTracks->E[ip]<<" and totalP: "<<totalP<<"\n";
            
            double invTotalP = 1.0 / totalP;
            gTrackV.fPV[indxP] = totalP; // update total P
            
            //TRY ROTATION HERE
            //double oldxdir=gTrackV.fXdirV[indxP];
            //double oldydir=gTrackV.fYdirV[indxP];
            //double oldzdir=gTrackV.fZdirV[indxP];
    
            gTrackV.fXdirV[indxP] = fPrimaryTracks->px[ip] * invTotalP; // update x-dir
            gTrackV.fYdirV[indxP] = fPrimaryTracks->py[ip] * invTotalP; // update y-dir
            gTrackV.fZdirV[indxP] = fPrimaryTracks->pz[ip] * invTotalP; // update z-dir
            
            //RotateNewTrack(oldxdir,oldydir, oldzdir, gTrackV, indxP);

    
            
           
            
        } else {                                                      // apply tracking cut:
            gTrackV.fEdepV[indxP] += kinE;                              // put ekin to energy depo
            gTrackV.fStatusV[indxP] = kKilled;                          // kill the primary track
            gTrackV.fEV[indxP] = gTrackV.fMassV[indxP];                 //added
            // should make sure that at-rest process is invoked if needed (not no ofc.)
        }
    }
    //std::cout<<"WriteBackTracks: go for the secondaries\n";
    
    // 2. go for the secondaries
    
    int numSecondaries = fSecondaryTracks->numTracks;
    if (numSecondaries < 1) // if there are no secondaries return immediately
    {   //std::cout<<"WriteBackTracks: NO secondaries\n";
        return 0;}
    //std::cout<<"WriteBackTracks: some secondaries: "<<numSecondaries<<"\n";
    // get the GeantPropagator
    GeantPropagator *propagator = GeantPropagator::Instance();
    
    // A process with the same secondary can define its base properties here
    // like Compton: the only possible secondary is e-
    const int secPDG = 11;   // e- PDG code, valid only for Compton
    //  const int secGVcode = TPartIndex::I()->PartIndex(secPDG); // e- GV code
    // the rest mass and charge of the secondary particle in a general way
    const double secMass = vecphys::electron_mass_c2/1000;  //mb: necessary conversion from MeV to GeV
    const double secCharge = vecphys::electron_charge;
    
    int numInsertedTracks = 0;
    int numKilledTracks = 0;
    for (int isec = 0; isec < numSecondaries; ++isec) {
        
        //std::cout<<"primary track indices in GeantTrack_v\n";
        //std::cout<<"isec: "<<isec<<"\n";
        //std::cout<<" fSecondaryTracks->parentId[isec]: "<<fSecondaryTracks->parentId[isec]<<"\n";
        //std::cout<<" fParentTrackIndices[fSecondaryTracks->parentId[isec]]: "<<fParentTrackIndices[fSecondaryTracks->parentId[isec]]<<"\n";
        // primary track indices in GeantTrack_v
        //int indxP = fParentTrackIndices[fSecondaryTracks->parentId[isec]]; // ---> TO be VERIFIED!!!!
        int indxP = fParentTrackIndices[isec]; //vale solo per il Compton
        
        // check the kinetic energy of the secondary
        //double kinE = fSecondaryTracks->E[isec] - secMass; // should be [GeV]
        //TotE= KinE+RestMass
        double kinE = fSecondaryTracks->E[isec];
        
        if(kinE<0)
        {
            std::cout<<"fSecondaryTracks->E[isec]: "<<fSecondaryTracks->E[isec]<<" [GeV]\nSecMass: "<<secMass<<" [GeV]\nfEnergyLimit: "<<fEnergyLimit<<"\n";
            std::cout<<"kinE for the secondary is negative, ERROR. Exiting, bye bye\n";
            std::cout<<"fSecondaryTracks->E["<<isec<<"]: "<<fSecondaryTracks->E[isec]<<" [GeV] \nfEnergyLimit: "<<fEnergyLimit<<"\n";
            exit(0);
        }
        if (kinE > fEnergyLimit) {
            
            // above tracking cut -> insert into GeantTrack_v
            // get a temporary GeantTrack from the propagator
            GeantTrack &gTrack = propagator->GetTempTrack(tid);
            //std::cout<<"WriteBackTracks:0\n";
            // set some general properties: initial values or same as parent
            SetGeantTrack(gTrack, gTrackV, indxP); //Qui copio i dati del parent nella nuova traccia
        
            /*
             track.fEvent = tracks.fEventV[t];                                  YES
             track.fEvslot = tracks.fEvslotV[t];                                YES
             //          track.fParticle = nTotSecPart;          //index of this particle ???
             track.fPDG = secPDG;    // PDG code of this particle               NO ---> YES
             track.fGVcode = pid[i]; // GV index of this particle               NO ---> YES
             track.fEindex = 0;                                                 YES (CHANGED FROM -1 to 0!)
             track.fCharge = secPartPDG->Charge(); // charge of this particle   NO ---> YES
             #ifndef USE_VECGEOM_NAVIGATOR
             track.fCharge /=3.;                                                NO ---> YES
             #endif
             track.fProcess = 0;                                                YES
             track.fVindex = tracks.fVindexV[t];                                YES
             track.fNsteps = 0;                                                 YES
             //          track.fSpecies  = 0;
             track.fStatus = kNew;           // status of this particle                     YES
             track.fMass = secMass;          // mass of this particle                       NO ---> YES
             track.fXpos = tracks.fXposV[t]; // rx of this particle (same as parent)        YES
             track.fYpos = tracks.fYposV[t]; // ry of this particle (same as parent)        YES
             track.fZpos = tracks.fZposV[t]; // rz of this particle (same as parent)        YES
             track.fXdir = px / secPtot;     // dirx of this particle (before transform.)   NO ---> YES
             track.fYdir = py / secPtot;     // diry of this particle before transform.)    NO ---> YES
             track.fZdir = pz / secPtot;     // dirz of this particle before transform.)    NO ---> YES
             track.fP = secPtot;             // momentum of this particle                   NO ---> YES
             track.fE = secEtot;             // total E of this particle                    NO ---> YES
             track.fTime = tracks.fTimeV[t]; // global time                                 YES
             track.fEdep = 0.;                              YES
             track.fPstep = 0.;                             YES
             track.fStep = 0.;                              YES
             track.fSnext = 0.;                             YES
             track.fSafety = tracks.fSafetyV[t];            YES
             track.fBoundary = tracks.fBoundaryV[t];        YES
             track.fPending = kFALSE;                       YES
             *track.fPath = *tracks.fPathV[t];              YES
             *track.fNextpath = *tracks.fPathV[t];          YES*/
            
            
            // set additional members of gTrack
            double secTotalP = std::sqrt(kinE * (kinE + 2. * secMass)); // secondary total P [GeV]
            double invSecTotalP = 1.0 / secTotalP;
            gTrack.fPDG = secPDG;                                     // PDG code ->(e-)
            gTrack.fGVcode = TPartIndex::I()->PartIndex(secPDG);      // corresponding GV code
            gTrack.fCharge = secCharge;                               // charge
#ifndef USE_VECGEOM_NAVIGATOR
            gTrack.fCharge /=3.;
#endif
            
            gTrack.fMass = secMass;                                   // rest mass [GeV]
            //gTrack.fXdir = fSecondaryTracks->px[isec] * secPtot;    // direction (x,y,z)
            //gTrack.fYdir = fSecondaryTracks->py[isec] * secPtot;
            //gTrack.fZdir = fSecondaryTracks->pz[isec] * secPtot;
            gTrack.fXdir = fSecondaryTracks->px[isec] * invSecTotalP; //mb: direction (x,y,z)
            gTrack.fYdir = fSecondaryTracks->py[isec] * invSecTotalP;
            gTrack.fZdir = fSecondaryTracks->pz[isec] * invSecTotalP;
            
            /*gTrack.fXdir = gTrackV.fXdirV[indxP]; //TRY: direction (x,y,z)
            gTrack.fYdir = gTrackV.fYdirV[indxP];
            gTrack.fZdir = gTrackV.fZdirV[indxP];*/
            gTrack.fP = secTotalP;                              // total momentum
            gTrack.fE = fSecondaryTracks->E[isec]+secMass;      // total energy
            //TotE= KinE+RestMass

            // insert the new track
            propagator->AddTrack(gTrack);
            gTrackV.AddTrack(gTrack);
            // increase number of inserted tracks
            ++numInsertedTracks;
        } else
        {
            // apply tracking cut:
            // do not insert this track into GeantTrack_v
            gTrackV.fEdepV[indxP] += kinE; // put ekin to energy depo
            numKilledTracks++;
            
            // should make sure that at-rest process is invoked if needed (not no ofc.)
        }
    }
    //std::cout<<"CheckDirectionUnitVector dentro WriteTracks\n";
    CheckDirectionUnitVector(gTrackV, *fPrimaryTracks, *fSecondaryTracks);
    //std::cout<<"Put fSecondaryTracks->numTracks to zero.\n";
    fSecondaryTracks->numTracks=0;
    //std::cout<<"WriteBackTracks: numPrimaryTracks: "<<numPrimaryTracks<<", numKilledTracks: "<<numKilledTracks<<", inserite n. "<<numInsertedTracks<<" tracks. End\n";
    //std::cout<<"*** ApplyPostStepProcess. END ***\n";
    
   
    return numInsertedTracks;
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::SetGeantTrack(GeantTrack &left, GeantTrack_v &right, int ip) {
    
    /*
     track.fEvent = tracks.fEventV[t]; YES
     track.fEvslot = tracks.fEvslotV[t]; YES
     //          track.fParticle = nTotSecPart;          //index of this particle ???
     track.fPDG = secPDG;    // PDG code of this particle NO
     track.fGVcode = pid[i]; // GV index of this particle NO
     track.fEindex = 0; YES (CHANGED FROM -1 TO 0!)
     track.fCharge = secPartPDG->Charge(); // charge of this particle NO
     #ifndef USE_VECGEOM_NAVIGATOR
     track.fCharge /=3.; NO
     #endif
     track.fProcess = 0; Yes
     track.fVindex = tracks.fVindexV[t]; YES
     track.fNsteps = 0; YES
     //          track.fSpecies  = 0;
     track.fStatus = kNew;           // status of this particle YES
     track.fMass = secMass;          // mass of this particle NO
     track.fXpos = tracks.fXposV[t]; // rx of this particle (same as parent) YES
     track.fYpos = tracks.fYposV[t]; // ry of this particle (same as parent) YES
     track.fZpos = tracks.fZposV[t]; // rz of this particle (same as parent) YES
     track.fXdir = px / secPtot;     // dirx of this particle (before transform.)   NO
     track.fYdir = py / secPtot;     // diry of this particle before transform.)    NO
     track.fZdir = pz / secPtot;     // dirz of this particle before transform.)    NO
     track.fP = secPtot;             // momentum of this particle                   NO
     track.fE = secEtot;             // total E of this particle                    NO
     track.fTime = tracks.fTimeV[t]; // global time                                 YES
     track.fEdep = 0.;                              YES
     track.fPstep = 0.;                             YES
     track.fStep = 0.;                              YES
     track.fSnext = 0.;                             YES
     track.fSafety = tracks.fSafetyV[t];            YES
     track.fBoundary = tracks.fBoundaryV[t];        YES
     track.fPending = kFALSE;                       YES
     *track.fPath = *tracks.fPathV[t];              YES
     *track.fNextpath = *tracks.fPathV[t];          YES
     
     */
    
    left.fEvent = right.fEventV[ip];   // same as parent
    left.fEvslot = right.fEvslotV[ip]; // same as parent
    //    left.fPDG      = secPDG;                              // will be set
    //    left.fGVcode   = TPartIndex::I()->PartIndex(secPDG);  // will be set
    left.fEindex = 0;                 // init, mb: why??? track.fEindex = 0;     //Eindex= 0
    //    left.fCharge   = secCharge;                           // will be set
    //left.fProcess = -1;                // init: mb: track.fProcess = 0; here?
    left.fProcess = 0;                // Transport
    //  left.fIzero = 0;               // init
    left.fVindex = right.fVindexV[ip]; //mb: added
    left.fNsteps = 0;                  // init
    left.fStatus = kNew;               // new track
    
    //    left.fMass     = secMass;                             // will be set
    left.fXpos = right.fXposV[ip];     // position (same as parent)
    left.fYpos = right.fYposV[ip];
    left.fZpos = right.fZposV[ip];
    //    left.fXdir     = px*inv_Ptot;                         // dir. will be set
    //    left.fYdir     = py*inv_Ptot;
    //    left.fZdir     = pz*inv_Ptot;
    //    left.fP        = secPtot;                             // will be set
    //    left.fE        = secEtot;                             // will be set
    left.fTime = right.fTimeV[ip]; // mb: global time
    left.fEdep = 0.;                     // init
    left.fPstep = 0.;                    // init
    left.fStep = 0.;                     // init
    left.fSnext = 0.;                    // init
    left.fSafety = right.fSafetyV[ip];   // init to (same as parent)
    
    left.fBoundary = right.fBoundaryV[ip]; //added mb!
    //  left.fFrombdr = right.fFrombdrV[ip]; // init to (same as parent)
    left.fPending = kFALSE;              // init
    *left.fPath = *right.fPathV[ip];     // init
    *left.fNextpath = *right.fPathV[ip]; // init
}
