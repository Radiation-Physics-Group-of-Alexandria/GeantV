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
fPrimaryTracks(0), fSecondaryTracks(0) {}

//------------------------------------------------------------------------------
VecPhysOrchestrator::VecPhysOrchestrator(int processId, double energyLimit, int numtracks)
: fProcessId(processId), fEnergyLimit(energyLimit), fPrimaryTracks(0), fSecondaryTracks(0) {
    fPrimaryTracks = new GUTrack_v();
    fSecondaryTracks = new GUTrack_v();
    Allocator(numtracks);                         // allocate numtracks GUTrack_v for primaries
    // and secondaries and so on; will be increased if needed
    fVComptonModel = new vecphys::ComptonKleinNishina(0,-1); // create the vector physics model
    
    std::cout << "------------ processId is: = " << fProcessId << std::endl;
}

//------------------------------------------------------------------------------
VecPhysOrchestrator::~VecPhysOrchestrator() {
    Deallocator();
    delete fPrimaryTracks;
    delete fSecondaryTracks;
    delete fVComptonModel;
}

void VecPhysOrchestrator::DebugTracksEnergies(GeantTrack_v &gTrackV, int numtracks, GUTrack_v &primaries, bool checkIdentity)
{
    //std::cout<<"DebugTracksEnergies for "<<numtracks<<" tracks\n";
    int primariesCompton=primaries.numTracks;
    for (int i = 0; i < numtracks; ++i)
    {
        for (int j=0; j<primariesCompton; ++j)
            if(primaries.parentId[j]==i)
            {
                //if((gTrackV.fEV[i]-gTrackV.fMassV[i])!=primaries.E[j] && checkIdentity)
                //{
                    std::cout<<"DebugTracksEnergies: fPrimaries.E["<<j<<"]:  "<<primaries.E[j];
                    std::cout<<" gTrackV.fEV["<<i<<"]: "<<gTrackV.fEV[i]<<"\n";
                    std::cout<<" gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<"\n";
                    //std::cout<<"DegugTracksEnergies ERROR, exiting. bye bye.\n";
                    //return;
                 //   exit(0);
                //}
                
            }
    
    }
    //std::cout<<"\n\n*DegugTracksEnergies for "<<numtracks<<" tracks-END\n*\n\n";
}

//------------------------------------------------------------------------------
int VecPhysOrchestrator::ApplyPostStepProcess(GeantTrack_v &gTrackV, int numtracks, int tid) {
    
    FilterPrimaryTracks(gTrackV, numtracks);
    if (fPrimaryTracks->numTracks == 0) // if there is no track with Compton -> return
        return 0;
    DebugTracksEnergies(gTrackV, numtracks, *fPrimaryTracks, 1); //[GeV]
    ConvertEnergiesToVecPhys(); //[GeV]->[MeV]
    PerformInteraction(); //call KleinNishina Compton
    ConvertEnergiesFromVecPhys(); //[MeV]->[GeV]
    return WriteBackTracks(gTrackV, tid);
    
    
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::FilterPrimaryTracks(GeantTrack_v &gTrackV, int numtracks)
{   //This method should scan all the tracks and check the "compton" ones
    // count number of tracks in GeantTrack_v with selected process = Compton
    
    int numRealTracks = 0, totNumTracks = 0;
    for (int i = 0; i < numtracks; ++i)
    {
        ++totNumTracks;
        if (gTrackV.fProcessV[i] == fProcessId)
        {
            ++numRealTracks;
            //mb: if possible here set a flag on the Compton tracks - so that the tabulated physics does not process them
        }
        
    }
    
    if (numRealTracks < 1) {
        fPrimaryTracks->numTracks = 0;
        return;
    }
    
     // form the input GUTrack_v with the corresponding primary track members
     fPrimaryTracks->numTracks = numRealTracks;
     int j = 0;
     for (int i = 0; i < numtracks && j < numRealTracks; ++i) {
          std::cout<<"DEBUG: gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<" and  gTrackV.fProcessV["<<i<<"]: "<<gTrackV.fProcessV[i]<<"\n";
         if (gTrackV.fProcessV[i] == fProcessId) {
                fParentTrackIndices[j] = i;    // store index of this track in GeantTrack_v to store it back after the                          interaction
                double momentum = gTrackV.fPV[i];                       // total momentum
                fPrimaryTracks->px[j] = momentum * gTrackV.fXdirV[i];   // 3-momentum (px, py, pz)
                fPrimaryTracks->py[j] = momentum * gTrackV.fYdirV[i];
                fPrimaryTracks->pz[j] = momentum * gTrackV.fZdirV[i];
             
             
                //fPrimaryTracks->E[j] is the KinEn
                //gTrackV.fEV[i] is the Total Energy
                // TotalEnergy= KinEnergy + RestMass
                fPrimaryTracks->E[j] = gTrackV.fEV[i] - gTrackV.fMassV[i] ; // Kinetic energy
                std::cout<<"FilterPrimaryTracks, track: "<<gTrackV.fParticleV[i]<<" with fPrimaryTracks->E["<<j<<"]: "<<fPrimaryTracks->E[j]<< " and gTrackV.fEV["<<i<<"]: "<<gTrackV.fEV[i]<<" and gTrackV.fMassV["<<i<<"]: "<<gTrackV.fMassV[i]<<"\n";
     
                fTargetElements[j] = gTrackV.fEindexV[i];               // Z of the target atom
                //fPrimaryTracks->id[j] = i; //need to initialize the id on the PrimaryTracks because it will be used in the interact method when storing the secondaries
             
             //std::cout<<"Preparing the fPrimaryTracks\n";
             fPrimaryTracks->id[j] = j; //fPrimaryTracks id is j<=i
             fPrimaryTracks->parentId[j] = i; //fPrimaryTracks parentId is the same as fParentTrackIndices
                ++j;
            }
     }
    //std::cout<<"fPrimaryTracks Ready\n";

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
    gutrack_v.id            = new int[size];  // ??
    gutrack_v.parentId = new int[size]; // corresponding parent index in another GUTrack_v
    //   gutrack_v.proc          = new int[size];  // process index (we don't need this now) ???
    //   gutrack_v.x             = new double[size];   // (x,y,z) position (we don't need this now) ???
    //   gutrack_v.y             = new double[size];
    //   gutrack_v.z             = new double[size];
    gutrack_v.px = new double[size]; // momentum (px,py,pz)
    gutrack_v.py = new double[size];
    gutrack_v.pz = new double[size];
    gutrack_v.E = new double[size]; // total energy (kinetic E would be good as well)
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


void VecPhysOrchestrator::PerformInteraction() {
    fVComptonModel->Interact<vecCore::backend::VcVector>(*fPrimaryTracks, fTargetElements, *fSecondaryTracks);
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
    int index=fPrimaryTracks->numTracks;
    for (int i=0; i<index; ++i)
        fPrimaryTracks->E[i]*=1000;
}
void VecPhysOrchestrator::ConvertEnergiesFromVecPhys(){
    /*
     GeantV unit is [GeV] -> [1 MeV] = 10^-3
     VecPhys unit is [MeV] -> [1 MeV] = 1
     To go back from VecPhys to GeantV we need to convert the energies from MeV to GeV
     */
    int indexPrimary=fPrimaryTracks->numTracks;
    for (int i=0; i<indexPrimary; ++i)
        fPrimaryTracks->E[i]/=1000;
    int indexSecondary=fSecondaryTracks->numTracks;
    for (int i=0; i<indexSecondary; ++i)
        fSecondaryTracks->E[i]/=1000;
    
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
            std::cout<<"kinE for the primary is negative, ERROR. Exiting, bye bye\n";
            exit(0);
        }
        
        
        if (kinE > fEnergyLimit) {
            // above tracking cut -> survived the phyisics interaction
            gTrackV.fEV[indxP]= fPrimaryTracks->E[ip]+gTrackV.fMassV[indxP];// update total E [GeV]
            double totalP = std::sqrt(kinE * (kinE + 2.0 * gTrackV.fMassV[indxP])); // total P in [GeV]
            double invTotalP = 1.0 / totalP;
            gTrackV.fPV[indxP] = totalP; // update total P
            
            // assume that the physics model has already rotated the track (in a vectorized way);
            gTrackV.fXdirV[indxP] = fPrimaryTracks->px[ip] * invTotalP; // update x-dir
            gTrackV.fYdirV[indxP] = fPrimaryTracks->py[ip] * invTotalP; // update y-dir
            gTrackV.fZdirV[indxP] = fPrimaryTracks->pz[ip] * invTotalP; // update z-dir
            
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
        int indxP = fParentTrackIndices[fSecondaryTracks->parentId[isec]]; //occhio qui
        
        // check the kinetic energy of the secondary
        //double kinE = fSecondaryTracks->E[isec] - secMass; // should be [GeV]
        //TotE= KinE+RestMass
        double kinE = fSecondaryTracks->E[isec];
        
        if(kinE<0)
        {
            std::cout<<"fSecondaryTracks->E[isec]: "<<fSecondaryTracks->E[isec]<<" [GeV]\nSecMass: "<<secMass<<" [GeV]\nfEnergyLimit: "<<fEnergyLimit<<"\n";
            std::cout<<"kinE for the secondary is negative, ERROR. Exiting, bye bye\n";
            exit(0);
        }
        if (kinE > fEnergyLimit) {
            
            /*
   
             track.fPDG = secPDG;    // PDG code of this particle
             track.fGVcode = pid[i]; // GV index of this particle
             
             track.fCharge = secPartPDG->Charge(); // charge of this particle
             
             #ifndef USE_VECGEOM_NAVIGATOR --> PERCHé?
             track.fCharge /=3.;
             #endif
             
             track.fProcess = 0;
             track.fMass = secMass;          // mass of this particle
             track.fXdir = px / secPtot;     // dirx of this particle (before transform.)
             track.fYdir = py / secPtot;     // diry of this particle before transform.)
             track.fZdir = pz / secPtot;     // dirz of this particle before transform.)
             track.fP = secPtot;             // momentum of this particle
             track.fE = secEtot;             // total E of this particle
             track.fTime = tracks.fTimeV[t]; // global time
             track.fSafety = tracks.fSafetyV[t];
             track.fBoundary = tracks.fBoundaryV[t];
             
             propagator->AddTrack(track);
             tracks.AddTrack(track);
             
             ++nTotSecPart;
             }
             
             */
            
            // above tracking cut -> insert into GeantTrack_v
            // get a temporary GeantTrack from the propagator
            GeantTrack &gTrack = propagator->GetTempTrack(tid);
            //std::cout<<"WriteBackTracks:0\n";
            // set some general properties: initial values or same as parent
            SetGeantTrack(gTrack, gTrackV, indxP); //Qui copio i dati del parent nella nuova traccia
            
            // set additional members of gTrack
            double secTotalP = std::sqrt(kinE * (kinE + 2. * secMass)); // secondary total P [GeV]
            double invSecTotalP = 1.0 / secTotalP;
            gTrack.fPDG = secPDG;                                     // PDG code ->(e-)
            gTrack.fGVcode = TPartIndex::I()->PartIndex(secPDG);      // corresponding GV code
            gTrack.fCharge = secCharge;                               // charge
            gTrack.fMass = secMass;                                   // rest mass [GeV]
            //gTrack.fXdir = fSecondaryTracks->px[isec] * secPtot;    // direction (x,y,z)
            //gTrack.fYdir = fSecondaryTracks->py[isec] * secPtot;
            //gTrack.fZdir = fSecondaryTracks->pz[isec] * secPtot;
            gTrack.fXdir = fSecondaryTracks->px[isec] * invSecTotalP; //mb: direction (x,y,z)
            gTrack.fYdir = fSecondaryTracks->py[isec] * invSecTotalP;
            gTrack.fZdir = fSecondaryTracks->pz[isec] * invSecTotalP;
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
    //std::cout<<"Azzero fSecondaryTracks->numTracks\n";
    fSecondaryTracks->numTracks=0;
    //std::cout<<"WriteBackTracks: numPrimaryTracks: "<<numPrimaryTracks<<", numKilledTracks: "<<numKilledTracks<<", inserite n. "<<numInsertedTracks<<" tracks. End\n";
    return numInsertedTracks;
}

//------------------------------------------------------------------------------
void VecPhysOrchestrator::SetGeantTrack(GeantTrack &left, GeantTrack_v &right, int ip) {
    
    left.fEvent = right.fEventV[ip];   // same as parent
    left.fEvslot = right.fEvslotV[ip]; // same as parent
    //    left.fPDG      = secPDG;                              // will be set
    //    left.fGVcode   = TPartIndex::I()->PartIndex(secPDG);  // will be set
    left.fEindex = -1;                 // init, mb: why??? track.fEindex = 0;     //Eindex= 0
    //    left.fCharge   = secCharge;                           // will be set
 /*
#ifndef USE_VECGEOM_NAVIGATOR
    track.fCharge /=3.;
#endif*/
    
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
    
    left.fBoundary = right.fBoundaryV[ip]; //added mb
    //  left.fFrombdr = right.fFrombdrV[ip]; // init to (same as parent)
    left.fPending = kFALSE;              // init
    *left.fPath = *right.fPathV[ip];     // init
    *left.fNextpath = *right.fPathV[ip]; // init
}
