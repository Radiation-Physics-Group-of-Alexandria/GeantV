
#include "KleinNishinaComptonModel.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"

#include "Spline.h"
#include "GLIntegral.h"
#include "AliasTable.h"

#include "PhysicsParameters.h"

#include "Gamma.h"
#include "Electron.h"

#include "LightTrack.h"
#include "PhysicsData.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>


namespace geantphysics {
    
    
    
    KleinNishinaComptonModel::KleinNishinaComptonModel(const std::string &modelname)
    : EMModel(modelname) {
        
        fNumSamplingGammaEnergies    = 71;// 171=>25 per decada//71; // between the min/max gamma kinetic energies ? mb: why this value? - need some tuning
        fNumSamplingElectronEnergies = 54;              // at each energy grid points
        fMinGammaEnergy              =  1.0*geant::keV; // minimum kinetic energy of the interacting gamma
        fMaxGammaEnergy              = 10.0*geant::GeV; // maximum kinetic energy of the interacting gamma
        fGammaEnLMin                 = 0.0;
        fGammaEnILDelta              = 1.0;
        fSamplingGammaEnergies       = nullptr;
        fLSamplingGammaEnergies      = nullptr;
        
        fNumMaterialCuts             = 0;
        fNumDifferentMaterialECuts   = 0;
        fGlobalMatGCutIndxToLocal    = nullptr;
        fAliasData                   = nullptr; //alias data for each matrial-gammacut pairs
        fAliasSampler                = nullptr;
        fSecondaryInternalCode       = Electron::Definition()->GetInternalCode();
        
    }
    
    KleinNishinaComptonModel::~KleinNishinaComptonModel() {
        
        if (fSamplingGammaEnergies)
            delete [] fSamplingGammaEnergies;
        if (fLSamplingGammaEnergies)
            delete [] fLSamplingGammaEnergies;
        
        if (fGlobalMatGCutIndxToLocal)
            delete [] fGlobalMatGCutIndxToLocal;
        
        if (fAliasData)
            for (int i=0; i<fNumDifferentMaterialECuts; ++i)
                for (int j=0; j<fNumSamplingGammaEnergies; ++j) {
                    int indx = i*fNumSamplingGammaEnergies+j;
                    if (fAliasData[indx]) {
                        delete [] fAliasData[indx]->fXdata;
                        delete [] fAliasData[indx]->fYdata;
                        delete [] fAliasData[indx]->fAliasW;
                        delete [] fAliasData[indx]->fAliasIndx;
                        delete fAliasData[indx];
                    }
                }
        
        if (fAliasSampler)
            delete fAliasSampler;
        
    }
        
    void KleinNishinaComptonModel::Initialize() {
        EMModel::Initialize();
        InitSamplingTables();
        fAliasSampler          = new AliasTable();
        fSecondaryInternalCode = Electron::Definition()->GetInternalCode();
    }
    
    //
    double KleinNishinaComptonModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                                const Particle * ) {
        //N.B. TO BE PROPERLY IMPLEMENTED
        double xsec = 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        const Material *mat =  matcut->GetMaterial();
        const double *cuts  =  matcut->GetProductionCutsInEnergy();
        double gammacut     =  cuts[0];
        xsec = ComputeXSectionPerVolume(mat, gammacut, kinenergy);
        return xsec;
    }
    
    //
    double KleinNishinaComptonModel::ComputeXSectionPerVolume(const Material *mat, double prodcutenergy, double particleekin) {
        
        //N.B. TO BE PROPERLY IMPLEMENTED
        double xsec = 0.0;
        if (particleekin<=prodcutenergy) {
            return xsec;
        }
        
        // we will need the element composition of this material
        const Vector_t<Element*> theElements = mat->GetElementVector();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        int   numElems = theElements.size();
        double sum    = 0.0;
        
        for (int ielem=0; ielem<numElems; ++ielem) {
            //double zet  = theElements[ielem]->GetZ();
            ComputeXSectionPerAtom(theElements[ielem], particleekin);
        }
        return xsec;
        
    }
    
    static const double
    d1= 2.7965e-1*geant::barn, d2=-1.8300e-1*geant::barn,
    d3= 6.7527   *geant::barn, d4=-1.9798e+1*geant::barn,
    e1= 1.9756e-5*geant::barn, e2=-1.0205e-2*geant::barn,
    e3=-7.3913e-2*geant::barn, e4= 2.7079e-2*geant::barn,
    f1=-3.9178e-7*geant::barn, f2= 6.8241e-5*geant::barn,
    f3= 6.0480e-5*geant::barn, f4= 3.0274e-4*geant::barn;
    
    //
    double KleinNishinaComptonModel::ComputeXSectionPerAtom(const Element *elem, double gammaenergy){
        
        double xsec = 0.0;
        if (gammaenergy<GetLowEnergyUsageLimit() || gammaenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        
        double z=elem->GetZ();
        
        static const double a = 20.0 , b = 230.0 , c = 440.0;
        
        double p1Z = z*(d1 + e1*z + f1*z*z);
        double p2Z = z*(d2 + e2*z + f2*z*z);
        double p3Z = z*(d3 + e3*z + f3*z*z);
        double p4Z = z*(d4 + e4*z + f4*z*z);
        
        double T0  = 15.0 * geant::keV;
        if (z < 1.5) { T0 = 40.0 * geant::keV; }
        
        double X   = std::max(gammaenergy, T0) * geant::kInvElectronMassC2;
        xsec = p1Z*std::log(1.+2.*X)/X + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
        
        //  modification for low energy. (special case for Hydrogen)
        if (gammaenergy < T0) {
            static const double dT0 = geant::keV;
            X = (T0+dT0) / geant::kElectronMassC2 ;
            double sigma = p1Z*std::log(1.+2*X)/X
            + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
            double   c1 = -T0*(sigma-xsec)/(xsec*dT0);
            double   c2 = 0.150;
            if (z > 1.5) { c2 = 0.375-0.0556*std::log(z); }
            double    y = std::log(gammaenergy/T0);
            xsec *= std::exp(-y*(c1+c2*y));
        }
        return xsec;
    }
    
    
    int KleinNishinaComptonModel::SampleSecondaries(LightTrack &track, std::vector<LightTrack> & /*sectracks*/,
                                                    Geant::GeantTaskData *td) {
        int    numSecondaries      = 0;
        double gammaekin0           = track.GetKinE();
        const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
        const double *cuts         = matCut->GetProductionCutsInEnergy();
        double electroncut         = cuts[0];
        
        // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is above fHighEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is below gamma production cut and do nothing if yes;
        if (gammaekin0<GetLowEnergyUsageLimit() || gammaekin0>GetHighEnergyUsageLimit() || gammaekin0<=electroncut) {
            return numSecondaries;
        }
        // sample gamma energy
        // here we need 3 random number + 1 later for photon direction phi sampling
        double *rndArray = td->fDblArray;
        td->fRndm->uniform_array(4, rndArray);
        
        
        //
        double cosTheta = 0.0;
        double sinTheta = 0.0;
        double phi      = 0.0;
        
        // store original gamma directions in the lab frame
        double gamDirX0=track.GetDirX();
        double gamDirY0=track.GetDirY();
        double gamDirZ0=track.GetDirZ();
        
        
        // sample gamma energy and gamma scattering angle in the scattering frame i.e. which z-dir points to the original gamma direction
        double gammaekin1  = SamplePhotonEnergy(matCut, gammaekin0, rndArray[0], rndArray[1], rndArray[2]);
        SamplePhotonDirection(gammaekin0, gammaekin1, sinTheta, cosTheta, phi, rndArray[3]);
        
        // Uncomment to sample with the composition rejection
        //double gammaekin1  = SamplePhotonEnergyAndDirection(matCut, gammaekin0, sinTheta, cosTheta,td);
        //phi      = geant::kTwoPi*(rndArray[3]);
        
        // new gamma direction in the scattering frame
        double gamDirX1  = sinTheta*std::cos(phi);
        double gamDirY1  = sinTheta*std::sin(phi);
        double gamDirZ1  = cosTheta;
        
        
        // rotate new gamma direction to the lab frame:
        RotateToLabFrame(gamDirX1, gamDirY1, gamDirZ1, gamDirX0, gamDirY0, gamDirZ0);
        
        double edep=0.;
        
        if(gammaekin1 > electroncut) {
            
            //Update primary track direction and energy
            track.SetDirX(gamDirX1);
            track.SetDirY(gamDirY1);
            track.SetDirZ(gamDirZ1);
            // update primary track kinetic energy
            track.SetKinE(gammaekin1);
            
        }
        else
        {
            track.SetTrackStatus(LTrackStatus::kKill);
            track.SetKinE(0.0);
            edep = gammaekin1;
        }
        
        double eDirX ;
        double eDirY ;
        double eDirZ ;
        double eekin = gammaekin0 - gammaekin1;
        
        
        //if the energy of the e- is greater than the cut, the secondary particle is effectively created
        if(eekin > electroncut) {
            
            eDirX = gammaekin0*gamDirX0 - gammaekin1*gamDirX1;
            eDirY = gammaekin0*gamDirY0 - gammaekin1*gamDirY1;
            eDirZ = gammaekin0*gamDirZ0 - gammaekin1*gamDirZ1;
            double norm  = 1.0/std::sqrt(eDirX*eDirX + eDirY*eDirY + eDirZ*eDirZ);
            
            // create the secondary partcile i.e. the electron
            numSecondaries = 1;
            // NO it can be dropped if we make sure that these secondary vectors are at least size 2
            //  PhysicsData *physData = td->fPhysicsData;
            // current capacity of secondary track container
            int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
            // currently used secondary tracks in the container
            int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();
            if (curSizeOfSecList-curNumUsedSecs<numSecondaries) {
                td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
            }
            int secIndx = curNumUsedSecs;
            curNumUsedSecs +=numSecondaries;
            td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
            
            std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
            
            
            // this is known since it is a secondary track
            //  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
            sectracks[secIndx].SetDirX(eDirX*norm); //mb: controllare questa operazione
            sectracks[secIndx].SetDirY(eDirY*norm);
            sectracks[secIndx].SetDirZ(eDirZ*norm);
            
            
            sectracks[secIndx].SetKinE(eekin);
            sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // electron GV code
            sectracks[secIndx].SetMass(geant::kElectronMassC2);
            sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
            
            // these are known from the parent GeantTrack
            //  sectracks[secIndx].SetMaterialCutCoupleIndex(track.GetMaterialCutCoupleIndex());
            //  sectracks[secIndx].SetNumOfInteractionLegthLeft(-1.); // i.e. need to sample in the step limit
            //  sectracks[secIndx].SetInvTotalMFP(0.);
            //  sectracks[secIndx].SetStepLength(0.);
            //  sectracks[secIndx].SetEnergyDeposit(0.);
            //  sectracks[secIndx].SetTime(??);
            //  sectracks[secIndx].SetWeight(??);
            //  sectracks[secIndx].SetProcessIndex(-1); // ???
            //  sectracks[secIndx].SetTargetZ(-1);
            //  sectracks[secIndx].SetTargetN(-1);
            //
            
            /* Forse bisogna usare questa per garantire la conservazione del momento e dell'energia
             // compute the primary e-/e+ post interaction direction: from momentum vector conservation
             //  double elInitTotalEnergy   = ekin+geant::kElectronMassC2;  // initial total energy of the e-/e+
             double elInitTotalMomentum = std::sqrt(ekin*(ekin+2.0*geant::kElectronMassC2));
             // final momentum of the e-/e+ in the lab frame
             double elDirX = elInitTotalMomentum*track.GetDirX() - gammaEnergy*gamDirX;
             double elDirY = elInitTotalMomentum*track.GetDirY() - gammaEnergy*gamDirY;
             double elDirZ = elInitTotalMomentum*track.GetDirZ() - gammaEnergy*gamDirZ;
             // normalisation
             double norm  = 1.0/std::sqrt(elDirX*elDirX + elDirY*elDirY + elDirZ*elDirZ);
             // update primary track direction
             track.SetDirX(elDirX*norm);
             track.SetDirY(elDirY*norm);
             track.SetDirZ(elDirZ*norm);
             // update primary track kinetic energy
             track.SetKinE(ekin-gammaEnergy);
             */
            
        } else {
            edep += eekin;
        }
        
        if(edep > 0.0) {
            track.SetEnergyDeposit(edep);
        }
        // return with number of secondaries i.e. 1 electron
        return numSecondaries;
    }
    
    double KleinNishinaComptonModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle*) const {
        double mine = (matcut->GetProductionCutsInEnergy())[0]; // gamma production cut in the given material-cuts
        return mine;
    }
    
    
    /**
     
     *
     *
     *
     */
    
    double KleinNishinaComptonModel::SamplePhotonEnergy(const MaterialCuts *matcut, double gammaekin, double r1, double r2, double r3){
        
        int mcindx       = matcut->GetIndex(); //GetIndex of the material cut
        int macindxlocal = fGlobalMatGCutIndxToLocal[mcindx]; //From global index to local index
        // the location of the first-gamma-energy lin-alias data of this mat-cut
        int indxstart    = macindxlocal*fNumSamplingGammaEnergies;//
        // determine gamma energy lower grid point
        double lgenergy  = std::log(gammaekin);//lower gamma energy grid point
        int genergyindx  = (int) ((lgenergy-fGammaEnLMin)*fGammaEnILDelta);
        
        if (genergyindx>=fNumSamplingGammaEnergies-1)
            genergyindx = fNumSamplingGammaEnergies-2;
        
        double plowgener = (fLSamplingGammaEnergies[genergyindx+1]-lgenergy)*fGammaEnILDelta;
        indxstart +=genergyindx;
        if (r1>plowgener|| !(fAliasData[indxstart]))
            ++indxstart;
        
        // sample the outgoing photon energy
        double gammaeout = fAliasSampler->SampleLinear(fAliasData[indxstart]->fXdata, fAliasData[indxstart]->fYdata,
                                                       fAliasData[indxstart]->fAliasW, fAliasData[indxstart]->fAliasIndx,
                                                       fNumSamplingElectronEnergies,r2,r3);
        return gammaeout;
        
    }
    
    
    void KleinNishinaComptonModel::SamplePhotonDirection(double gammaEnIn, double gammaEnOut, double &sinTheta, double &cosTheta, double &phi, double rnd){
        
        double E0_m = gammaEnIn * geant::kInvElectronMassC2 ;
        double epsilon = gammaEnOut/gammaEnIn;
        double onecost = (1.- epsilon)/(epsilon*E0_m);
        
        cosTheta= 1. - onecost;
        // check cosTheta limit
        if (cosTheta>1.0) {
            cosTheta = 1.0;
        }
        sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));
        phi      = geant::kTwoPi*(rnd);
        
    }
    
    
    
    //Added from geant4
    //static const double
    //d1= 2.7965e-1*geant::barn, d2=-1.8300e-1*geant::barn,
    //d3= 6.7527   *geant::barn, d4=-1.9798e+1*geant::barn,
    //e1= 1.9756e-5*geant::barn, e2=-1.0205e-2*geant::barn,
    //e3=-7.3913e-2*geant::barn, e4= 2.7079e-2*geant::barn,
    //f1=-3.9178e-7*geant::barn, f2= 6.8241e-5*geant::barn,
    //f3= 6.0480e-5*geant::barn, f4= 3.0274e-4*geant::barn;
    static const int nlooplim = 1000;
    
    // Sample photon energy and direction with the composition rejection
    double KleinNishinaComptonModel::SamplePhotonEnergyAndDirection(double gammaekin, double &sinTheta, double &cosTheta, Geant::GeantTaskData *td){
        
        double *rndArray = td->fDblArray;
        double E0_m = gammaekin / geant::kElectronMassC2 ;
        
        // sample the energy rate of the scattered gamma
        double epsilon, epsilonsq, greject, sint2, onecost ;
        double eps0       = 1./(1. + 2.*E0_m);
        double epsilon0sq = eps0*eps0;
        double alpha1     = - std::log(eps0);
        double alpha2     = alpha1 + 0.5*(1.- epsilon0sq);
        int nloop = 0;
        do {
            ++nloop;
            // false interaction if too many iterations
            if(nloop > nlooplim) {
                return -1; }
            td->fRndm->uniform_array(3, rndArray);
            
            if ( alpha1 > alpha2*rndArray[0] ) {
                epsilon   = std::exp(-alpha1*rndArray[1]);
                epsilonsq = epsilon*epsilon;
                
            } else {
                epsilonsq = epsilon0sq + (1.- epsilon0sq)*rndArray[1];
                epsilon   = sqrt(epsilonsq);
            };
            
            onecost = (1.- epsilon)/(epsilon*E0_m);
            sint2   = onecost*(2.-onecost);
            greject = 1. - epsilon*sint2/(1.+ epsilonsq);
            
        } while (greject < rndArray[2]);
        
        // scattered gamma angles. ( Z - axis along the parent gamma)
        if(sint2 < 0.0) { sint2 = 0.0; }
        cosTheta = 1. - onecost;
        sinTheta = sqrt (sint2);
        double gamEnergy1=epsilon*gammaekin;
        return gamEnergy1;
    }
    
    
    void KleinNishinaComptonModel::InitSamplingTables() {
        // set up the common gamma energy grid
        if (fSamplingGammaEnergies) {
            delete [] fSamplingGammaEnergies;
            fSamplingGammaEnergies = nullptr;
        }
        fSamplingGammaEnergies  = new double[fNumSamplingGammaEnergies];
        fLSamplingGammaEnergies = new double[fNumSamplingGammaEnergies];
        
        fGammaEnLMin    = std::log(fMinGammaEnergy);
        double delta = std::log(fMaxGammaEnergy/fMinGammaEnergy)/(fNumSamplingGammaEnergies-1.0);
        fGammaEnILDelta = 1.0/delta;
        fSamplingGammaEnergies[0]  = fMinGammaEnergy;
        fLSamplingGammaEnergies[0] = fGammaEnLMin;
        fSamplingGammaEnergies[fNumSamplingGammaEnergies-1]  = fMaxGammaEnergy;
        fLSamplingGammaEnergies[fNumSamplingGammaEnergies-1] = std::log(fMaxGammaEnergy);
        
        /*
         std::cout<<"\n\n*****\nfSamplingGammaEnergies[0]=fMinGammaEnergy= "<< fMinGammaEnergy<<std::endl;
         std::cout<<"fMaxGammaEnergy: "<< fMaxGammaEnergy<<std::endl;
         std::cout<<"fNumSamplingGammaEnergies: "<< fNumSamplingGammaEnergies<<std::endl;
         std::cout<<"delta energy: "<< delta<<std::endl;
         std::cout<<"fGammaEnILDelta: "<< 1.0/delta<<std::endl;
         std::cout<<"fGammaEnLMin: "<< std::log(fMinGammaEnergy)<<std::endl;
         std::cout<<"fLSamplingGammaEnergies[0] = fGammaEnLMin= "<< fGammaEnLMin<<std::endl;
         std::cout<<"fSamplingGammaEnergies["<<fNumSamplingGammaEnergies-1<<"]: "<< std::log(fMaxGammaEnergy)<<std::endl;
         */
        for (int i=1; i<fNumSamplingGammaEnergies-1; ++i) {
            fLSamplingGammaEnergies[i] = fGammaEnLMin+i*delta;
            fSamplingGammaEnergies[i]  = std::exp(fGammaEnLMin+i*delta);
            //std::cerr<<" E("<<i<<") = "<<fSamplingGammaEnergies[i]/geant::GeV<<std::endl;
        }  // fMinElecEnergy fMaxElecEnergy az fLoadDCSElectronEnergyGrid[0] es fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1]
        //std::cerr<<" E("<<fNumSamplingGammaEnergies-1<<") = "<<fSamplingGammaEnergies[fNumSamplingGammaEnergies-1]/geant::GeV<<std::endl;
        
        // - get number of different material-gammacut pairs
        // - allocate space and fill the global to local material-cut index map
        const std::vector<MaterialCuts*> theMaterialCutsTable = MaterialCuts::GetTheMaterialCutsTable();
        int fNumMaterialCuts = theMaterialCutsTable.size();
        if (fGlobalMatGCutIndxToLocal) {
            delete [] fGlobalMatGCutIndxToLocal;
            fGlobalMatGCutIndxToLocal = nullptr;
        }
        fGlobalMatGCutIndxToLocal = new int[fNumMaterialCuts];
        //std::cerr<<" === Number of global Material-Cuts = "<<fNumMaterialCuts<<std::endl;
        
        // count diffenet material-gammacut pairs and set to global to local mat-cut index map
        int oldnumDif = fNumDifferentMaterialECuts;
        int oldnumSEE = fNumSamplingGammaEnergies;
        fNumDifferentMaterialECuts = 0;
        for (int i=0; i<fNumMaterialCuts; ++i) {
            // if the current MaterialCuts does not belong to the current active regions
            if (!IsActiveRegion(theMaterialCutsTable[i]->GetRegionIndex())) {
                fGlobalMatGCutIndxToLocal[i] = -1;
                continue;
            }
            bool isnew = true;
            int j = 0;
            for (; j<fNumDifferentMaterialECuts; ++j) {
                if (theMaterialCutsTable[i]->GetMaterial()->GetIndex()==theMaterialCutsTable[j]->GetMaterial()->GetIndex() &&
                    theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]==theMaterialCutsTable[j]->GetProductionCutsInEnergy()[0]) {
                    isnew = false;
                    break;
                }
            }
            if (isnew) {
                fGlobalMatGCutIndxToLocal[i] = fNumDifferentMaterialECuts;
                ++fNumDifferentMaterialECuts;
            } else {
                fGlobalMatGCutIndxToLocal[i] = j;
            }
        }
        std::cerr<<" === Number of local Material-Cuts = "<<fNumDifferentMaterialECuts<<std::endl;
        
        // allocate space for the material-gcut sampling tables and init these pointers to null
        if (fAliasData)
            for (int i=0; i<oldnumDif; ++i)
                for (int j=0; j<oldnumSEE; ++j) {
                    int indx = i*oldnumSEE+j;
                    if (fAliasData[indx]) {
                        delete [] fAliasData[indx]->fXdata;
                        delete [] fAliasData[indx]->fYdata;
                        delete [] fAliasData[indx]->fAliasW;
                        delete [] fAliasData[indx]->fAliasIndx;
                        fAliasData[indx] = nullptr;
                    }
                }
        
        int *isdone = new int[fNumDifferentMaterialECuts]();
        int  idum   = fNumDifferentMaterialECuts*fNumSamplingGammaEnergies;
        fAliasData = new LinAlias*[idum];
        for (int i=0; i<idum; ++i)
            fAliasData[i] = nullptr;
        
        for (int i=0; i<fNumMaterialCuts; ++i) {
            //std::cerr<<"   See if Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
            int localindx = fGlobalMatGCutIndxToLocal[i];
            if (localindx<0) {
                continue;
            }
            int ialias    = localindx*fNumSamplingGammaEnergies;
            if (!isdone[localindx]) { // init for this material-gammacut pair if it has not been done yet
                //       std::cerr<<"   -> Will init for Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
                BuildOneLinAlias(ialias, theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]);
                isdone[localindx] = 1;
            }
        }
        delete [] isdone;
        // test
        for (int i=0; i<fNumMaterialCuts; ++i)
            std::cerr<<"     --> Global MatCut-indx = "<< i << " local indx = "<<fGlobalMatGCutIndxToLocal[i] <<std::endl;
    }
    
    
    
    double KleinNishinaComptonModel::CalculateDiffCrossSection(double energy0, double energy1)
    {
        double E0_m = energy0 / geant::kElectronMassC2;
        double epsilon = energy1 / energy0;
        
        double onecost = (1. - epsilon) / (epsilon * E0_m);
        double sint2 = onecost * (2. - onecost);
        double greject = 1. - epsilon * sint2 / (1. + epsilon * epsilon);
        double dsigma = (epsilon + 1. / epsilon) * greject;
        
        return dsigma;
    }
    
    // for the common particle kinetic energy grid points
    //this loop is iterating through all the GammaIn energies and building the corresponding AliasTables - given a certain GammaE_In0 is building the Alias table calculating the p.d.f. - differential cross section
    // igener                           = indexGammaEnergy
    // gener                            = gammaEnergy
    // fNumSamplingGammaEnergies        = total number of gamma sampling energies
    // fNumSamplingElectronEnergies     = total number of electron sampling energies (gamma out sampling energies)
    // fSamplingGammaEnergies[]         = vector of dimension fNumSamplingGammaEnergies storing all the gamma sampling energies
    //
    void KleinNishinaComptonModel::BuildOneLinAlias(int ialias, double ecut){
        
        for (int igener=0; igener<fNumSamplingGammaEnergies; ++igener) {
            
            double gEnergy = fSamplingGammaEnergies[igener];
            
            if (gEnergy>ecut) { // otherwise no electron production at that gamma energy so let the sampling table to be nullptr
                
                // create the alias data struct
                fAliasData[ialias] = new LinAlias();
                fAliasData[ialias]->fXdata     = new double[fNumSamplingElectronEnergies]();
                fAliasData[ialias]->fYdata     = new double[fNumSamplingElectronEnergies]();
                fAliasData[ialias]->fAliasW    = new double[fNumSamplingElectronEnergies]();
                fAliasData[ialias]->fAliasIndx = new int[fNumSamplingElectronEnergies]();
                
                // Algorithm:
                // Constant: 1. number of Tables (Gamma input energies = fNumSamplingGammaEnergies
                // Constant: 2. number of Bins (how many bins for every table - adaptive binning)
                // Procedure:
                // a) For every table (for every igener) calculate the range of outgoing energies [Eg_out_min, Eg_out_max]
                // b) Calculate the middlepoint and calculate the differential cross section in the three points (initial, medium and final)
                // c) Calculate the linear interpolation between those points li(xi) li(xj)
                // d) Calculate the differential crossSection in the middle points xsec(xi) xsec(xj)
                // f) Depending on where is the greater error between |xsec(xi)-li(xi)| and |xsec(xj)-li(xj)| put the next bin
                // g) Iterate the procedure until the maximum number of bins is reached
                // h) At the end a pdf is produced with the corresponding energy points (x - coordinates)
                
                
                double kappa= gEnergy/geant::kElectronMassC2;
                
                double gEnergyFractionMin=(1/(2*kappa+1));
                //double gEnergyFractionMax=(1.);
                
                double gEnergyOutMin=gEnergyFractionMin*gEnergy;
                double gEnergyOutMax=gEnergy;
                
                double middlepoint=0.5*(gEnergyOutMax+gEnergyOutMin);
                
                // fill the first, last and middle value (NB: we should convert to 0-1)
                fAliasData[ialias]->fXdata[0] =  gEnergyOutMin;
                fAliasData[ialias]->fYdata[0] =  CalculateDiffCrossSection(gEnergy, gEnergyOutMin);
                
                fAliasData[ialias]->fXdata[1] =  middlepoint;
                fAliasData[ialias]->fYdata[1] =  CalculateDiffCrossSection(gEnergy, middlepoint);
                
                fAliasData[ialias]->fXdata[2] =  gEnergyOutMax;
                fAliasData[ialias]->fYdata[2] =  CalculateDiffCrossSection(gEnergy, gEnergyOutMax);
                
                int numdata=3;
                /*
                 std::cout<<"fAliasData["<<ialias<<"]->fXdata[0]: "<<fAliasData[ialias]->fXdata[0]<<std::endl;
                 std::cout<<"fAliasData["<<ialias<<"]->fYdata[0]: "<<fAliasData[ialias]->fYdata[0]<<std::endl;
                 std::cout<<"fAliasData["<<ialias<<"]->fXdata[1]: "<<fAliasData[ialias]->fXdata[1]<<std::endl;
                 std::cout<<"fAliasData["<<ialias<<"]->fYdata[1]: "<<fAliasData[ialias]->fYdata[1]<<std::endl;
                 std::cout<<"fAliasData["<<ialias<<"]->fXdata[2]: "<<fAliasData[ialias]->fXdata[2]<<std::endl;
                 std::cout<<"fAliasData["<<ialias<<"]->fYdata[2]: "<<fAliasData[ialias]->fYdata[2]<<std::endl;
                 */
                
                // expand the data up to maximum ---> dynamic binning based on the linear interpolator + differential cross-section
                while(numdata<fNumSamplingElectronEnergies)
                {
                    // find the lower index of the bin, where we have the biggest linear interp. error compared to Diff cross section value
                    double maxerr     = 0.0; // value of the current maximum error
                    double thexval    = 0.0;
                    double theyval    = 0.0;
                    int    maxerrindx = 0;   // the lower index of the corresponding bin
                    
                    for (int i=0; i<numdata-1; ++i)
                    {
                        double xx = 0.5*(fAliasData[ialias]->fXdata[i]+fAliasData[ialias]->fXdata[i+1]); // mid point
                        double yy = 0.5*(fAliasData[ialias]->fYdata[i]+fAliasData[ialias]->fYdata[i+1]); // lin func val at the mid point
                        
                        double xsecVal=CalculateDiffCrossSection(gEnergy, xx);
                        double err   = std::fabs(yy-xsecVal);
                        if (err>maxerr) {
                            maxerr     = err;
                            maxerrindx = i;
                            thexval    = xx;
                            theyval    = xsecVal;
                        }
                        
                    }
                    // extend x,y data by putting the calculated differential cross section value at the mid point of the highest error bin
                    // first shift all values to the right
                    for (int j=numdata; j>maxerrindx+1; --j) {
                        fAliasData[ialias]->fXdata[j] = fAliasData[ialias]->fXdata[j-1];
                        fAliasData[ialias]->fYdata[j] = fAliasData[ialias]->fYdata[j-1];
                    }
                    // fill x mid point
                    fAliasData[ialias]->fXdata[maxerrindx+1] = thexval;
                    fAliasData[ialias]->fYdata[maxerrindx+1] = theyval;
                    // increase number of data
                    ++numdata;
                }//end of dynamic binning
                
                // set up a linear alias sampler on this data
                AliasTable *alst = new AliasTable();
                alst->PreparLinearTable(fAliasData[ialias]->fXdata, fAliasData[ialias]->fYdata,
                                        fAliasData[ialias]->fAliasW, fAliasData[ialias]->fAliasIndx,
                                        fNumSamplingElectronEnergies);
                
                delete alst;
                
            }
            ++ialias;
        } //end of for
        
    }
    
    /*  void KleinNishinaComptonModel::LoadDCSData()
     {
     using geant::MeV;
     using geant::millibarn;
     
     // get the path to the main physics data directory
     char *path = std::getenv("GEANT_PHYSICS_DATA");
     if (!path) {
     std::cerr<<"******   ERROR in KleinNishinaComptonModel::LoadDCSData() \n"
     <<"         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
     <<"         environmental variable to the location of Geant data directory!\n"
     <<std::endl;
     exit(1);
     }
     
     char baseFilename[512];
     switch (fDataFileIndx) {
     case 0: sprintf(baseFilename,"%s/brems/NIST_BREM1/nist_brems_",path);
     break;
     case 1: sprintf(baseFilename,"%s/brems/NIST_BREM/nist_brems_",path);
     break;
     case 2: sprintf(baseFilename,"%s/brems/NRC_BREM/nist_brems_",path);
     break;
     default: sprintf(baseFilename,"%s/brems/NIST_BREM1/nist_brems_",path);
     }
     
     FILE *f = nullptr;
     char filename[512];
     sprintf(filename,"%sgrid",baseFilename);
     f = fopen(filename,"r");
     if (!f) {
     std::cerr<<"******   ERROR in KleinNishinaComptonModel::LoadDCSData() \n"
     <<"         "<< filename << " could not be found!\n"
     <<std::endl;
     exit(1);
     }
     // before we take the fDCSMaxZet make sure that the fDCSForElements if free
     if (fLoadDCSForElements) {
     for (int i=0; i<fDCSMaxZet; ++i)
     if (fLoadDCSForElements[i]) {
     delete [] fLoadDCSForElements[i];
     fLoadDCSForElements[i] = nullptr;
     }
     delete [] fLoadDCSForElements;
     }
     
     fscanf(f,"%d%d%d",&fDCSMaxZet,&fLoadDCSNumElectronEnergies,&fLoadDCSNumReducedPhotonEnergies);
     // allocate space for the elemental DCS data, for the electron energy and reduced photon energy grids and load them
     fLoadDCSForElements = new double*[fDCSMaxZet];
     for (int i=0; i<fDCSMaxZet; ++i)
     fLoadDCSForElements[i] = nullptr;
     if (fLoadDCSElectronEnergyGrid) {
     delete [] fLoadDCSElectronEnergyGrid;
     fLoadDCSElectronEnergyGrid = nullptr;
     }
     fLoadDCSElectronEnergyGrid = new double[fLoadDCSNumElectronEnergies];
     if (fLoadDCSReducedPhotonEnergyGrid) {
     delete [] fLoadDCSReducedPhotonEnergyGrid;
     fLoadDCSReducedPhotonEnergyGrid = nullptr;
     }
     fLoadDCSReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
     for (int i=0;i<fLoadDCSNumElectronEnergies;++i) {
     fscanf(f,"%lf",&(fLoadDCSElectronEnergyGrid[i]));
     fLoadDCSElectronEnergyGrid[i] *= MeV;   // change to internal energy units
     }
     for (int i=0;i<fLoadDCSNumReducedPhotonEnergies;++i)
     fscanf(f,"%lf",&(fLoadDCSReducedPhotonEnergyGrid[i]));
     fclose(f);
     //
     // for (int i=0;i<fDCSNumElectronEnergies;++i)
     //   printf("%.6e\n",fDCSElectronEnergyGrid[i]);
     // for (int i=0;i<fDCSNumReducedPhotonEnergies;++i)
     //   printf("%.6e\n",fDCSReducedPhotonEnergyGrid[i]);
     
     // now go for each element that we have in the global element table and load data for them
     int numDCSdataPerElement = fLoadDCSNumElectronEnergies*fLoadDCSNumReducedPhotonEnergies;
     const Vector_t<Element*> theElements= Element::GetTheElementTable();
     
     int numElements = theElements.size();
     for (int i=0; i<numElements; ++i) {
     int zet = std::lrint(theElements[i]->GetZ());
     sprintf(filename,"%s%d",baseFilename,zet);
     f = fopen(filename,"r");
     if (!f) {
     std::cerr<<"******   ERROR in KleinNishinaComptonModel::LoadDCSData() \n"
     <<"         "<< filename << " could not be found!\n"
     <<std::endl;
     exit(1);
     }
     // allocate space for this elemental DCS
     fLoadDCSForElements[zet-1] = new double[numDCSdataPerElement];
     for (int j=0; j<numDCSdataPerElement; ++j) {
     fscanf(f,"%lf",&(fLoadDCSForElements[zet-1][j]));
     fLoadDCSForElements[zet-1][j] *= millibarn;   // change to internal units
     }
     fclose(f);
     }
     }*/
    
}   // namespace geantphysics
