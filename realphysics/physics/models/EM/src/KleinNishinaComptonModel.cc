
#include "KleinNishinaComptonModel.h"

//#include "PhysicalConstants.h"

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
    
    void KleinNishinaComptonModel::Initialize() {
        EMModel::Initialize();
        Initialise();
        // if it needs element selector: particle is coded in the model in case of this model and due to the cut dependence
        // we we need element selectors per MaterialCuts
        //  InitialiseElementSelectors(this, nullptr, false);
    }
    
    /*
    double KleinNishinaComptonModel::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle * ,
                                                 bool istotal){
        const Material *mat =  matcut->GetMaterial();
        const double *cuts  =  matcut->GetProductionCutsInEnergy();
        double gammacut     =  cuts[0];
        if (istotal) {
            // for the total stopping power we just need a gamma production cut >=kinenergy
            gammacut = 1.01*kinenergy;
        }
        return ComputeDEDXPerVolume(mat, gammacut, kinenergy);
    }
    */
    
    
    /*
    double KleinNishinaComptonModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                                const Particle * /*particle*/ /*) {
     
                                                                    double xsec = 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        const Material *mat =  matcut->GetMaterial();
        const double *cuts  =  matcut->GetProductionCutsInEnergy();
        double gammacut     =  cuts[0];
        xsec = ComputeXSectionPerVolume(mat, gammacut, kinenergy);
        return xsec;
    }*/
    
    //mb: ok modified
    /*double KleinNishinaComptonModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle*) {
        
        double xsec = 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        const Material *mat =  matcut->GetMaterial();
        const double *cuts  =  matcut->GetProductionCutsInEnergy();
        double electroncut     =  cuts[0];
        xsec = ComputeXSectionPerAtom(elem, mat, electroncut, kinenergy);
        return xsec;
    }*/
    
    
    int KleinNishinaComptonModel::SampleSecondaries(LightTrack &track, std::vector<LightTrack> & /*sectracks*/,
                                                    Geant::GeantTaskData *td) {
        //std::cout<<"SampleSecondaries"<<std::endl;
        int    numSecondaries      = 0;
        double gammaekin0           = track.GetKinE();
        //std::cout<<"gamEnergy0: "<<gammaekin0<<std::endl;
        const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
        const double *cuts         = matCut->GetProductionCutsInEnergy();
        double electroncut         = cuts[0];
        //std::cout<<"electroncut: "<<electroncut<<std::endl;
        // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is above fHighEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is below gamma production cut and do nothing if yes;
        if (gammaekin0<GetLowEnergyUsageLimit() || gammaekin0>GetHighEnergyUsageLimit() || gammaekin0<=electroncut) {
            //std::cout<<"fuori soglia!"<<std::endl;
            return numSecondaries;
        }
        // sample gamma energy
        // here we need 3 random number + 2 later for photon direction theta phi sampling
        double *rndArray = td->fDblArray;
        td->fRndm->uniform_array(5, rndArray);
        
        
        //
        double cosTheta = 0.0;
        double sinTheta = 0.0;
        
        //store old gamma directions in the lab frame
        double gamDirX0=track.GetDirX();
        double gamDirY0=track.GetDirY();
        double gamDirZ0=track.GetDirZ();
        
        //std::cout<<"\n\n********\nOriginal Directions: "<<gamDirX0<<" "<< gamDirY0 << " "<<gamDirZ0<<std::endl;
        
       // sample gamma energy and gamma scattering angle in the scattering frame i.e. which z-dir points to the original gamma direction
       //double gammaekin2  = SamplePhotonEnergy(matCut, gammaekin0,rndArray[0], rndArray[1], rndArray[2]);
        
       double gammaekin1  = SamplePhotonEnergyAndDirection(matCut, gammaekin0, sinTheta, cosTheta,td);
        
       //std::cout<<"gammaekin2: "<< gammaekin2<< ", gammaekin1: "<<gammaekin1<<std::endl;
        
       double phi      = geant::kTwoPi*(rndArray[4]);
       
        // new gamma direction in the scattering frame
        double gamDirX1  = sinTheta*std::cos(phi);
        double gamDirY1  = sinTheta*std::sin(phi);
        double gamDirZ1  = cosTheta;
        
        /*
        std::cout<<"cosTheta: "<<cosTheta<<std::endl;
        std::cout<<"sinTheta: "<<sinTheta<<std::endl;
        std::cout<<"Phi: "<<phi<<std::endl;
        std::cout<<"New Directions: "<<gamDirX1<<" "<< gamDirY1 << " "<<gamDirZ1<<std::endl;
        */
        
        
        // rotate new gamma direction to the lab frame:
        RotateToLabFrame(gamDirX1, gamDirY1, gamDirZ1, gamDirX0, gamDirY0, gamDirZ0);
        
        //std::cout<<"After rotations\n New Gamma Dir: "<<gamDirX1<<" "<< gamDirY1 << " "<<gamDirZ1<<std::endl;
        //std::cout<<"Original Gamma Dir: "<<gamDirX0<<" "<< gamDirY0 << " "<<gamDirZ0<<std::endl;
        
        double edep=0.;
        //electroncut=0.000002;
        //if(gammaekin1 >GetLowEnergyUsageLimit()) {
        //std::cout<<"if("<<gammaekin1<<" > "<<electroncut<<"): survive"<<std::endl;
        if(gammaekin1 > electroncut) {
            
            //std::cout<<"Survive\n";
            //Update primary track direction and energy
            track.SetDirX(gamDirX1);
            track.SetDirY(gamDirY1);
            track.SetDirZ(gamDirZ1);
            // update primary track kinetic energy
            track.SetKinE(gammaekin1);
            
        }
        else
        {
            //std::cout<<"Kill\n";
            track.SetTrackStatus(LTrackStatus::kKill);
            track.SetKinE(0.0);
            edep = gammaekin1;
        }
        
        double eDirX ;
        double eDirY ;
        double eDirZ ;
        double eekin = gammaekin0 - gammaekin1;
        
        //std::cout<<"electroncut: "<<electroncut<<". GetLowEnergyUsageLimit(): "<<GetLowEnergyUsageLimit()<<std::endl;
        
        //if the energy of the e- is greater than the cut, the secondary particle is effectively created
        //std::cout<<"if("<<eekin<<" > "<<electroncut<<"): create electron"<<std::endl;
        if(eekin > electroncut) {
            //std::cout<<"electron created\n";
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
            // this is known since it is a secondary track
            //  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
            std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
            
            // this is known since it is a secondary track
            //  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
            /*std::cout<<"Gamma original energy: "<< gammaekin0<< std::endl;
            std::cout<<"Gamma final energy:    "<< gammaekin1<< std::endl;
            std::cout<<"cosTheta:              "<< cosTheta<< std::endl;
            std::cout<<"sinTheta:              "<< sinTheta<< std::endl;
            std::cout<<"Electron energy:       "<< eekin<< std::endl;
*/
            sectracks[secIndx].SetDirX(eDirX*norm); //mb: controllare questa operazione
            sectracks[secIndx].SetDirY(eDirY*norm);
            sectracks[secIndx].SetDirZ(eDirZ*norm);
            //std::cout<<"Norm: "<<norm<<std::endl;
            /*std::cout<<"eDirX: "<<eDirX<<std::endl;
            std::cout<<"eDirY: "<<eDirY<<std::endl;
            std::cout<<"eDirZ: "<<eDirZ<<std::endl;
            std::cout<<"eDirX*norm: "<<eDirX*norm<<std::endl;
            std::cout<<"eDirY*norm: "<<eDirY*norm<<std::endl;
            std::cout<<"eDirZ*norm: "<<eDirZ*norm<<std::endl;
            std::cout<<"Module: "<<sqrt(eDirX*norm*eDirX*norm+eDirY*norm*eDirY*norm+eDirZ*norm*eDirZ*norm)<<std::endl;
            */
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
    
    
    
    KleinNishinaComptonModel::KleinNishinaComptonModel(int datafileindx, const std::string &modelname)
    : EMModel(modelname), /*fIsElectron(iselectron),*/ fDataFileIndx(datafileindx) {
        fDCSMaxZet                       = 0;
        fLoadDCSNumElectronEnergies      = 0;
        fLoadDCSReducedPhotonEnergyGrid  = 0;
        fLoadDCSElectronEnergyGrid       = nullptr;
        fLoadDCSReducedPhotonEnergyGrid  = nullptr;
        fLoadDCSForElements              = nullptr;
        
        
        fNumSamplingGammaEnergies = 71;// 171=>25 per decada//71; // between the min/max gamma kinetic energies ? mb: why this value?
        fNumSamplingPhotEnergies = 54; // at each energy grid points
        fMinGammaEnergy           =  1.0*geant::keV; // minimum kinetic energy of the interacting gamma
        fMaxGammaEnergy           = 10.0*geant::GeV; // maximum kinetic energy of the interacting gamma
        fGammaEnLMin                = 0.0;
        fGammaEnILDelta             = 1.0;
        fSamplingGammaEnergies    = nullptr;
        fLSamplingGammaEnergies   = nullptr;
        
        fNumMaterialCuts           = 0;
        fNumDifferentMaterialECuts = 0;
        fGlobalMatGCutIndxToLocal = nullptr;
        fAliasData                = nullptr; //alias data for each matrial-gammacut pairs
        fAliasSampler             = nullptr;
        
        fSecondaryInternalCode    = Electron::Definition()->GetInternalCode();
        
        //Initialise();
    }
    
    KleinNishinaComptonModel::~KleinNishinaComptonModel() {
        //if (fLoadDCSElectronEnergyGrid)
        //    delete [] fLoadDCSElectronEnergyGrid;
        
        //if (fLoadDCSReducedPhotonEnergyGrid)
        //    delete [] fLoadDCSReducedPhotonEnergyGrid;
        
        //if (fLoadDCSForElements) {
        //    for (int i=0; i<fDCSMaxZet; ++i)
        //        if (fLoadDCSForElements[i])
        //            delete [] fLoadDCSForElements[i];
        //    delete [] fLoadDCSForElements;
        //}
        
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
    
    
    
    void KleinNishinaComptonModel::Initialise() {
        LoadDCSData();
        InitSamplingTables();
        fAliasSampler          = new AliasTable();
        fSecondaryInternalCode = Electron::Definition()->GetInternalCode();
        //  std::cerr<<"  ----> SB model = " << GetName() << "  is initialized!"<< std::endl;
    }
    
    /**
     *  The sampling is based on the sampling tables prepared at initialisation. Statistical interpolation is used to
     *  select one of the incident particle kinetic energy grid points out of \f$ E_i \leq E_{kin} < E_{i+1}\f$ (linear
     *  interpolation in log kinetic energy) at a given primary particle kinetic energy \f$E_{kin}\f$. Then the transformed
     *  variable \f$\xi\in[0,1]\f$ is sampled from the sampling table (prepared at initialisation) that belongs to the
     *  selected incident particle kinetic energy grid point. The emitted photon energy \f$k\f$ then is obtained by
     *  applying the following transformation:
     *  \f[
     *     k = \sqrt{[k_c^2+k_p^2]\exp{ \left[ \xi\ln\frac{E_{kin}^2+k_p^2}{k_c^2+k_p^2}\right]}-k_p^2}
     *  \f]
     *  where \f$E_{kin}\f$ is the current incident particle (i.e. e-/e+) kinetic energy, \f$k_c\f$ is the current gamma
     *  particle kinetic energy production threshold and \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E/(mc^2)\f$,
     *  (\f$E\f$ is the total energy of the incident particle)
     *  \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$ (\f$n_e\f$ is the electron density of the current material and
     *  \f$r_e\f$ is the classical electron radius).
     */
    
    double KleinNishinaComptonModel::SamplePhotonEnergy(const MaterialCuts *matcut, double gammaekin, double r1, double r2, double r3){
        
        std::cout<<"TBI"<<std::endl;
        //constexpr double  mgdl  = 4.0*geant::kPi*geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght*geant::kRedElectronComptonWLenght;
        //constexpr double  ddum0 = 2.0*geant::kElectronMassC2;
        //constexpr double  ddum1 = geant::kElectronMassC2*geant::kElectronMassC2;
        
        //double densityFactor = matcut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*mgdl;
        //double ptot2         = gammaekin*(gammaekin+ddum0);
        //double etot2         = ptot2+ddum1;
        //double densityCor    = densityFactor*etot2;
        //double gcut          = matcut->GetProductionCutsInEnergy()[0];
        
        int mcindx       = matcut->GetIndex();
        int macindxlocal = fGlobalMatGCutIndxToLocal[mcindx];
        // the location of the first-gamma-energy lin-alias data of this mat-cut
        int indxstart    = macindxlocal*fNumSamplingGammaEnergies;
        // determine gamma energy lower grid point
        double leenergy  = std::log(gammaekin);
        int eenergyindx  = (int) ((leenergy-fGammaEnLMin)*fGammaEnILDelta);
        
        if (eenergyindx>=fNumSamplingGammaEnergies-1)
            eenergyindx = fNumSamplingGammaEnergies-2;
        
        double ploweener = (fLSamplingGammaEnergies[eenergyindx+1]-leenergy)*fGammaEnILDelta;
        
        indxstart +=eenergyindx;
        if (r1>ploweener|| !(fAliasData[indxstart]))
            ++indxstart;
        
        // sample the transformed variable
        double gammae = fAliasSampler->SampleLinear(fAliasData[indxstart]->fXdata, fAliasData[indxstart]->fYdata,
                                                    fAliasData[indxstart]->fAliasW, fAliasData[indxstart]->fAliasIndx,
                                                    fNumSamplingPhotEnergies,r2,r3);
        
        std::cout<<"gammae: "<<gammae<<std::endl;
        return gammae;

        
    }
    void KleinNishinaComptonModel::SamplePhotonDirection(double elenergy, double sint2, double onecost, double &sinTheta, double &cosTheta, double rndm){
    
    //TBI
    
    }
    
    
    //mb: added from geant4
    static const double
    d1= 2.7965e-1*geant::barn, d2=-1.8300e-1*geant::barn,
    d3= 6.7527   *geant::barn, d4=-1.9798e+1*geant::barn,
    e1= 1.9756e-5*geant::barn, e2=-1.0205e-2*geant::barn,
    e3=-7.3913e-2*geant::barn, e4= 2.7079e-2*geant::barn,
    f1=-3.9178e-7*geant::barn, f2= 6.8241e-5*geant::barn,
    f3= 6.0480e-5*geant::barn, f4= 3.0274e-4*geant::barn;
    static const int nlooplim = 1000;
    
    //mb: Sample photon energy
    //double KleinNishinaComptonModel::SamplePhotonEnergyAndDirection(const MaterialCuts *matcut, double gammaekin, double &sinTheta, double &cosTheta, double r1, double r2, double r3){
    double KleinNishinaComptonModel::SamplePhotonEnergyAndDirection(const MaterialCuts *matcut, double gammaekin, double &sinTheta, double &cosTheta, Geant::GeantTaskData *td){
        
        
        double *rndArray = td->fDblArray;
        
        
        double E0_m = gammaekin / geant::kElectronMassC2 ;
        //std::cout<<"gamEnergy0: "<<gammaekin<<std::endl;
        //std::cout<<"E0_m: "<<E0_m<<std::setprecision(8)<<std::endl;
        
        //
        // sample the energy rate of the scattered gamma
        //
        
        double epsilon, epsilonsq, greject, sint2, onecost ;
        double eps0       = 1./(1. + 2.*E0_m);
        double epsilon0sq = eps0*eps0;
        double alpha1     = - std::log(eps0);
        double alpha2     = alpha1 + 0.5*(1.- epsilon0sq);
        /*std::cout<<"eps0: "<<eps0<<std::endl;
        std::cout<<"epsilon0sq: "<<epsilon0sq<<std::endl;
        std::cout<<"alpha1: "<<alpha1<<std::endl;
        std::cout<<"alpha2: "<<alpha2<<std::endl;

        std::cout<<"rndm[0]: "<<r1<<std::endl;
        std::cout<<"rndm[1]: "<<r2<<std::endl;
        std::cout<<"rndm[2]: "<<r3<<std::endl;*/
        
        int nloop = 0;
        do {
            ++nloop;
            // false interaction if too many iterations
            if(nloop > nlooplim) {
                //std::cout<<"occhio"<<std::endl;
                return -1; }
            td->fRndm->uniform_array(3, rndArray);
            
            //if ( alpha1 > alpha2*r1 ) {
            if ( alpha1 > alpha2*rndArray[0] ) {
                epsilon   = std::exp(-alpha1*rndArray[1]);   // eps0**r
                epsilonsq = epsilon*epsilon;
                
            } else {
                epsilonsq = epsilon0sq + (1.- epsilon0sq)*rndArray[1];
                epsilon   = sqrt(epsilonsq);
            };
            
            onecost = (1.- epsilon)/(epsilon*E0_m);
            sint2   = onecost*(2.-onecost);
            greject = 1. - epsilon*sint2/(1.+ epsilonsq);
            /*std::cout<<"onecost: "<<onecost<<std::endl;
            std::cout<<"sint2: "<<sint2<<std::endl;
            std::cout<<"greject: "<<greject<<std::endl;*/
            
        } while (greject < rndArray[2]);
        //if(nloop>2)
        //std::cout<<nloop<<std::endl;
        
        //std::cout<<"AFTER LOOP\nonecost: "<<onecost<<std::endl;
        //std::cout<<"sint2: "<<sint2<<std::endl;

        
        // scattered gamma angles. ( Z - axis along the parent gamma)
        if(sint2 < 0.0) { sint2 = 0.0; }
        cosTheta = 1. - onecost;
        sinTheta = sqrt (sint2);
        double gamEnergy1=epsilon*gammaekin;
        /*std::cout<<"cosTheta: "<<cosTheta<<", sinTheta: "<<sinTheta<<". And Gammaout: "<<epsilon*gammaekin<<std::endl;
        std::cout<<"cosTeta: "<<cosTheta<<std::endl;
        std::cout<<"sinTeta: "<<sinTheta<<std::endl;
        
        std::cout<<"gamEnergy1: "<<gamEnergy1<<"*******\n\n\n"<<std::endl;*/
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
        for (int i=1; i<fNumSamplingGammaEnergies-1; ++i) {
            fLSamplingGammaEnergies[i] = fGammaEnLMin+i*delta;
            fSamplingGammaEnergies[i]  = std::exp(fGammaEnLMin+i*delta);
            //    std::cerr<<" E("<<i<<") = "<<fSamplingGammaEnergies[i]/geant::GeV<<std::endl;
        }  // fMinElecEnergy fMaxElecEnergy az fLoadDCSElectronEnergyGrid[0] es fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1]
        
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
        //std::cerr<<" === Number of local Material-Cuts = "<<fNumDifferentMaterialECuts<<std::endl;
        
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
                BuildOneLinAlias(ialias, theMaterialCutsTable[i]->GetMaterial(), theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]);
                isdone[localindx] = 1;
            }
        }
        delete [] isdone;
        // test
        //  for (int i=0; i<fNumMaterialCuts; ++i)
        //    std::cerr<<"     --> Global MatCut-indx = "<< i << " local indx = "<<fGlobalMatGCutIndxToLocal[i] <<std::endl;
    }

    
    /**
     *  The distribution of the energy of the emitted bremsstrahlung photons \f$ k \f$ is determined by the differential
     *  cross section
     *  \f[
     *    p(k) \propto \frac{\mathrm{d}\sigma}{\mathrm{d}k} \frac{1}{k\Gamma(k)}
     *  \f]
     *  where \f$ \Gamma(k) = (1+k_p^2/k^2) \f$ is the main factor related to dielectric suppression. The dielectric
     *  suppression removes the infrared divergence and strongly suppresses the low photon energy tail of the distribution
     *  which results in a non-ideal shape from sampling point of view. Furthermore, since sampling tables are built only
     *  at discrete \f$ E_{i}^{kin} \f$ incident particle kinetic energies and during the simulation we always have incident
     *  particle with \f$ E_{i}^{kin} < E^{kin} < E_{i+1}^{kin}\f$, we need to transform the prepared distributions to the
     *  common range i.e. [0,1] in order to guarante that the sampled photon energy is always within the proper kinematic
     *  limits i.e. \f$ k \in [k_c,E^{kin}] \f$ where \f$ k_c \f$ is the kinetic energy threshold for gamma particle
     *  production. So we apply the following variable transformations:
     *
     *  - since the Seltzer-Berger DCS are available in a "scalled" form i.e.
     *    \f[
     *         \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{Z^2}{\beta^2} \frac{1}{k}\chi(\kappa; E,Z)
     *    \f]
     *    we transform
     *    \f$ k \to \kappa(k) = k/E,\; \frac{\mathrm{d}k}{\mathrm{d}\kappa} = E \f$ and by using
     *    \f$p(k)\mathrm{d}k \propto p(\kappa)\mathrm{d}\kappa \f$ one can get
     *    \f$ p(\kappa) \propto p(k) \frac{\mathrm{d}k}{\mathrm{d}\kappa}
     *    \propto \chi(\kappa; E,Z) \frac{1}{\kappa} \frac{\kappa^2E^2}{\kappa^2E^2+k_p^2}\f$
     *  - then we can get rid of the dielectric suppression tail by transforming
     *    \f$ \kappa \to u(\kappa)=\ln(\kappa^2E^2+k_p^2),\; \kappa=\sqrt{\frac{e^u-k_p^2}{E^2}},\;
     *    \frac{\mathrm{d}\kappa}{\mathrm{d}u}=\frac{\kappa^2E^2+k_p^2}{2E^2\kappa} \f$ and one can use
     *    \f$ p(\kappa)\mathrm{d}\kappa \propto p(u)\mathrm{d}u \f$ to get
     *    \f$ p(u) \propto p(\kappa)\frac{\mathrm{d}\kappa}{\mathrm{d}u} \propto \chi(\kappa;E,Z) \f$
     *  - then we appaly the \f$ u \to \xi(u) = [u-\ln(k_c^2+k_p^2)]/\ln[(E^2+k_p^2)/(k_c^2+k_p^2)] \in [0,1]\f$
     *    transformation that by taking into account that \f$ u = \xi\ln[(E^2+k_p^2)/(k_c^2+k_p^2)]+\ln[k_c^2+k_p^2],\;
     *    \frac{\mathrm{d}\xi}{\mathrm{d}u} = \ln[(E^2+k_p^2)/(k_c^2+k_p^2)]\f$ and by using \f$ p(u)\mathrm{d}u \propto
     *    p(\xi)\mathrm{d}\xi\f$ one can get \f$ p(\xi) \propto p(u)\frac{\mathrm{d}u}{\mathrm{d}\xi} \propto
     *    \ln[(E^2+k_p^2)/(k_c^2+k_p^2)] \chi(\kappa;E,Z) \propto \chi(\kappa;E,Z)\f$ since
     *    \f$ \ln[(E^2+k_p^2)/(k_c^2+k_p^2)]\f$ is just a constant.
     *
     *  When this transformed p.d.f. are prepared the numerical "scalled" Seltzer-Berger DCS are interpolated similarly
     *  like in the case of SeltzerBergerBremsModel::ComputeXSectionPerAtom. The variable \f$ \phi=\ln[1-\kappa+10^{-12}]\f$
     *  can be obtained at a given \f$ \xi \f$ value by \f$ \kappa = \sqrt{\frac{e^u-k_p^2}{E^2}}
     *  =\frac{1}{E}\sqrt{[k_c^2+k_p^2]e^{\xi\ln[(E^2+k_p^2)/(k_c^2+k_p^2)]}-k_p^2} \f$ and during the sampling, the emitted
     *  gamma photon energy \f$k = \kappa E = \sqrt{[k_c^2+k_p^2]e^{\xi\ln[(E^2+k_p^2)/(k_c^2+k_p^2)]}-k_p^2} \f$  for a
     *  given sampled transformed variable value \f$\xi\f$ (where \f$E\f$ is the actual run-time primary kinetic energ)y so
     *  the sampled photon energy will always be within the actual kinematic limits i.e. \f$ k_c<k \leq E\f$) .
     */
 
    
    void KleinNishinaComptonModel::BuildOneLinAlias(int ialias, const Material *mat, double gcut){
        // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
        double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
        // we will need the element composition of this material
        const std::vector<Element*> theElements = mat->GetElementVector();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        //  int   numElems = theElements.size();
        
        //
        // for the common particle kinetic energy grid points
        //
        for (int ieener=0; ieener<fNumSamplingGammaEnergies; ++ieener) {
            double eener = fSamplingGammaEnergies[ieener];
            if (eener>gcut) { // otherwise no gamma production at that e- energy so let the sampling table to be nullptr
                //      std::cerr<<"         indx ="<<ialias<<std::endl;
                // find the e- energy index in the available e- energies grid such that we are just above
                int eenerindx = 0;
                for (; eenerindx<fLoadDCSNumElectronEnergies; ++eenerindx)
                    if (eener<fLoadDCSElectronEnergyGrid[eenerindx])
                        break;
                if (eenerindx==fLoadDCSNumElectronEnergies)
                    --eenerindx;
                //eenerindx is the upper index now
                
                double ptot2   = eener*(eener+2.0*geant::kElectronMassC2);
                double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
                //      double ibeta2  = etot2/ptot2;
                
                double densityCor = densityFactor*etot2;
                
                double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies];
                double *logReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
                // ln(x)-ln(x1)
                double dum0 = std::log(eener/fLoadDCSElectronEnergyGrid[eenerindx-1]);
                // ln(x2)-ln(x1)
                double dum1 = std::log(fLoadDCSElectronEnergyGrid[eenerindx]/fLoadDCSElectronEnergyGrid[eenerindx-1]);
                
                for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
                    logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
                    int indxdcsh = eenerindx*fLoadDCSNumReducedPhotonEnergies + irpener;
                    int indxdcsl = (eenerindx-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
                    theDCS[irpener] = 0.0;
                    for (unsigned long ielem=0; ielem<theElements.size(); ++ielem) {
                        double zet  = theElements[ielem]->GetZ();
                        int    izet = lrint(zet);
                        // ln(y2) -ln(y1)
                        double dum2 = std::log(fLoadDCSForElements[izet-1][indxdcsh]/fLoadDCSForElements[izet-1][indxdcsl]);
                        double dcs  = dum2/dum1*dum0+std::log(fLoadDCSForElements[izet-1][indxdcsl]); //this is ln(dcs)
                        dcs = std::exp(dcs);
                        dcs *= theAtomicNumDensityVector[ielem]*zet*zet;
                        // correction for positrons
                        //if (!fIsElectron) {
                            //dcs *= PositronCorrection1(eener, fLoadDCSReducedPhotonEnergyGrid[irpener], gcut, zet);
                            //             dcs *= PositronCorrection(eener, ibeta2, fLoadDCSReducedPhotonEnergyGrid[irpener], zet);
                        //}
                        theDCS[irpener] += dcs;//dcs/(fLoadDCSReducedPhotonEnergyGrid[irpener]+1e-12);
                    }
                    //         if (theDCS[irpener]<=0.0) // e+ dcs is zero at kappa=1
                    //           theDCS[irpener] = 1.0e-13;
                    //         theDCS[irpener] = std::log(theDCS[irpener]);
                }
                
                // set up a spline on the log(1-kappa)::log(dcs)
                Spline     *sp = new Spline();
                sp->SetUpSpline(logReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);
                
                //
                // fill in the initial x,y data
                //
                // create the alias data struct
                fAliasData[ialias] = new LinAlias();
                fAliasData[ialias]->fXdata     = new double[fNumSamplingPhotEnergies]();
                fAliasData[ialias]->fYdata     = new double[fNumSamplingPhotEnergies]();
                fAliasData[ialias]->fAliasW    = new double[fNumSamplingPhotEnergies]();
                fAliasData[ialias]->fAliasIndx = new int[fNumSamplingPhotEnergies]();
                // find reduced foton energy index such that we are just above
                double kappac = gcut/eener;
                int kappaindx = 0;
                for (; kappaindx<fLoadDCSNumReducedPhotonEnergies; ++kappaindx)
                    if (kappac<fLoadDCSReducedPhotonEnergyGrid[kappaindx])
                        break;
                if (std::abs(1.0-kappac/fLoadDCSReducedPhotonEnergyGrid[kappaindx])<1.e-12) // handle some possible rounding problems: if kappac is very close
                    ++kappaindx;
                
                if (kappaindx>=fLoadDCSNumReducedPhotonEnergies)
                    kappaindx = fLoadDCSNumReducedPhotonEnergies-1;
                
                // fill the first initial value; convert kappa scale to 0-1
                int numdata = 1;
                fAliasData[ialias]->fXdata[0] = 0.0;
                double kappaconv = std::log(1.0-kappac+1.0e-12);
                //       fAliasData[ialias]->fYdata[0] = std::exp(sp->GetValueAt(kappaconv));
                fAliasData[ialias]->fYdata[0] = sp->GetValueAt(kappaconv);
                
                
                // fill possible values if any
                for (int k=kappaindx; k<fLoadDCSNumReducedPhotonEnergies-1; ++k) {
                    double thekappa = fLoadDCSReducedPhotonEnergyGrid[k];
                    kappaconv = std::log( (thekappa*thekappa*eener*eener+densityCor)/(gcut*gcut+densityCor) )/
                    std::log( (eener*eener+densityCor)/(gcut*gcut+densityCor)); // coverted to 0-1
                    fAliasData[ialias]->fXdata[numdata] = kappaconv;  // xi
                    //         fAliasData[ialias]->fYdata[numdata] = std::exp(theDCS[k]);
                    fAliasData[ialias]->fYdata[numdata] = theDCS[k];
                    
                    ++numdata;
                }
                // and the last point
                fAliasData[ialias]->fXdata[numdata] = 1.0;
                //       fAliasData[ialias]->fYdata[numdata] = std::exp(theDCS[fLoadDCSNumReducedPhotonEnergies-1]);
                fAliasData[ialias]->fYdata[numdata] = theDCS[fLoadDCSNumReducedPhotonEnergies-1];
                ++numdata;
                
                // expand the data up to maximum
                while(numdata<fNumSamplingPhotEnergies) {
                    // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
                    double maxerr     = 0.0; // value of the current maximum error
                    double thexval    = 0.0;
                    double theyval    = 0.0;
                    int    maxerrindx = 0;   // the lower index of the corresponding bin
                    
                    for (int i=0; i<numdata-1; ++i) {
                        double xx = 0.5*(fAliasData[ialias]->fXdata[i]+fAliasData[ialias]->fXdata[i+1]); // mid point
                        double yy = 0.5*(fAliasData[ialias]->fYdata[i]+fAliasData[ialias]->fYdata[i+1]); // lin func val at the mid point
                        
                        double dum0  = (gcut*gcut+densityCor);
                        double dum1  = (eener*eener+densityCor)/dum0;
                        double    u  = dum0*std::exp(xx*std::log(dum1));
                        double thekappa = std::sqrt( u-densityCor)/eener; // kappa
                        
                        double conv  = std::log(1.0-thekappa+1.0e-12);
                        //           double spval = std::exp(sp->GetValueAt(conv)); // spline intp. val. at mid point
                        double spval = sp->GetValueAt(conv); // spline intp. val. at mid point
                        double err   = std::fabs(yy-spval);// works better than fabs(1-yy/spval) might try constrained spline?
                        if (err>maxerr) {
                            maxerr     = err;
                            maxerrindx = i;
                            thexval    = xx;
                            theyval    = spval;
                        }
                    }
                    // extend x,y data by puting a spline interp.ted value at the mid point of the highest error bin
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
                }
                
                // set up a linear alias smapler on this data
                AliasTable *alst = new AliasTable();
                alst->PreparLinearTable(fAliasData[ialias]->fXdata, fAliasData[ialias]->fYdata,
                                        fAliasData[ialias]->fAliasW, fAliasData[ialias]->fAliasIndx,
                                        fNumSamplingPhotEnergies);
                
                delete alst;
                delete sp;
                delete [] theDCS;
                delete [] logReducedPhotonEnergyGrid;
            }
            ++ialias;
        }
    }
    
    
    /**
     * The restricted atomic cross section for bremsstrahlung photon emission for target element with atomic number
     * \f$Z\f$, gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according to
     * \f[
     *   \sigma(E;k_c,Z) = \int_{k_c}^{E} \frac{1}{\Gamma}\frac{\mathrm{d}\sigma}{\mathrm{d}k}\mathrm{d}k
     * \f]
     * if \f$E>k_c\f$ and immediate return with \f$0\f$ otherwise. (The \f$1/\Gamma\f$ factor is the main dielectric
     * suppression factor and \f$\Gamma = (1+k_p^2/k^2)\f$ where
     * \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E_t/(mc^2)\f$, \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$
     * (\f$n_e\f$ is the electron density) see more details at RelativisticBremsModel::ComputeURelDXSecPerAtom ).
     * The Seltzer-Berger numerical DCS are available in the from of "scalled" DCS as
     * \f[
     *      \chi(\kappa;E,Z) =  \frac{\beta^2}{Z^2}k \frac{\mathrm{d}\sigma}{\mathrm{d}k}
     * \f]
     * where \f$k\f$ is the emitted photon energy, \f$\kappa=k/E\f$ is the reduced photon energy. The above integral can be
     * written now with the "scalled" DCS as
     * \f[
     *   \sigma(E;k_c,Z) = \int_{k_c}^{E} \frac{1}{\Gamma}\frac{\mathrm{d}\sigma}{\mathrm{d}k}\mathrm{d}k =
     *                     \frac{Z^2}{\beta^2}\int_{\kappa_c}^{1} \frac{1}{\kappa\Gamma}\chi(\kappa;E,Z)\mathrm{d}\kappa
     * \f]
     * where \f$\kappa_c=k_c/E\f$.
     * Since the "scalled" DCS are available at fixed electron kinetic energies \f$\{E_i\}_i^N\f$, linear interpolation in
     * log electron energies is applied to obtain \f$\chi(\kappa;E,Z)\f$ from \f$\chi(\kappa;E_i,Z)\f$ and
     * \f$\chi(\kappa;E_{i+1},Z)\f$ such that \f$E_i \leq E < E_{i+1}\f$ over the available fixed set of reduced photon
     * energies \f$\{\kappa_j\}_j^M\f$ where \f$ \kappa_j \in [0,1]\; \forall\; j\f$. During the interpolation, the reduced
     * photon energy grid i.e. \f$\{\kappa_j\}_j^M\f$ is transformed to \f$ \phi = \ln[1-\kappa+10^{-12}] \f$ for getting a
     * more accurate interpolation later when the integral is computed (with 64-points Gauss-Legendre quadrature using cubic
     * spline interpolation of DCS values over the \f$\phi\f$ grid).
     *
     * The integral is computed by 64-points Gauss-Legendre quadrature after applying the following transformations
     * - first the reduced photon energy is transformed to \f$\kappa \to u=\ln(\kappa) \f$
     * - then we apply the following transformation \f$ u \to \xi = [u-\ln(\kappa_c)]/\ln(1/\kappa_c) \in [0,1] \f$
     *
     * The transformed integral
     * \f[
     *   \sigma(E;k_c,Z) = \frac{Z^2}{\beta^2}\int_{\kappa_c}^{1} \frac{1}{\kappa\Gamma}\chi(\kappa;E,Z)\mathrm{d}\kappa
     *                   = \frac{Z^2}{\beta^2}\int_{\ln(\kappa_c)}^{0} \frac{1}{\Gamma}\chi(\kappa;E,Z)\mathrm{d}u
     *                   = \frac{Z^2}{\beta^2}\ln\frac{1}{\kappa_c}\int_{0}^{1} \frac{1}{\Gamma}\chi(\phi;E,Z)\mathrm{d}\xi
     * \f]
     * where \f$\Gamma\f$ must be evaluated at \f$k=E\kappa_c e^{\xi\ln(1/\kappa_c)}\f$ and \f$ \chi(\phi;E,Z) \f$ must be
     * evaluated (interpolated) at \f$\phi =\ln[1-\kappa_c e^{\xi\ln(1/\kappa_c)}+10^{-12}]\f$ at a given value of \f$\xi\f$.
     */
    
    
    
    //mb: ok modified
    /*
    double KleinNishinaComptonModel::ComputeXSectionPerAtom(const Element *elem, const Material *mat, double electronprodcutenergy,
                                                            double gammaekin){
        double xSection = 0.0;
        static const double a = 20.0 , b = 230.0 , c = 440.0;
        double Z=elem->GetZ();
        
        double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
        p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);
        
        double T0  = 15.0*geant::keV;
        if (Z < 1.5) { T0 = 40.0*geant::keV; }
        
        double X   = std::max(gammaekin, T0) / geant::kElectronMassC2;
        xSection = p1Z * std::log(1.+2.*X)/X
        + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
        
        //  modification for low energy. (special case for Hydrogen)
        if (gammaekin < T0) {
            static const double dT0 = geant::keV;
            X = (T0+dT0) / geant::kElectronMassC2 ;
            double sigma = p1Z*std::log(1.+2*X)/X
            + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
            double   c1 = -T0*(sigma-xSection)/(xSection*dT0);
            double   c2 = 0.150;
            if (Z > 1.5) { c2 = 0.375-0.0556*std::log(Z); }
            double    y = std::log(gammaekin/T0);
            xSection *= std::exp(-y*(c1+c2*y));
        }
        // G4cout<<"e= "<< GammaEnergy<<" Z= "<<Z<<" cross= " << xSection << G4endl;
        return xSection;
        /*
        
        /*
         if (electronekin<=gammaprodcutenergy || electronekin<fLoadDCSElectronEnergyGrid[0]
         || electronekin>fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1])
         return xsec;
         
         // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
         double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
         *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
         
         double zet     = elem->GetZ();
         int    izet    = std::lrint(zet);
         double ptot2   = electronekin*(electronekin+2.0*geant::kElectronMassC2);
         double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
         double ibeta2  = etot2/ptot2;
         double densityCor = densityFactor*etot2;//electronekin*electronekin;  // this is k_p^2
         double kappacr    = gammaprodcutenergy/electronekin;
         double logikappacr = std::log(1./kappacr);
         
         // find the electron energy index that
         int ieener = 0;
         for (ieener=0; ieener<fLoadDCSNumElectronEnergies; ++ieener)
         if (electronekin<fLoadDCSElectronEnergyGrid[ieener])
         break;
         if (ieener==fLoadDCSNumElectronEnergies)
         --ieener;
         
         double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies];
         double *logReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
         // ln(x)-ln(x1)
         double dum0 = std::log(electronekin/fLoadDCSElectronEnergyGrid[ieener-1]);
         // ln(x2)-ln(x1)
         double dum1 = std::log(fLoadDCSElectronEnergyGrid[ieener]/fLoadDCSElectronEnergyGrid[ieener-1]);
         //  std::cerr<<" -->"<<electronekin<<" "<<ieener<<std::endl;
         
         for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
         logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
         int indxdcsh = ieener*fLoadDCSNumReducedPhotonEnergies + irpener;
         int indxdcsl = (ieener-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
         double dcsl = fLoadDCSForElements[izet-1][indxdcsl];
         double dcsh = fLoadDCSForElements[izet-1][indxdcsh];
         // ln(y2) -ln(y1)
         double dum2 = std::log(dcsh/dcsl);
         double dcs  = dum2/dum1*dum0+std::log(dcsl);
         dcs  = std::exp(dcs);
         // correction for positrons
         //    if (!fIsElectron)
         //      dcs *= PositronCorrection(electronekin, ibeta2, fLoadDCSReducedPhotonEnergyGrid[irpener], zet);
         theDCS[irpener] = dcs;
         }
         
         // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
         int ngl = 64;
         GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
         const std::vector<double> glW = gl->GetWeights();
         const std::vector<double> glX = gl->GetAbscissas();
         // we will need a natural cubic spile for the integral
         Spline     *sp = new Spline();
         sp->SetUpSpline(logReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);
         
         double integral = 0.0;
         for (int i=0;i<ngl;++i) {
         double dumx = 1.0-std::exp(glX[i]*logikappacr)*kappacr; // ez a felso x -nek az exp(x)-e
         double x = std::log(dumx+1.0e-12);
         double egamma = (1.0-dumx)*electronekin;
         double poscor = 1.0;
         if (!fIsElectron)
         poscor *= PositronCorrection(electronekin, ibeta2, 1.0-dumx, zet);
         
         integral+= glW[i]*poscor*sp->GetValueAt(x)/(1.+densityCor/(egamma*egamma));
         }
         
         delete [] theDCS;
         delete [] logReducedPhotonEnergyGrid;
         delete sp;
         delete gl;
         
         return logikappacr*zet*zet*ibeta2*integral;*/
   // }
    
    /**
     *   The restricted macroscopic cross section for bremsstrahlung photon emission for the given target material,
     *   gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according to
     *  \f[
     *      \Sigma(E;k_c,\mathrm{material}) = \sum_i n_i \sigma_i(E;k_c,Z_i)
     *  \f]
     *  if \f$ E>k_c\f$ otherwise immediate return with \f$0\f$. The summation goes over the elements the matrial is
     *  composed from. \f$\sigma_i(E;k_c,Z_i)\f$ is the restricted atomic cross
     *  secion for the \f$i\f$-th element of the material with atomic number of \f$Z_i \f$ (computed similarly like
     *  SeltzerBergerBremsModel::ComputeXSectionPerAtom()) and \f$n_i\f$ is the number of atoms per unit volume of
     *  \f$i \f$-th element of the material that is \f$ n_i = \mathcal{N}\rho w_i/A_i \f$ where \f$\mathcal{N}\f$ is the
     *  Avogadro number, \f$\rho\f$ is the material density, \f$w_i\f$ is the proportion by mass (or mass fraction) of the
     *  \f$i \f$-th element and \f$A_i\f$ is the molar mass of the \f$i \f$-th element. The corresponding mean free path
     *  is \f$\lambda = 1/\Sigma \f$.
     */
    
    /*
    double KleinNishinaComptonModel::ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin){
        double xsec = 0.0;
        if (electronekin<=gammaprodcutenergy || electronekin<fLoadDCSElectronEnergyGrid[0]
            || electronekin>fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1])
            return xsec;
        
        // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
        double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
        
        //  double zet     = elem->GetZ();
        //  int    izet    = std::lrint(zet);
        double ptot2   = electronekin*(electronekin+2.0*geant::kElectronMassC2);
        double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
        double ibeta2  = etot2/ptot2;
        double densityCor = densityFactor*etot2;//electronekin*electronekin;  // this is k_p^2
        double kappacr    = gammaprodcutenergy/electronekin;
        double logikappacr = std::log(1./kappacr);
        
        // find the electron energy index that
        int ieener = 0;
        for (ieener=0; ieener<fLoadDCSNumElectronEnergies; ++ieener)
            if (electronekin<fLoadDCSElectronEnergyGrid[ieener])
                break;
        if (ieener==fLoadDCSNumElectronEnergies)
            --ieener;
        
        
        // we will need the element composition of this material
        const std::vector<Element*> theElements = mat->GetElementVector();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        int   numElems = theElements.size();
        
        double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies*numElems];
        double *logReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
        // ln(x)-ln(x1)
        double dum0 = std::log(electronekin/fLoadDCSElectronEnergyGrid[ieener-1]);
        // ln(x2)-ln(x1)
        double dum1 = std::log(fLoadDCSElectronEnergyGrid[ieener]/fLoadDCSElectronEnergyGrid[ieener-1]);
        
        for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
            logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
            int indxdcsh = ieener*fLoadDCSNumReducedPhotonEnergies + irpener;
            int indxdcsl = (ieener-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
            for (unsigned long ielem=0; ielem<theElements.size(); ++ielem) {
                double zet  = theElements[ielem]->GetZ();
                int    izet = lrint(zet);
                // ln(y2) -ln(y1)
                double dum2 = std::log(fLoadDCSForElements[izet-1][indxdcsh]/fLoadDCSForElements[izet-1][indxdcsl]);
                double dcs  = dum2/dum1*dum0+std::log(fLoadDCSForElements[izet-1][indxdcsl]); //this is ln(dcs)
                //      //dcs  = std::exp(dcs);
                //      dcs += std::log(theAtomicNumDensityVector[ielem]*zet*zet);
                dcs  = std::exp(dcs);
                dcs *= theAtomicNumDensityVector[ielem]*zet*zet;
                
                // correction for positrons
                //      if (!fIsElectron) {
                //        dcs *= PositronCorrection(electronekin, ibeta2, fLoadDCSReducedPhotonEnergyGrid[irpener], zet);
                //      }
                theDCS[ielem*fLoadDCSNumReducedPhotonEnergies+irpener] = dcs;
            }
            //    theDCS[irpener] = std::log(theDCS[irpener]);
        }
        
        // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
        int ngl = 64;
        GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
        const std::vector<double> glW = gl->GetWeights();
        const std::vector<double> glX = gl->GetAbscissas();
        // we will need a natural cubic spile for the integral
        //  Spline     *sp = new Spline();
        //  sp->SetUpSpline(logReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);
        // we will need as many natural cubic spile as elements
        Spline  **sp = new Spline*[numElems];
        for (int i=0; i<numElems; ++i) {
            sp[i] = new Spline();
            sp[i]->SetUpSpline(logReducedPhotonEnergyGrid, &theDCS[i*fLoadDCSNumReducedPhotonEnergies], fLoadDCSNumReducedPhotonEnergies);
        }
        
        double integral = 0.0;
        for (int i=0;i<ngl;++i) {
            double dumx = 1.0-std::exp(glX[i]*logikappacr)*kappacr; // ez a felso x -nek az exp(x)-e
            double x = std::log(dumx+1.0e-12);
            double egamma = (1.0-dumx)*electronekin;
            double sum = 0.0;
            for (int ielem=0; ielem<numElems; ++ielem) {
                //      double val = std::exp(sp[ielem]->GetValueAt(x));
                double val = sp[ielem]->GetValueAt(x);
                if (!fIsElectron) {
                    double zet  = theElements[ielem]->GetZ();
                  //mb  val *= PositronCorrection(electronekin, ibeta2, (1.0-dumx), zet);
                }
                sum += val;
            }
            integral+= glW[i]*sum/(1.+densityCor/(egamma*egamma));
        }
        
        delete [] theDCS;
        delete [] logReducedPhotonEnergyGrid;
        //  delete sp;
        for (int i=0; i<numElems; ++i)
            delete sp[i];
        delete [] sp;
        delete gl;
        return logikappacr*ibeta2*integral;
    }
     /*
    
    /**
     *  The stopping power, i.e. average energy loss per unit path length, from bremsstrahlung photon emission is computed
     *  for the given e-/e+ kinetic energy \f$E\f$, the given material and gamma production threshold energy \f$k_c\f$
     *  \f[
     *      S(E;k_c,\mathrm{material})=\int_{0}^{\eta} k \sum_i n_i \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \mathrm{d}k
     *  \f]
     *  where
     *  \f[
     *     \eta =
     *      \begin{cases}
     *            E    & \quad \mathrm{if}\; E<k_c \\
     *            k_c  & \quad \mathrm{if}\; E \geq k_c
     *      \end{cases}
     *  \f]
     *  the summation goes over the elements the matrial is composed from. \f$ \mathrm{d}\sigma_i \mathrm{d}k\f$ is the
     *  differential cross section for bremsstrahlung photon emission for for the \f$i\f$-th element of the material with
     *  atomic number of \f$Z_i\f$ and \f$n_i\f$ is the number of atoms per unit volume of \f$i\f$-th element of the
     *  material that is \f$n_i=\mathcal{N}\rho w_i/A_i\f$ where \f$\mathcal{N}\f$ is the Avogadro number, \f$\rho\f$ is
     *  the material density, \f$w_i\f$ is the proportion by mass (or mass fraction) of the \f$i\f$-th element and \f$A_i\f$
     *  is the molar mass of the \f$i\f$-th element.
     *
     *  The Seltzer-Berger atomic DCS are used and interpolated similarly like in case of atomic cross section computation
     *  (SeltzerBergerBremsModel::ComputeXSectionPerAtom).  The Seltzer-Berger numerical atomic DCS are available in
     *  the from of "scalled" DCS as
     *  \f[
     *      \chi(\kappa;E,Z) =  \frac{\beta^2}{Z^2}k \frac{\mathrm{d}\sigma}{\mathrm{d}k}
     *  \f]
     *  where \f$k\f$ is the emitted photon energy, \f$\kappa=k/E \in [0,1]\f$ is the reduced photon energy. The above
     *  integral can be written now with the "scalled" DCS as
     *  \f[
     *      S(E;k_c,\mathrm{material})= \int_{0}^{\eta} k \sum_i n_i \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \mathrm{d}k
     *    = \frac{E}{\beta^2}\int_{0}^{\eta/E}\frac{1}{\Gamma}\sum_i n_i Z_i^2 \chi(\kappa;E,Z_i) \mathrm{d}\kappa
     *  \f]
     *  (The \f$1/\Gamma\f$ factor is the main dielectric
     *  suppression factor and \f$\Gamma = (1+k_p^2/k^2)\f$ where
     *  \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E_t/(mc^2)\f$, \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$
     *  (\f$n_e\f$ is the electron density) see more details at RelativisticBremsModel::ComputeURelDXSecPerAtom ).
     *
     *  The integral is computed by 64-points Gauss-Legendre quadrature after the following transformation
     *  - the reduced photon energy is transformed \f$ \kappa \to \xi = \kappa/(\eta/E) \in [0,1] \f$.
     *
     *  The integral then becomes
     *  \f[
     *   S(E;k_c,\mathrm{material})
     *    = \frac{E}{\beta^2}\int_{0}^{\eta/E}\frac{1}{\Gamma}\sum_i n_i Z_i^2 \chi(\kappa;E,Z_i) \mathrm{d}\kappa
     *    = \frac{\eta}{\beta^2}\int_{0}^{1}\frac{1}{\Gamma}\sum_i n_i Z_i^2 \chi(\kappa;E,Z_i) \mathrm{d}\xi
     *  \f]
     *  where \f$\Gamma\f$ must be evaluated at \f$k=\xi\eta\f$ and \f$\chi\f$ at \f$\kappa=\xi\eta/E\f$ at a given value of
     *  \f$\xi\f$.
     */
    /*
    double KleinNishinaComptonModel::ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin){
        double dedx = 0.0;
        
        if (electronekin<fLoadDCSElectronEnergyGrid[0] || electronekin>fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1])
            return dedx;
        
        // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
        double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
        
        double ptot2   = electronekin*(electronekin+2.0*geant::kElectronMassC2);
        double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
        double ibeta2  = etot2/ptot2;
        double densityCor = densityFactor*etot2;//*electronekin*electronekin;  // this is k_p^2
        
        // find the index of the electron energy that is just below the requested one
        int ieener = 0;
        for (ieener=0; ieener<fLoadDCSNumElectronEnergies; ++ieener)
            if (electronekin<fLoadDCSElectronEnergyGrid[ieener])
                break;
        // handle the case when electronekin=fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1]
        if (ieener==fLoadDCSNumElectronEnergies) {
            --ieener;
        }
        
        // we will need the element composition of this material
        const std::vector<Element*> theElements = mat->GetElementVector();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        int   numElems = theElements.size();
        
        double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies*numElems];
        // ln(x)-ln(x1)
        double dum0 = std::log(electronekin/fLoadDCSElectronEnergyGrid[ieener-1]);
        // ln(x2)-ln(x1)
        double dum1 = std::log(fLoadDCSElectronEnergyGrid[ieener]/fLoadDCSElectronEnergyGrid[ieener-1]);
        for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
            int indxdcsh = ieener*fLoadDCSNumReducedPhotonEnergies + irpener;
            int indxdcsl = (ieener-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
            for (unsigned long ielem=0; ielem<theElements.size(); ++ielem) {
                double zet  = theElements[ielem]->GetZ();
                int    izet = lrint(zet);
                // ln(y2) -ln(y1)
                double dum2 = std::log(fLoadDCSForElements[izet-1][indxdcsh]/fLoadDCSForElements[izet-1][indxdcsl]);
                double dcs  = dum2/dum1*dum0+std::log(fLoadDCSForElements[izet-1][indxdcsl]); //this is ln(dcs)
                dcs  = std::exp(dcs);
                dcs *= theAtomicNumDensityVector[ielem]*zet*zet;
                //      dcs += std::log(theAtomicNumDensityVector[ielem]*zet*zet);
                theDCS[ielem*fLoadDCSNumReducedPhotonEnergies+irpener] = dcs;
            }
        }
        
        // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
        int ngl = 64;
        GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
        const std::vector<double> glW = gl->GetWeights();
        const std::vector<double> glX = gl->GetAbscissas();
        // we will need a natural cubic spile for the integral
        //  Spline     *sp = new Spline();
        //  sp->SetUpSpline(fLoadDCSReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);
        // we will need as many natural cubic spile as elements
        Spline  **sp = new Spline*[numElems];
        for (int i=0; i<numElems; ++i) {
            sp[i] = new Spline();
            sp[i]->SetUpSpline(fLoadDCSReducedPhotonEnergyGrid, &theDCS[i*fLoadDCSNumReducedPhotonEnergies], fLoadDCSNumReducedPhotonEnergies);
        }
        
        // integrate
        double integral    = 0.0;
        double upperlimit  = gammaprodcutenergy;
        if (upperlimit>electronekin)
            upperlimit = electronekin;
        double kappacr = upperlimit/electronekin;
        for (int i=0; i<ngl; ++i) {
            double t = glX[i]*kappacr;          // kappa
            double egamma = glX[i]*upperlimit;
            double sum = 0.0;
            for (int ielem=0; ielem<numElems; ++ielem) {
                //      double val = std::exp(sp[ielem]->GetValueAt(t));
                double val = sp[ielem]->GetValueAt(t);
                if (!fIsElectron) {
                    double zet  = theElements[ielem]->GetZ();
                 // mb   val *= PositronCorrection(electronekin, ibeta2, egamma/electronekin, zet);
                }
                sum += val;
            }
            
            integral += glW[i]*sum/(1.+densityCor/(egamma*egamma));
            // x 1/(1+k_p^2/k^2) i.e. density effect correction
        }
        
        delete [] theDCS;
        //  delete sp;
        for (int i=0; i<numElems; ++i)
            delete sp[i];
        delete [] sp;
        delete gl;
        
        return upperlimit*ibeta2*integral;
    }/*
    
/*
    // correction for positrons : DCS must be multiplied by this for positrons
    // ephoton is the reduced photon energy
    double KleinNishinaComptonModel::PositronCorrection(double ekinelectron, double ibeta2electron,
                                                        double ephoton, double z) {
        using geant::kElectronMassC2;
        constexpr double dum1 = geant::kTwoPi*geant::kFineStructConst;
        
        double poscor = 0.0;
        double ibeta1   = std::sqrt(ibeta2electron);
        double e2       = ekinelectron * (1.0 - ephoton);
        if (e2 > 0.0) {
            double ibeta2 = (e2 + kElectronMassC2)/std::sqrt(e2*(e2+2.0*kElectronMassC2));
            double dum0   = dum1*z*(ibeta1-ibeta2);
            if (dum0<-12.0) {
                poscor = 0.0;
            } else {
                poscor = std::exp(dum0);
            }
        } else {
            poscor = 0.0;
        }
        return poscor;
    }
 */
    
  /*
    
    //ephoton is the reduced photon energy
    double KleinNishinaComptonModel::PositronCorrection1(double ekinelectron, double ephoton, double gcutener, double z) {
        using geant::kElectronMassC2;
        constexpr double dum1 = geant::kTwoPi*geant::kFineStructConst;
        
        double poscor = 0.0;
        double e1     = ekinelectron-gcutener;  // here is the dif.
        double ibeta1 = (e1+kElectronMassC2)/std::sqrt(e1*(e1+2.0*kElectronMassC2));
        double e2     = ekinelectron * (1.0 - ephoton);
        double ibeta2 = (e2+kElectronMassC2)/std::sqrt(e2*(e2+2.0*kElectronMassC2));
        double ddum   = dum1*z*(ibeta1-ibeta2);
        if (ddum<-12.0) {
            poscor = 0.0;
        } else {
            poscor = std::exp(ddum);
        }
        return poscor;
    }*/
    
    void KleinNishinaComptonModel::LoadDCSData() {
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
        const std::vector<Element*> theElements = Element::GetTheElementTable();
        //std::cout<<theElements;
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
    }
    
}   // namespace geantphysics
