
#include "SauterGavrilaPhotoElectricModel.h"

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
    
    void SauterGavrilaPhotoElectricModel::Initialize() {
        EMModel::Initialize();
        Initialise();

    }
    

    double SauterGavrilaPhotoElectricModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle*) const {
        double mine = (matcut->GetProductionCutsInEnergy())[0]; // gamma production cut in the given material-cuts
        return mine;
    }
    
    
    SauterGavrilaPhotoElectricModel::SauterGavrilaPhotoElectricModel(int datafileindx, const std::string &modelname)
    : EMModel(modelname){
        
        fNumSamplingGammaEnergies = 71;// 171=>25 per decada//71; // between the min/max gamma kinetic energies ? mb: why this value?
        fNumSamplingElectronEnergies = 54;              // at each energy grid points
        fMinGammaEnergy           =  0.007*geant::keV;  // minimum kinetic energy of the interacting gamma
        fMaxGammaEnergy           =  26.0*geant::MeV;   //10.0*geant::GeV; // maximum kinetic energy of the interacting gamma
        fGammaEnLMin                = 0.0;
        fGammaEnILDelta             = 1.0;
        fSamplingGammaEnergies    = nullptr;
        fLSamplingGammaEnergies   = nullptr;
        
        fNumMaterialCuts           = 0;
        fNumDifferentMaterialECuts = 0;
        fGlobalMatGCutIndxToLocal = nullptr;
        fAliasData                = nullptr; //alias data for each material-gammacut pairs
        fAliasSampler             = nullptr;
        
        fSecondaryInternalCode    = Electron::Definition()->GetInternalCode();
        
    }
    
    SauterGavrilaPhotoElectricModel::~SauterGavrilaPhotoElectricModel() {
        
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
    
    
    void SauterGavrilaPhotoElectricModel::Initialise() {
        InitSamplingTables();
        fAliasSampler          = new AliasTable();
        fSecondaryInternalCode = Electron::Definition()->GetInternalCode();
    }
    
//____________________
    //NB: cosTheta is supposed to contain the dirZ of the incoming photon -> used in case that gamma>5
    void SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Rejection(double gammaEnIn, double &sinTheta, double &cosTheta, double &phi, Geant::GeantTaskData *td){
        
        double *rndArray = td->fDblArray;
        td->fRndm->uniform_array(1, rndArray);
        phi     = geant::kTwoPi * rndArray[0];
        
        double tau = gammaEnIn/geant::kElectronMassC2;
        static const double taulimit = 50.0;
        
        if (tau > taulimit) {
            sinTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
        
        } else
        {
            
            // algorithm according Penelope 2008 manual and
            // F.Sauter Ann. Physik 9, 217(1931); 11, 454(1931).
            
            double gamma = tau + 1;
            double beta  = std::sqrt(tau*(tau + 2))/gamma;
            double A     = (1 - beta)/beta;
            double Ap2   = A + 2;
            double B     = 0.5*beta*gamma*(gamma - 1)*(gamma - 2);
            double grej  = 2*(1 + A*B)/A;
            double z, g;
            do {
                td->fRndm->uniform_array(2, rndArray);
                double q = rndArray[0];
                z = 2*A*(2*q + Ap2*std::sqrt(q))/(Ap2*Ap2 - 4*q);
                g = (2 - z)*(1.0/(A + z) + B);
                
            } while(g < rndArray[1]*grej);
            
            cosTheta = 1 - z;
            sinTheta = std::sqrt(z*(2 - z));
        }
        return;
    }

    
//____________________
    double SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Alias(const MaterialCuts *matcut, double gammaEnIn, double r1, double r2, double r3){
        
        int mcindx       = matcut->GetIndex(); //GetIndex of the material cut
        int macindxlocal = fGlobalMatGCutIndxToLocal[mcindx]; //From global index to local index
        // the location of the first-gamma-energy lin-alias data of this mat-cut
        int indxstart    = macindxlocal*fNumSamplingGammaEnergies;//
        // determine gamma energy lower grid point
        double lgenergy  = std::log(gammaEnIn);//lower gamma energy grid point
        int genergyindx  = (int) ((lgenergy-fGammaEnLMin)*fGammaEnILDelta);
        
        if (genergyindx>=fNumSamplingGammaEnergies-1)
            genergyindx = fNumSamplingGammaEnergies-2;
        
        double plowgener = (fLSamplingGammaEnergies[genergyindx+1]-lgenergy)*fGammaEnILDelta;
        indxstart +=genergyindx;
        if (r1>plowgener|| !(fAliasData[indxstart]))
            ++indxstart;
        
        // sample the outgoing electron cosTheta
        double ecosTheta = fAliasSampler->SampleLinear(fAliasData[indxstart]->fXdata, fAliasData[indxstart]->fYdata,
                                                       fAliasData[indxstart]->fAliasW, fAliasData[indxstart]->fAliasIndx,
                                                       fNumSamplingElectronEnergies,r2,r3);
        return ecosTheta;
        
    }
    
    
    double SauterGavrilaPhotoElectricModel::CrossSectionPerVolume(const const MaterialCuts *matcut,
                                         LightTrack &track,
                                         double kineticEnergy,
                                         double cutEnergy,
                                         double maxEnergy)
    {
        //Dummy for now
        std::cout<<"DUMMY\n";
        return 0.0;
    }
    

    //NB: At the moment is only sampling the direction of the photoElectron using SauterGavrila distribution for the k-Shell differential cross section. All the rest is dummy (binding energies - selection of ionization shell - deexcitation)
    int SauterGavrilaPhotoElectricModel::SampleSecondaries(LightTrack &track,
                                                           std::vector<LightTrack> &sectracks,
                                                           Geant::GeantTaskData *td){
        int    numSecondaries      = 0;
        double gammaekin0          = track.GetKinE();
        const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
        const double *cuts         = matCut->GetProductionCutsInEnergy();
        double gammacut         = cuts[0];
        

    
        // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is above fHighEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is below gamma production cut and do nothing if yes;
        if (gammaekin0<GetLowEnergyUsageLimit() || gammaekin0>GetHighEnergyUsageLimit() || gammaekin0<=gammacut)
        {
            return numSecondaries;
        }
        // sample gamma energy
        // here we need 3 random number + 1 later for photon direction phi sampling
        double *rndArray = td->fDblArray;
        td->fRndm->uniform_array(4, rndArray);
        
    
        //commented out for now
        //SetCurrentCouple(couple);
        //const Material* aMaterial = couple->GetMaterial();
        
        //double energy = track.GetKinE();
        
        // select randomly one element constituing the material.
        //const Element* anElement = SelectRandomAtom(aMaterial,theGamma,energy);
        
        //
        // Photo electron
        //
        double cosTheta = track.GetDirZ();
        double sinTheta = 0.0;
        double phi      = 0.0;
        
        
        // Select atomic shell
        //dummy for now
        int nShells = 44;//anElement->GetNbOfAtomicShells();
        int i = 0;
        //for(; i<nShells; ++i) {
            
            //if(energy >= anElement->GetAtomicShell(i)) { break; }
        //}
        
        double edep = gammaekin0;
        
        // Normally one shell is available
        if (i < nShells)
        {
            
            //bindingEnergy dummy value for now
            double bindingEnergy = 0;//anElement->GetAtomicShell(i);
            edep = bindingEnergy;
            double esec = 0.0;
            
            // sample deexcitation
            //
            //MISSING for now
            
            // create photo electron
            //
            double elecKineEnergy = gammaekin0 - bindingEnergy;
            if (elecKineEnergy > gammacut) {
                
                // store original gamma directions in the lab frame
                double gamDirX0=track.GetDirX();
                double gamDirY0=track.GetDirY();
                double gamDirZ0=track.GetDirZ();
                
                double eDirX1;
                double eDirY1;
                double eDirZ1;
                
                double tau = gammaekin0/geant::kElectronMassC2;
                if (tau > 50.){
                    eDirX1=gamDirX0;
                    eDirY1=gamDirY0;
                    eDirZ1=gamDirZ0;
                }else
                {
                    //START REJECTION SAMPLING
                    //SamplePhotoElectronDirection_Rejection(gammaekin0, sinTheta, cosTheta, phi, td);
                    //END REJECTION SAMPLING
                    
                    //START ALIAS SAMPLING
                    cosTheta=SamplePhotoElectronDirection_Alias(matCut,gammaekin0,rndArray[0], rndArray[1], rndArray[2]);
                    phi = geant::kTwoPi * rndArray[3];
                    //END ALIAS SAMPLING
                    
                    sinTheta=std::sqrt((1 - cosTheta)*(1 + cosTheta));
                
                    // new photoelectron direction in the scattering frame
                    eDirX1  = sinTheta*std::cos(phi);
                    eDirY1  = sinTheta*std::sin(phi);
                    eDirZ1  = cosTheta;
                    
                    // rotate new photoelectron direction to the lab frame:
                    RotateToLabFrame(eDirX1, eDirY1, eDirZ1, gamDirX0, gamDirY0, gamDirZ0);
                }
                
        
                // create the secondary particle i.e. the photoelectron
                numSecondaries = 1;
                
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
        
                sectracks[secIndx].SetDirX(eDirX1);
                sectracks[secIndx].SetDirY(eDirY1);
                sectracks[secIndx].SetDirZ(eDirZ1);
                sectracks[secIndx].SetKinE(elecKineEnergy);
                sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // electron GV code
                sectracks[secIndx].SetMass(geant::kElectronMassC2);
                sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
            }
            else {
                edep += elecKineEnergy;
                elecKineEnergy = 0.0;
            }
            
            if(fabs(gammaekin0 - elecKineEnergy - esec - edep) > geant::eV) {
                std::cout << "### G4PEffectFluoModel dE(eV)= "
                << (gammaekin0 - elecKineEnergy - esec - edep)/geant::eV
                << " shell= " << i
                << "  E(keV)= " << gammaekin0/geant::keV
                << "  Ebind(keV)= " << bindingEnergy/geant::keV
                << "  Ee(keV)= " << elecKineEnergy/geant::keV
                << "  Esec(keV)= " << esec/geant::keV
                << "  Edep(keV)= " << edep/geant::keV
                << std::endl;
            }
        }
        
        //always kill primary photon
        track.SetTrackStatus(LTrackStatus::kKill);
        track.SetKinE(0.0);

        if(edep > 0.0) {
            //MISSING
            //fParticleChange->ProposeLocalEnergyDeposit(edep);
        }
        // return with number of secondaries i.e. 1 photoelectron
        return numSecondaries;
    }
    
    //static const int nlooplim = 1000;
    
    void SauterGavrilaPhotoElectricModel::InitSamplingTables() {
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
        
        /* UNCOMMENT TO DEBUG
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
    
    
    
    double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSection(int /*Zelement*/, double energy0, double cosTheta)
    {
       
        // based on Geant4 : G4SauterGavrilaAngularDistribution
        // input  : energy   (incomming photon energy)
        //          cosTheta (cons(theta) of photo-electron)
        // output : dsigmaK  (differential cross section, K-shell only)
        
        double tau = energy0 / geant::kElectronMassC2;
        
        double g = tau + 1.0;
        double invgamma = 1.0 / (tau + 1.0);
        double beta = std::sqrt(tau * (tau + 2.0)) * invgamma;
        
        double z = 1 - beta * cosTheta;
        double z2 = z * z;
        double z4 = z2 * z2;
        double y = 1 - cosTheta * cosTheta;
        
        double dsigma = (y / z4) * (1 + 0.5 * g * (g - 1) * (g - 2) * z);
        
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
    void SauterGavrilaPhotoElectricModel::BuildOneLinAlias(int ialias, double ecut){
        
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


                double cosMin=-1;
                double cosMax=1;
                double middlepoint=0.0;
                
                // fill the first, last and middle value (NB: we should convert to 0-1)
                fAliasData[ialias]->fXdata[0] =  cosMin;
                fAliasData[ialias]->fYdata[0] =  CalculateDiffCrossSection(1, gEnergy, cosMin);
                
                fAliasData[ialias]->fXdata[1] =  middlepoint;
                fAliasData[ialias]->fYdata[1] =  CalculateDiffCrossSection(1, gEnergy, middlepoint);
                
                fAliasData[ialias]->fXdata[2] =  cosMax;
                fAliasData[ialias]->fYdata[2] =  CalculateDiffCrossSection(1, gEnergy, cosMax);
                
                int numdata=3;
                
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
                        
                        double xsecVal=CalculateDiffCrossSection(1, gEnergy, xx);
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
    
}   // namespace geantphysics
