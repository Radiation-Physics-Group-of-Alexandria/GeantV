
#include "SauterGavrilaPhotoElectricModel.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"

//#include "Spline.h"
#include "GLIntegral.h"
#include "AliasTable.h"
#include "XSectionsVector.h"

#include "PhysicsParameters.h"

#include "Gamma.h"
#include "Electron.h"

#include "LightTrack.h"
#include "PhysicsData.h"


#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>



// from geantV
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"


using namespace std;
namespace geantphysics {
    
    std::vector<double>*  SauterGavrilaPhotoElectricModel::fParamHigh[] = {nullptr};
    std::vector<double>*  SauterGavrilaPhotoElectricModel::fParamLow[] = {nullptr};
    
    SauterGavrilaPhotoElectricModel::ShellData  ** SauterGavrilaPhotoElectricModel::fShellCrossSection = {nullptr};
    //SauterGavrilaPhotoElectricModel::CrossSectionsVector ** SauterGavrilaPhotoElectricModel::fLECSVector33 = {nullptr};
    //SauterGavrilaPhotoElectricModel::CrossSectionsVector ** SauterGavrilaPhotoElectricModel::fCSVector33 = {nullptr};
    
    int                   SauterGavrilaPhotoElectricModel::fNShells[] = {0};
    int                   SauterGavrilaPhotoElectricModel::fNShellsUsed[] = {0};
    
    Material*             SauterGavrilaPhotoElectricModel::fWater = nullptr;
    double                SauterGavrilaPhotoElectricModel::fWaterEnergyLimit = 0.0;
    
    
    
    
    SauterGavrilaPhotoElectricModel::SauterGavrilaPhotoElectricModel(const std::string &modelname)
    : EMModel(modelname){
        
        //all these values need to be seen and set with the optimal values
        fNumSamplingPrimEnergiesPerDecade = 20;    // should be set/get and must be done before init
        fNumSamplingPrimEnergies = 60;
        fNumSamplingAngles = 60;                        // at each energy grid points
        fMinPrimEnergy           =  0.1*geant::keV;     // minimum kinetic energy of the interacting gamma
        fMaxPrimEnergy           =  100*geant::GeV;     //10.0*geant::GeV; // maximum kinetic energy of the interacting gamma
        
        fPrimEnLMin                = 0.;       // will be set in InitSamplingTables if needed
        fPrimEnILDelta             = 0.;       // will be set in InitSamplingTables if needed
        fSamplingPrimEnergies    = nullptr;    // will be set in InitSamplingTables if needed
        fLSamplingPrimEnergies   = nullptr;    // will be set in InitSamplingTables if needed
        
        fAliasData                = nullptr;    // will be set in InitSamplingTables if needed
        fAliasSampler             = nullptr;
    
        fShellCrossSection        = nullptr;
        fCrossSection             = nullptr;
        fCrossSectionLE           = nullptr;
        //fCSVector33                = nullptr;
        //fLECSVector33               = nullptr;
        
    }
    
    SauterGavrilaPhotoElectricModel::~SauterGavrilaPhotoElectricModel() {
        
        //CLEANING fParamHigh
        for(int i=0; i<gMaxSizeData; ++i) {
            delete fParamHigh[i];
            fParamHigh[i] = 0;
            
        }
        //CLEANING fParamLow
        for(int i=0; i<gMaxSizeData; ++i) {
            delete fParamLow[i];
            fParamLow[i] = 0;
        }
        
        /*//CLEANING fLECSVector33
        if (fLECSVector33) {
            for(int i=0; i<gMaxSizeData; i++)
            {
                if(fLECSVector33[i])
                {
                    fLECSVector33[i]->fBinVector = std::vector<double>();
                    fLECSVector33[i]->fDataVector = std::vector<double>();
                    //delete fLECSVector[i]->sp;
                    delete fLECSVector33[i];
                }
            }
            
            delete [] fLECSVector33;
        }
        
        //CLEANING fCSVector33
        if (fCSVector33) {
            for(int i=0; i<gMaxSizeData; i++)
            {
                if(fCSVector33[i])
                {
                    fCSVector33[i]->fBinVector = std::vector<double>();
                    fCSVector33[i]->fDataVector = std::vector<double>();
                    //delete fCSVector[i]->sp;
                    delete fCSVector33[i];
                }
            }
            
            delete [] fCSVector33;
        }*/
        
        //CLEANING fShellCrossSection
        if (fShellCrossSection) {
            for(int i=0; i<gMaxSizeData; i++){
                if(fShellCrossSection[i])
                {
                    for(int j=0; j<fNShellsUsed[i]; j++)
                    {
                        fShellCrossSection[i]->fCompBinVector[j] = std::vector<double>();
                        fShellCrossSection[i]->fCompDataVector[j] = std::vector<double>();
                    }
                    
                    delete fShellCrossSection[i]->fCompID;
                    delete fShellCrossSection[i]->fCompLength;
                    delete fShellCrossSection[i];
                }
            }
            
            delete [] fShellCrossSection;
        }
        
        
        if (fSamplingPrimEnergies)
            delete [] fSamplingPrimEnergies;
        if (fLSamplingPrimEnergies)
            delete [] fLSamplingPrimEnergies;
        
        // NO ALIAS FOR THE MOMENT
         if (fAliasData) {
         for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
         if (fAliasData[i]) {
         delete [] fAliasData[i]->fXdata;
         delete [] fAliasData[i]->fYdata;
         delete [] fAliasData[i]->fAliasW;
         delete [] fAliasData[i]->fAliasIndx;
         delete fAliasData[i];
         }
         }
         delete [] fAliasData;
         }
         
         if (fAliasSampler)
         delete fAliasSampler;//*/
        
    }
    
    void SauterGavrilaPhotoElectricModel::Initialize() {
        EMModel::Initialize();
        fSecondaryInternalCode    = Electron::Definition()->GetInternalCode();
        InitializeModel();
    }
    
    void SauterGavrilaPhotoElectricModel::InitializeModel() {
        
        
        fMinPrimEnergy                    = GetLowEnergyUsageLimit();
        fMaxPrimEnergy                    = GetHighEnergyUsageLimit();
    
        
        
        //ALLOCATION fCrossSection
        if (fCrossSection) {
            delete [] fCrossSection;
            fCrossSection = nullptr;
        }
        
        fCrossSection = new bool [gMaxSizeData];
        for (int i=0; i<gMaxSizeData; ++i) {
            fCrossSection[i] = false;
        }
    
        
        //ALLOCATION fCrossSectionLE
        if (fCrossSectionLE) {
            delete [] fCrossSectionLE;
            fCrossSectionLE = nullptr;
            
        }
        
        fCrossSectionLE = new bool [gMaxSizeData];
        for (int i=0; i<gMaxSizeData; ++i) {
            fCrossSectionLE[i] = false;
        }
        
        if (fShellCrossSection) {
            for (int i=0; i<99 ; i++){
                
                delete fLE[i];
                fLE[i] = 0;
                delete fCS[i];
                fLE[i] = 0;
                
            }
        }
     
        
        
        //ALLOCATION fShellCrossSection
        if (fShellCrossSection) {
            for (int i=0; i<gMaxSizeData; ++i)
                for (int j=0; j<gNShellLimit; ++j){
                    if (fShellCrossSection[i]) {
                        fShellCrossSection[i]->fCompBinVector[j]=std::vector<double>();
                        fShellCrossSection[i]->fCompDataVector[j]=std::vector<double>();
                        delete fShellCrossSection[i]->fCompID;
                        delete fShellCrossSection[i]->fCompLength;
                        delete fShellCrossSection[i];
                    }
                }
            delete [] fShellCrossSection;
            fShellCrossSection = nullptr;
        }
        
        fShellCrossSection= new ShellData*[gMaxSizeData];
        
        for (int i=0; i<gMaxSizeData; ++i)
            for (int j=0; j<fNShellsUsed[i]; ++j){
                fShellCrossSection[i] = new ShellData;
                fShellCrossSection[i]->fCompBinVector[j].clear();
                fShellCrossSection[i]->fCompDataVector[j].clear();
                fShellCrossSection[i]->fCompID=nullptr;
                fShellCrossSection[i]->fCompLength=nullptr;
                
            }
        /*
        //ALLOCATION fLECSVector
        if (fLECSVector33) {
            for(int i=0; i<gMaxSizeData; i++)
            {
                if(fLECSVector33[i])
                {
                    fLECSVector33[i]->fBinVector = std::vector<double>();
                    fLECSVector33[i]->fDataVector = std::vector<double>();
                    //delete fLECSVector[i]->sp;
                    delete fLECSVector33[i];
                }
            }
            
            delete [] fLECSVector33;
            fLECSVector33 = nullptr;
            
        }
        
        fLECSVector33= new CrossSectionsVector*[gMaxSizeData];
        
        for (int i=0; i<gMaxSizeData; i++)
        {
            fLECSVector33[i]=new CrossSectionsVector;
            //fLECSVector[i]->sp=nullptr;
            fLECSVector33[i]->fBinVector.clear();
            fLECSVector33[i]->fDataVector.clear();
        }
        
        
        
        
        
        //ALLOCATION fCSVector
        if (fCSVector33) {
            for(int i=0; i<gMaxSizeData; i++)
            {
                if(fCSVector33[i])
                {
                    fCSVector33[i]->fBinVector = std::vector<double>();
                    fCSVector33[i]->fDataVector = std::vector<double>();
                    //delete fCSVector[i]->sp;
                    delete fCSVector33[i];
                }
            }
            
            delete [] fCSVector33;
            fCSVector33 = nullptr;
            
        }
        
        fCSVector33= new CrossSectionsVector*[gMaxSizeData];
        
        for (int i=0; i<gMaxSizeData; i++)
        {
            fCSVector33[i]=new CrossSectionsVector;
            //fCSVector[i]->sp=nullptr;
            fCSVector33[i]->fBinVector.clear();
            fCSVector33[i]->fDataVector.clear();
        }
         */
        
        fVerboseLevel=1;
        LoadData();
        //NOT NEEDED FOR THE MOMENT - NO ALIAS USED
        //InitSamplingTables();
        
    }
    
    void SauterGavrilaPhotoElectricModel::LoadData()
    {
        
        int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
        // get list of active region
        std::vector<bool> isActiveInRegion = GetListActiveRegions();
        for (int i=0; i<numMatCuts; ++i) {
            const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
            // if this MaterialCuts belongs to a region where this model is active:
            if (isActiveInRegion[matCut->GetRegionIndex()]) {
                // get the list of elements
                const Vector_t<Element*> theElements =  matCut->GetMaterial()->GetElementVector();
                int numElems = theElements.size();
                for (int j=0; j<numElems; ++j) {
                    // if RatinAliasDataPerElement has not been built for this element yet => build
                    double zet = theElements[j]->GetZ();
                    int elementIndx = std::lrint(zet);
                    //std::cout<<"\n***Calling ReadData on: "<<elementIndx<<" ....\n";
                    ReadData(elementIndx);
                }
            }
        }
    }
    
    void SauterGavrilaPhotoElectricModel::ReadData(int Z)
    {
        //bool debug= false;
        using geant::MeV;
        using geant::barn;
        
        if (fVerboseLevel > 1)
        {
            std::cout << "Calling ReadData() of SauterGavrilaPhotoElectricModel"
            << std::endl;
        }
        //to change
        //if( ( (fCrossSection[Z] && Z<23) || (!fCrossSection[Z] && Z>22 )) && ( (fCrossSectionLE[Z] && Z>2) || (!fCrossSectionLE[Z] && Z<3)) ) {std::cout<<"Data loaded before!\n"; return;}
        //else std::cout<<"Data not loaded before!\n";
        
        if(  (fCrossSection[Z]) && ( (fCrossSectionLE[Z] && Z>2) || (!fCrossSectionLE[Z] && Z<3)) )
        {//std::cout<<"Data loaded before!\n";
            return;}
        //else std::cout<<"Data not loaded before!\n";
        
        
        
        // get the path to the main physics data directory
        char *path = std::getenv("GEANT_PHYSICS_DATA");
        if (!path) {
            std::cerr<<"******   ERROR in SauterGavrilaPhotoElectricModel::ReadData() \n"
            <<"         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
            <<"         environment variable to the location of Geant data directory!\n"
            <<std::endl;
            exit(1);
        }
        
        // spline for photoeffect total x-section above K-shell - but below the parameterized ones
        
       
        
        fCrossSection[Z] =true;
        //fCrossSection[Z]->SetSpline(true); //IMPORTANT
        
        std::ostringstream ost;
        ost << path << "/livermore/phot_epics2014/pe-cs-" << Z <<".dat";
        std::ifstream fin(ost.str().c_str());
        if( !fin.is_open()){
            std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str()
            << "> is not opened!" << std::endl;
            
            return;
        } else {
            if(fVerboseLevel > 3) { std::cout << "File " << ost.str().c_str()
                << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;}
            
            
            fCSVector[Z]= new XSectionsVector;
            fin >> fCSVector[Z]->edgeMin >> fCSVector[Z]->edgeMax >> fCSVector[Z]->numberOfNodes;
            
            int siz=0;
            fin >> siz;
            if (fin.fail() || siz<=0)
            {
                std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str()
                << "> is not opened!" << std::endl;
            }
            
            fCSVector[Z]->fBinVector.reserve(siz);
            fCSVector[Z]->fDataVector.reserve(siz);
            
            
            double vBin, vData;
            
            for(int i = 0; i < siz ; i++)
            {
                vBin = 0.;
                vData= 0.;
                fin >> vBin >> vData;
                
                if (fin.fail())  {
                    std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str()
                    << "> is not opened!" << std::endl;
                }
                
                fCSVector[Z]->fBinVector.push_back(vBin*MeV);
                fCSVector[Z]->fDataVector.push_back(vData*barn);
            }
            
            // to remove any inconsistency
            fCSVector[Z]->numberOfNodes = siz;
            fCSVector[Z]->edgeMin = fCSVector[Z]->fBinVector[0];
            fCSVector[Z]->edgeMax = fCSVector[Z]->fBinVector[fCSVector[Z]->numberOfNodes-1];
            
            
            //to do: add
            //fCSVector[Z]->sp= new Spline(&(fCSVector[Z]->fBinVector[0]),&(fCSVector[Z]->fDataVector[0]),fCSVector[Z]->numberOfNodes);
            fin.close();
        }
        
        //}
        
        
        // read high-energy fit parameters
        fParamHigh[Z] = new std::vector<double>;
        
        int n1 = 0;
        int n2 = 0;
        double x;
        std::ostringstream ost1;
        ost1 << path << "/livermore/phot_epics2014/pe-high-" << Z <<".dat";
        std::ifstream fin1(ost1.str().c_str());
        if( !fin1.is_open()) {
            std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost1.str().c_str()
            << "> is not opened!" << std::endl;
            return;
        } else {
            if(fVerboseLevel > 3) {
                std::cout << "File " << ost1.str().c_str()
                << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
            }
            fin1 >> n1;
            if(fin1.fail()) { return; }
            if(0 > n1 || n1 >= INT_MAX) { n1 = 0; }
            
            fin1 >> n2;
            if(fin1.fail()) { return; }
            if(0 > n2 || n2 >= INT_MAX) { n2 = 0; }
            
            fin1 >> x;
            if(fin1.fail()) { return; }
            
            fNShells[Z] = n1;
            fParamHigh[Z]->reserve(7*n1+1);
            fParamHigh[Z]->push_back(x*MeV);
            for(int i=0; i<n1; ++i) {
                for(int j=0; j<7; ++j) {
                    fin1 >> x;
                    if(0 == j) { x *= MeV; }
                    else       { x *= barn; }
                    fParamHigh[Z]->push_back(x);
                }
            }
            fin1.close();
        }
        
        //std::cout<<"pe-high parameterization for ["<<Z<<"], loaded\n";
        
        // read low-energy fit parameters
        fParamLow[Z] = new std::vector<double>;
        int n1_low = 0;
        int n2_low = 0;
        double x_low;
        std::ostringstream ost1_low;
        ost1_low << path << "/livermore/phot_epics2014/pe-low-" << Z <<".dat";
        std::ifstream fin1_low(ost1_low.str().c_str());
        if( !fin1_low.is_open()) {
            std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost1_low.str().c_str()
            << "> is not opened!" << std::endl;
            return;
        } else {
            if(fVerboseLevel > 3) {
                std::cout << "File " << ost1_low.str().c_str()
                << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
            }
            fin1_low >> n1_low;
            if(fin1_low.fail()) { return; }
            if(0 > n1_low || n1_low >= INT_MAX) { n1_low = 0; }
            
            fin1_low >> n2_low;
            if(fin1_low.fail()) { return; }
            if(0 > n2_low || n2_low >= INT_MAX) { n2_low = 0; }
            
            fin1_low >> x_low;
            if(fin1_low.fail()) { return; }
            
            fNShells[Z] = n1_low;
            fParamLow[Z]->reserve(7*n1_low+1);
            fParamLow[Z]->push_back(x_low*MeV);
            for(int i=0; i<n1_low; ++i) {
                for(int j=0; j<7; ++j) {
                    fin1_low >> x_low;
                    if(0 == j) { x_low *= MeV; }
                    else       { x_low *= barn; }
                    fParamLow[Z]->push_back(x_low);
                }
            }
            fin1_low.close();
        }
        
        //std::cout<<"pe-low parameterization for ["<<Z<<"], loaded\n";
        
        
        // there is a possibility to use only main shells
        if(gNShellLimit < n2)
        {
            n2 = gNShellLimit;
        }
        fNShellsUsed[Z]=n2;
        
        fShellCrossSection[Z] = new ShellData;
        fShellCrossSection[Z]->fCompDataVector= new std::vector<double>[n2];
        fShellCrossSection[Z]->fCompBinVector= new std::vector<double>[n2];
        fShellCrossSection[Z]->fCompID=new int[n2];
        fShellCrossSection[Z]->fCompLength = new size_t [n2];
        
        for (int i=0 ; i<n2 ; i++)
        {
            fShellCrossSection[Z]->fCompDataVector[i].clear();
            fShellCrossSection[Z]->fCompBinVector[i].clear();
            fShellCrossSection[Z]->fCompLength[i] = 0;
        }
        
        fNShellsUsed[Z] = n2;
        
        //If more than one shell is used -> Read sub-shells cross section data
        if(1 < n2) {
            std::ostringstream ost2;
            ost2 << path << "/livermore/phot_epics2014/pe-ss-cs-" << Z <<".dat";
            std::ifstream fin2(ost2.str().c_str());
            if( !fin2.is_open()) {
                std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost2.str().c_str()
                << "> is not opened!" << std::endl;
                return;
            } else {
                if(fVerboseLevel > 3) {
                    std::cout << "File " << ost2.str().c_str()
                    << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
                }
                
                int n3, n4;
                double y;
                
                for(int i=0; i<n2; ++i)
                {
                    fin2 >> x >> y >> n3 >> n4;
                    fShellCrossSection[Z]->fCompBinVector[i].clear();
                    fShellCrossSection[Z]->fCompDataVector[i].clear();
                    fShellCrossSection[Z]->fCompBinVector[i].reserve(n3);
                    fShellCrossSection[Z]->fCompDataVector[i].reserve(n3);
                    
                    for(int j=0; j<n3; ++j)
                    {
                        fin2 >> x >> y;
                        fShellCrossSection[Z]->fCompBinVector[i].push_back(x*MeV);
                        fShellCrossSection[Z]->fCompDataVector[i].push_back(y*barn);
                        ++fShellCrossSection[Z]->fCompLength[i];
                    }
                    
                    fShellCrossSection[Z]->fCompID[i]=n4;
                    
                }
                
                fin2.close();
            }
        }
        //std::cout<<"pe-ss-cs- cross sections for ["<<Z<<"], loaded\n";
        
        // no spline for photoeffect total x-section below K-shell
        if(1 < fNShells[Z]) {
            fCrossSectionLE[Z] = true;
            std::ostringstream ost3;
            ost3 << path << "/livermore/phot_epics2014/pe-le-cs-" << Z <<".dat";
            std::ifstream fin3(ost3.str().c_str());
            if( !fin3.is_open()) {
                std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str()
                << "> is not opened!" << std::endl;
                return;
            } else {
                if(fVerboseLevel > 3) {
                    std::cout << "File " << ost3.str().c_str()
                    << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
                }
    
                // binning
                fLECSVector[Z]= new XSectionsVector;
                //fLECSVector[Z]->fBinVector.clear();
                //fLECSVector[Z]->fDataVector.clear();
                
                fin3 >> fLECSVector[Z]->edgeMin >> fLECSVector[Z]->edgeMax >> fLECSVector[Z]->numberOfNodes;
                int siz=0;
                fin3 >> siz;
                if (fin3.fail() || siz<=0)
                {
                    std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str()
                    << "> is not opened!" << std::endl;
                }
                
                fLECSVector[Z]->fBinVector.reserve(siz);
                fLECSVector[Z]->fDataVector.reserve(siz);
                
                double vBin, vData;
                
                for(int i = 0; i < siz ; i++)
                {
                    vBin = 0.;
                    vData= 0.;
                    fin3 >> vBin >> vData;
                    
                    //Scale vector
                    vBin  *= MeV;
                    vData *= barn;
                    
                    if (fin3.fail())  {
                        std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str()
                        << "> is not opened!" << std::endl;
                    }
                    
                    fLECSVector[Z]->fBinVector.push_back(vBin);
                    fLECSVector[Z]->fDataVector.push_back(vData);
                    
                }
                
                // to remove any inconsistency
                fLECSVector[Z]->numberOfNodes = siz;
                fLECSVector[Z]->edgeMin = fLECSVector[Z]->fBinVector[0];
                fLECSVector[Z]->edgeMax = fLECSVector[Z]->fBinVector[fLECSVector[Z]->numberOfNodes-1];
                
                //Here we will use LINEAR interpolation
                fin3.close();
            }
            //std::cout<<"pe-le-cs- cross sections for ["<<Z<<"], loaded\n";
        }
    }
    
    
    //____________________
    //NB: cosTheta is supposed to contain the dirZ of the incoming photon
    void SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Rejection(double gammaEnIn, double &sinTheta, double &cosTheta, double &phi, Geant::GeantTaskData *td){
        
        
        double *rndArray = td->fDblArray; //needed?
        td->fRndm->uniform_array(1, rndArray);
        phi     = geant::kTwoPi * rndArray[0];
        
        if (gammaEnIn > 1*geant::GeV) {
            
            sinTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
            
        } else
        {
            
            //1) initialize energy-dependent variables
            // Variable naming according to Eq. (2.24) of Penelope Manual
            // (pag. 44)
            double gamma = 1.0 + gammaEnIn/geant::kElectronMassC2;
            double gamma2 = gamma*gamma;
            double beta = std::sqrt((gamma2-1.0)/gamma2);
            
            // ac corresponds to "A" of Eq. (2.31)
            //
            double ac = (1.0/beta) - 1.0;
            double a1 = 0.5*beta*gamma*(gamma-1.0)*(gamma-2.0);
            double a2 = ac + 2.0;
            // gtmax = maximum of the rejection function according to Eq. (2.28), obtained for tsam=0
            double gtmax = 2.0*(a1 + 1.0/ac);
            
            double tsam = 0;
            double gtr = 0;
            
            //2) sampling. Eq. (2.31) of Penelope Manual
            // tsam = 1-std::cos(theta)
            // gtr = rejection function according to Eq. (2.28)
            do{
                //double rand = G4UniformRand();
                td->fRndm->uniform_array(2, rndArray);
                tsam = 2.0*ac * (2.0*rndArray[0] + a2*std::sqrt(rndArray[0])) / (a2*a2 - 4.0*rndArray[0]);
                gtr = (2.0 - tsam) * (a1 + 1.0/(ac+tsam));
            }while(rndArray[1]*gtmax > gtr);
            
            cosTheta = 1.0-tsam;
            sinTheta = std::sqrt(1-cosTheta*cosTheta);
            
        }
        return;
    }
    
    
    //____________________
    double SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Alias(double primekin, double r1, double r2, double r3){
        
        
        // determine primary energy lower grid point
        if (primekin > 1*geant::GeV) {
            
            return 1.;
        } else
        {
            
            double lGammaEnergy  = std::log(primekin);
            int gammaEnergyIndx  = (int) ((lGammaEnergy-fPrimEnLMin)*fPrimEnILDelta);
            //
            if (gammaEnergyIndx>=fNumSamplingPrimEnergies-1)
                gammaEnergyIndx = fNumSamplingPrimEnergies-2;
            //
            double pLowerGammaEner = (fLSamplingPrimEnergies[gammaEnergyIndx+1]-lGammaEnergy)*fPrimEnILDelta;
            if (r1>pLowerGammaEner) {
                ++gammaEnergyIndx;
            }
            // sample the outgoing electron cosTheta
            double ecosTheta = fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx]->fXdata, fAliasData[gammaEnergyIndx]->fYdata,fAliasData[gammaEnergyIndx]->fAliasW, fAliasData[gammaEnergyIndx]->fAliasIndx,fAliasData[gammaEnergyIndx]->fNumdata,r2,r3);
            
            //This have to be seen (Transformation)
            //double xsi=fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx]->fXdata, fAliasData[gammaEnergyIndx]->fYdata,fAliasData[gammaEnergyIndx]->fAliasW, fAliasData[gammaEnergyIndx]->fAliasIndx,fAliasData[gammaEnergyIndx]->fNumdata,r2,r3);
            //double ecosTheta= 1-std::exp(xsi);
            
            return ecosTheta;
            
        }
        
    }
    
    double SauterGavrilaPhotoElectricModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy,
                                                                   const Particle*) {
        
        double xsec  = 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        // compute the parameterized atomic cross section: depends only on target Z and gamma energy.
        xsec = ComputeXSectionPerAtom(elem->GetZ(), kinenergy);
        return xsec;
    }
    
    
    double SauterGavrilaPhotoElectricModel::ComputeXSectionPerAtom(double zeta, double energy)
    {
        
        int verboseLevel = 0;
        using geant::keV;
        using geant::MeV;
        using geant::barn;
        
        
        if (verboseLevel > 3) {
            std::cout << "G4LivermorePhotoElectricModel_new::ComputeCrossSectionPerAtom():"
            << " Z= " << zeta << "  R(keV)= " << energy/keV << std::endl;
        }
        double cs = 0.0;
        int Z = std::lrint(zeta);
        if(Z < 1 || Z >= gMaxSizeData) { return cs; }
        
        
        //EXTRA CHECK --
        //if( !fCrossSection[Z] && !fCrossSectionLE[Z] )
        //{
        //    ReadData(Z);
        //    if(!fCrossSectionLE[Z] && !fCrossSection[Z]) { return cs; }
        //}
        
        int idx = fNShells[Z]*7 - 5; //7: rows in the parameterization file; 5: number of parameters
        
        if (energy < (*(fParamHigh[Z]))[idx-1]) { energy = (*(fParamHigh[Z]))[idx-1];}
        
        double x1 = (MeV)/(energy);
        double x2 = x1*x1;
        double x3 = x2*x1;
        
        // high energy parameterisation
        if(energy >= (*(fParamHigh[Z]))[0]) {
            
            double x4 = x2*x2;
            double x5 = x4*x1;
            
            cs = x1*((*(fParamHigh[Z]))[idx] + x1*(*(fParamHigh[Z]))[idx+1]
                     + x2*(*(fParamHigh[Z]))[idx+2] + x3*(*(fParamHigh[Z]))[idx+3]
                     + x4*(*(fParamHigh[Z]))[idx+4]+ x5*(*(fParamHigh[Z]))[idx+5]);
            
            
        }
        // low energy parameterisation
        else if(energy >= (*(fParamLow[Z]))[0]) {
            
            double x4 = x2*x2;
            double x5 = x4*x1;//this variable usage can probably be optimized
            cs = x1*((*(fParamLow[Z]))[idx] + x1*(*(fParamLow[Z]))[idx+1]
                     + x2*(*(fParamLow[Z]))[idx+2] + x3*(*(fParamLow[Z]))[idx+3]
                     + x4*(*(fParamLow[Z]))[idx+4]+ x5*(*(fParamLow[Z]))[idx+5]);
        }
        
        // Tabulated values above k-shell ionization energy
        else if(energy >= (*(fParamHigh[Z]))[1]) {
            
            //size_t index=0;
            //double value;
            //size_t index= GetIndex(energy,fCSVector, Z);
            ///THIS MUST BE SUBSTITUTED WITH SPLINE INTERPOLATOR
            //double value = LinearInterpolation(energy, fCSVector[Z]->fBinVector, fCSVector[Z]->fDataVector,  index);
            size_t index=0;
            double value=fCSVector[Z]->GetValue(energy, index);
            
            //double value=GetValuefCSVector(energy,Z);
            cs=x3*value;
            
        }
        // Tabulated values below k-shell ionization energy
        else
        {
            //index=0;
            //double value;
            //size_t index= GetIndex(energy,fLECSVector, Z);
            //double value = LinearInterpolation(energy, fLECSVector[Z]->fBinVector, fLECSVector[Z]->fDataVector,  index);
            //double value=GetValuefLECSVector(energy,Z);
            size_t index=0;
            double value=fLECSVector[Z]->GetValue(energy,index);
            
            cs=x3*value;
            
        }
        if (verboseLevel > 1) {
            std::cout << "LivermorePhotoElectricModel: E(keV)= " << energy/keV
            << " Z= " << Z << " cross(barn)= " << cs/barn << std::endl;
        }
        
        return cs;
        
    }
    
    double SauterGavrilaPhotoElectricModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                                       const Particle* /*particle*/) {
        
        double xsec = 0.0;
        //double xsecTemp= 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        // compute the macroscopic cross section as the sum of the atomic cross sections weighted by the number of atoms in
        // in unit volume.
        const Material *mat =  matcut->GetMaterial();
        double       egamma = kinenergy;
        // we need the element composition of this material
        const Vector_t<Element*> theElements    = mat->GetElementVector();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        int   numElems = theElements.size();
        for (int iel=0; iel<numElems; ++iel) {
            //xsecTemp=ComputeXSectionPerAtom(theElements[iel]->GetZ(), egamma);
            xsec += theAtomicNumDensityVector[iel]*ComputeXSectionPerAtom(theElements[iel]->GetZ(), egamma);
        }
        return xsec;
    }
    
    
    int SauterGavrilaPhotoElectricModel::SampleTargetElementIndex (const MaterialCuts *matCut, double gammaekin0, Geant::GeantTaskData *td)
    {
        int index=0;
        //retrieve the elements vector
        const Vector_t<Element*> theElements = matCut->GetMaterial()->GetElementVector();
        //retrieve the number of elements in the material
        int num    = matCut->GetMaterial()->GetNumberOfElements();
        if (num > 1)
        {
            double macxsec=ComputeMacroscopicXSection(matCut,gammaekin0, Gamma::Definition());
            double rnd=macxsec * td->fRndm->uniform();
            
            const Material *mat =  matCut->GetMaterial();
            const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
            double cumxsec=0.;
            for(; index<num-1; ++index)
            {
                double xsec = theAtomicNumDensityVector[index]* ComputeXSectionPerAtom(theElements[index]->GetZ(), gammaekin0);
                cumxsec += xsec;
                if (rnd <= cumxsec)
                {
                    break;
                }
            }
        }
        return index;
    }
    
    void SauterGavrilaPhotoElectricModel::TestSampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td){
        
        std::cout<<"testSampleTargetElementIndex\n";
        int index=0;
        double sum=0;
        //retrieve the elements vector
        const Vector_t<Element*> theElements = matcut->GetMaterial()->GetElementVector();
        //retrieve the number of elements in the material
        int num    = matcut->GetMaterial()->GetNumberOfElements();
        double xsec[num];
        double xsecSampled[num];
        
        const Material *mat =  matcut->GetMaterial();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        for (int i=0; i<num; i++)
        {
            xsec[i] = theAtomicNumDensityVector[i]* ComputeXSectionPerAtom(theElements[i]->GetZ(), energy);
            xsecSampled[i]=0.;
            sum+=xsec[i];
        }
        for (int i=0; i<1000000000; i++){
            index= SampleTargetElementIndex(matcut, energy, td);
            xsecSampled[index]++;
        }
        for (int i=0; i<num; i++)
        {
            xsec[i]/=sum;
            xsecSampled[i]/=1000000000;
        }
        
        char filename[521];
        sprintf(filename,"SampleTargetElementIndexTest_%s", (matcut->GetMaterial()->GetName()).c_str() );
        FILE *f     = fopen(filename,"w");
        
        for (int i=0; i<num; ++i) {
            fprintf(f,"%d\t%.8g\t%.8g\n",i,xsec[i],xsecSampled[i]);
        }
        fclose(f);
        
    }
    
    
    int SauterGavrilaPhotoElectricModel::SampleSecondaries(LightTrack &track,
                                                           Geant::GeantTaskData *td){
        
        
        using geant::MeV;
        //int    numSecondaries      = 0;
        double gammaekin0          = track.GetKinE();
        
        // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is above fHighEnergyUsageLimit and do nothing if yes;
        
        if (gammaekin0<GetLowEnergyUsageLimit() || gammaekin0>GetHighEnergyUsageLimit())
        {
            return 0; //numSecondaries
        }
        
        //interaction is possible so sample target element
        MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[track.GetMaterialCutCoupleIndex()];
        const Vector_t<Element*> theElements = matCut->GetMaterial()->GetElementVector();
        
        int targetElemIndx = 0;
        if (theElements.size()>1) {
            //uncomment the following line to test SampleTargetElementIndex
            //testSampleTargetElementIndex(matCut, gammaekin0, td );
            targetElemIndx = SampleTargetElementIndex(matCut, gammaekin0, td);
        }
        double  zeta  = theElements[targetElemIndx]->GetZ();
        int     Z = std::lrint(zeta);
    
        // if element was not initialised, gamma should be absorbed
        if(!fCrossSectionLE[Z] && !fCrossSection[Z]) {
            track.SetEnergyDeposit(gammaekin0);
            //std::cout<<"Model not initialized, Exiting!\n";
            return 0;
        }
    
        //SAMPLING OF THE SHELL
        size_t shellIdx = 0;
        size_t nn = fNShellsUsed[Z];
        bool barbra= false;
        if(nn > 1)
        {
            // sample gamma energy
            double *rndArray = td->fDblArray;
            td->fRndm->uniform_array(1, rndArray);
            
            if(gammaekin0 >= (*(fParamHigh[Z]))[0])
            {
                
                double x1 = (MeV)/gammaekin0;
                double x2 = x1*x1;
                double x3 = x2*x1;
                double x4 = x3*x1;
                double x5 = x4*x1;
                int idx   = nn*7 - 5;
                // when do sampling common factors are not taken into account
                // so cross section is not real
                double cs0 = rndArray[0]*(     (*(fParamHigh[Z]))[idx]
                                          + x1*(*(fParamHigh[Z]))[idx+1]
                                          + x2*(*(fParamHigh[Z]))[idx+2]
                                          + x3*(*(fParamHigh[Z]))[idx+3]
                                          + x4*(*(fParamHigh[Z]))[idx+4]
                                          + x5*(*(fParamHigh[Z]))[idx+5]);
                
                for(shellIdx=0; shellIdx<nn; ++shellIdx)
                {
                    idx = shellIdx*7 + 2;
                    if(gammaekin0 > (*(fParamHigh[Z]))[idx-1])
                    {
                        double cs =
                        (*(fParamHigh[Z]))[idx]
                        + x1*(*(fParamHigh[Z]))[idx+1]
                        + x2*(*(fParamHigh[Z]))[idx+2]
                        + x3*(*(fParamHigh[Z]))[idx+3]
                        + x4*(*(fParamHigh[Z]))[idx+4]
                        + x5*(*(fParamHigh[Z]))[idx+5];
                        if(cs >= cs0) {break;}
                    }
                }
                if(shellIdx >= nn) { shellIdx = nn-1; }
                
            }
            else if(gammaekin0 >= (*(fParamLow[Z]))[0])
            {
                double x1 = MeV/gammaekin0;//Precision lost?
                double x2 = x1*x1;
                double x3 = x2*x1;
                double x4 = x3*x1;
                double x5 = x4*x1;
                int idx   = nn*7 - 5;
                // when do sampling common factors are not taken into account
                // so cross section is not real
                double cs0 = rndArray[0]*((*(fParamLow[Z]))[idx]
                                          + x1*(*(fParamLow[Z]))[idx+1]
                                          + x2*(*(fParamLow[Z]))[idx+2]
                                          + x3*(*(fParamLow[Z]))[idx+3]
                                          + x4*(*(fParamLow[Z]))[idx+4]
                                          + x5*(*(fParamLow[Z]))[idx+5]);
                for(shellIdx=0; shellIdx<nn; ++shellIdx)
                {
                    idx = shellIdx*7 + 2;
                    if(gammaekin0 > (*(fParamHigh[Z]))[idx-1])
                    {
                        double cs = (*(fParamLow[Z]))[idx] + x1*(*(fParamLow[Z]))[idx+1]
                        + x2*(*(fParamLow[Z]))[idx+2] + x3*(*(fParamLow[Z]))[idx+3]
                        + x4*(*(fParamLow[Z]))[idx+4]+ x5*(*(fParamLow[Z]))[idx+5];
                        if(cs >= cs0) { break; }
                    }
                }
                if(shellIdx >= nn) {shellIdx = nn-1;}
            }
            else
            {
                
                // when do sampling common factors are not taken into account
                // so cross section is not real
                double cs = rndArray[0];
                
                if(gammaekin0 >= (*(fParamHigh[Z]))[1]) {
                    //above K-shell binding energy
                    //std::cout<<"*********** ****** *** *** ***\n";
                    //size_t index= GetIndex(gammaekin0,fCSVector, Z);
                    //THIS MUST BE SUBSTITUTED WITH THE SPLINE INTERPOLATION
                    //std::cout<<"Index: "<<index<<" - gamma energy: "<<gammaekin0<<std::endl;
                    //double value = LinearInterpolation(gammaekin0, fCSVector[Z]->fBinVector, fCSVector[Z]->fDataVector,  index);
                    //double value=GetValuefCSVector(gammaekin0,Z);
                    size_t index=0;
                    double value=fCSVector[Z]->GetValue(gammaekin0, index);
                    cs*=value;
                    //std::cout<<"Intepolated value: "<< value <<" - new cs: " <<cs<<std::endl;
                    
                }
                
                else
                {
         
                    //below K-shell binding energy
                    size_t index=0;
                    if(barbra)std::cout<<"START \n 4. below K-shell binding energy\n";
                    //double value= GetValuefLECSVector(gammaekin0,Z);
                    
                    double value=fLECSVector[Z]->GetValue(gammaekin0, index);
                    //size_t index= GetIndex(gammaekin0,fLECSVector, Z);
                    //if(barbra) std::cout<<"Index: "<<index<<" - gamma energy: "<<gammaekin0<<std::endl;
                    //double value = LinearInterpolation(gammaekin0, fLECSVector[Z]->fBinVector, fLECSVector[Z]->fDataVector,  index);
                    cs*=value;
                    if(barbra) std::cout<<"Calculated CS is "<<cs<<" because value was: "<<value<<std::endl;
                }
                size_t j=0;
                //clock_t start_time = clock();
                
                for(j=0; j<nn; ++j)
                {
                    
                    shellIdx = (size_t)fShellCrossSection[Z]->fCompID[j];
                    if(barbra)std::cout<<"\nshellIdx: "<<shellIdx<<std::endl;
                    if(barbra)std::cout<<"gammaekin0: "<<gammaekin0<<" is > than "<<(*(fParamLow[Z]))[7*shellIdx+1]<<"?"<<std::endl;
                    if(gammaekin0 > (*(fParamLow[Z]))[7*shellIdx+1]) {
                        //size_t bin= FindCSBinLocation(gammaekin0,fShellCrossSection[Z]->fCompLength[shellIdx], fShellCrossSection[Z]->fCompBinVector[shellIdx]);
                        //double value = LinearInterpolation (gammaekin0, fShellCrossSection[Z]->fCompBinVector[shellIdx], fShellCrossSection[Z]->fCompDataVector[shellIdx],  bin);
                        cs-=GetValue(gammaekin0, Z, shellIdx);
                        //if(barbra)std::cout<<"Yes, and cs now: cs - "<<value<<" = "<<cs<<"\n";
                    }
                    if(cs <= 0.0 || j+1 == nn) {if(barbra) std::cout<<"I am breaking because: "<<cs<< "-- "<< j+1<< " --- "<<nn<<std::endl;
                        break;}
                }
            }
        }
        
        //std::cout<<"END OF SAMPLING\n **********\n";
        // END: SAMPLING OF THE SHELL
        
        //Retrieving ionized shell bindingEnergy
        double bindingEnergy = (*(fParamHigh[Z]))[shellIdx*7 + 1];
        
        if(gammaekin0 < bindingEnergy) {
            track.SetEnergyDeposit(gammaekin0);
            return 0; //numSecondaries
        }
        //since edep is equal to bindingenergy I get rid of it
        //double edep = bindingEnergy;
     
        //  deexcitation is MISSING for now
        /*
         DEEXCITATION
         */
        //
        
        
        // Create the secondary particle: the photoelectron
        double elecKineEnergy = gammaekin0 - bindingEnergy;
        double cosTheta = track.GetDirZ();
        double sinTheta = 0.0;
        double phi      = 0.0;
        
        //*************START REJECTION SAMPLING
        //
        SamplePhotoElectronDirection_Rejection(gammaekin0, sinTheta, cosTheta, phi, td);
        //
        //*************END REJECTION SAMPLING
        
        
        //*************START ALIAS SAMPLING
        //
        //td->fRndm->uniform_array(4, rndArray);
        //double *rndArray2 = td->fDblArray;
        //td->fRndm->uniform_array(4, rndArray2);
        //cosTheta=SamplePhotoElectronDirection_Alias(gammaekin0, rndArray2[0], rndArray2[1], rndArray2[2]);
        //phi = geant::kTwoPi * rndArray2[3];
        //sinTheta=std::sqrt((1 - cosTheta)*(1 + cosTheta));
        //*************END ALIAS SAMPLING
        
        
        // new photoelectron direction in the scattering frame
        double eDirX1  = sinTheta*std::cos(phi);
        double eDirY1  = sinTheta*std::sin(phi);
        double eDirZ1  = cosTheta;
        
        
        // store original gamma directions in the lab frame
        //double gamDirX0=track.GetDirX();
        //double gamDirY0=track.GetDirY();
        //double gamDirZ0=track.GetDirZ();
        //RotateToLabFrame(eDirX1, eDirY1, eDirZ1, track.GetDirX(), gamDirY0, gamDirZ0);
        
        // rotate new photoelectron direction to the lab frame:
        RotateToLabFrame(eDirX1, eDirY1, eDirZ1,track.GetDirX(), track.GetDirY(), track.GetDirZ());
        
        
        // create the secondary particle i.e. the photoelectron
        //numSecondaries = 1;
        
        // current capacity of secondary track container
        int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
        // currently used secondary tracks in the container
        int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();
        if (curSizeOfSecList-curNumUsedSecs<1) {
            td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
        }
        int secIndx = curNumUsedSecs;
        curNumUsedSecs +=1;
        td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
        
        std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
        
        // this is known since it is a secondary track
        //  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
        
        sectracks[secIndx].SetDirX(eDirX1);
        sectracks[secIndx].SetDirY(eDirY1);
        sectracks[secIndx].SetDirZ(eDirZ1);
        sectracks[secIndx].SetKinE(elecKineEnergy);
        //std::cout<<"sectracks: elecKineEnergy: "<< elecKineEnergy<<"\n";
        sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // electron GV code
        sectracks[secIndx].SetMass(geant::kElectronMassC2);
        sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
        
        
        /*if(fabs(gammaekin0 - elecKineEnergy - esec - edep) > geant::eV) {
         std::cout << "### SauterGavrilaPhotoElectricModel dE(eV)= "
         << (gammaekin0 - elecKineEnergy - esec - edep)/geant::eV
         << "  shell= " << shellIdx
         << "  E(keV)= " << gammaekin0/geant::keV
         << "  Ebind(keV)= " << bindingEnergy/geant::keV
         << "  Ee(keV)= " << elecKineEnergy/geant::keV
         << "  Esec(keV)= " << esec/geant::keV
         << "  Edep(keV)= " << edep/geant::keV
         << std::endl;
         }*/
        
        
        //always kill primary photon
        track.SetTrackStatus(LTrackStatus::kKill);
        track.SetKinE(0.0);
        //edep is = bindingEnergy
        //if(edep > 0.0) {
        if(bindingEnergy > 0.0) {
            track.SetEnergyDeposit(bindingEnergy);
            //track.SetEnergyDeposit(edep);
        }
        // return with number of secondaries i.e. 1 photoelectron
        return 1;
    }
    
    
    void SauterGavrilaPhotoElectricModel::InitSamplingTables() {
        
        // set number of primary gamma energy grid points
        // keep the prev. value of primary energy grid points.
        int oldNumGridPoints = fNumSamplingPrimEnergies;
        
        fNumSamplingPrimEnergies = fNumSamplingPrimEnergiesPerDecade*std::lrint(std::log10(fMaxPrimEnergy/fMinPrimEnergy))+1;
        if (fNumSamplingPrimEnergies<2) {
            fNumSamplingPrimEnergies = 2;
        }
        
        // set up the initial gamma energy grid
        if (fSamplingPrimEnergies) {
            delete [] fSamplingPrimEnergies;
            delete [] fLSamplingPrimEnergies;
            fSamplingPrimEnergies  = nullptr;
            fLSamplingPrimEnergies = nullptr;
        }
        fSamplingPrimEnergies  = new double[fNumSamplingPrimEnergies];
        fLSamplingPrimEnergies = new double[fNumSamplingPrimEnergies];
        fPrimEnLMin    = std::log(fMinPrimEnergy);
        double delta   = std::log(fMaxPrimEnergy/fMinPrimEnergy)/(fNumSamplingPrimEnergies-1.0);
        fPrimEnILDelta = 1.0/delta;
        fSamplingPrimEnergies[0]  = fMinPrimEnergy;
        fLSamplingPrimEnergies[0] = fPrimEnLMin;
        fSamplingPrimEnergies[fNumSamplingPrimEnergies-1]  = fMaxPrimEnergy;
        fLSamplingPrimEnergies[fNumSamplingPrimEnergies-1] = std::log(fMaxPrimEnergy);
        for (int i=1; i<fNumSamplingPrimEnergies-1; ++i) {
            fLSamplingPrimEnergies[i] = fPrimEnLMin+i*delta;
            fSamplingPrimEnergies[i]  = std::exp(fPrimEnLMin+i*delta);
        }
        //
        // build the sampling tables at each primary gamma energy grid point.
        //
        // prepare the array that stores pointers to sampling data structures
        if (fAliasData) {
            for (int i=0; i<oldNumGridPoints; ++i) {
                if (fAliasData[i]) {
                    delete [] fAliasData[i]->fXdata;
                    delete [] fAliasData[i]->fYdata;
                    delete [] fAliasData[i]->fAliasW;
                    delete [] fAliasData[i]->fAliasIndx;
                    delete fAliasData[i];
                }
            }
            delete [] fAliasData;
        }
        // create new fAliasData array
        fAliasData = new LinAlias*[fNumSamplingPrimEnergies];
        for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
            fAliasData[i] = nullptr;
        }
        // create one sampling data structure at each primary gamma energy grid point:
        // -first create an AliasTable object
        if (fAliasSampler) {
            delete fAliasSampler;
        }
        // -the prepare each table one-by-one
        for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
            double egamma = fSamplingPrimEnergies[i];
            double kappa  = egamma/geant::kElectronMassC2;
            //std::cout<<"Creating AliasTable for: index "<<i<<" - energyGamma: "<<egamma<<" and egamma/geant::kElectronMassC2: "<<kappa<<std::endl;
            BuildOneLinAlias(i,kappa);
        }
    }
    
    //This method is calculating the differential cross section in the transformed variable xsi
    double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSectionLog(double tau, double xsi)
    {
        
        // Based on Geant4 : G4SauterGavrilaAngularDistribution
        // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
        // input  : energy0  (incoming photon energy)
        // input  : cosTheta (cons(theta) of photo-electron)
        // output : dsigma   (differential cross section, K-shell only)
        
        //double tau = energy0 / geant::kElectronMassC2;
        
        double cosTheta= 1-std::exp(xsi);
        //if(cosTheta>(1-1.e-12)) cosTheta=1.e-12;
        
        //gamma and beta: Lorentz factors of the photoelectron
        double gamma = tau + 1.0;
        double beta = std::sqrt(tau * (tau + 2.0)) / gamma;
        
        double z = 1 - beta * cosTheta;
        double z2 = z * z;
        double z4 = z2 * z2;
        double y = 1 - cosTheta * cosTheta; //sen^2(theta)
        
        double dsigmadcostheta = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
        double dsigmadxsi=dsigmadcostheta*(-std::exp(xsi));
        
        return dsigmadxsi;
        
    }
    
    
    double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSection(double tau, double cosTheta)
    {
        
        // Based on Geant4 : G4SauterGavrilaAngularDistribution
        // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
        // input  : energy0  (incoming photon energy)
        // input  : cosTheta (cons(theta) of photo-electron)
        // output : dsigma   (differential cross section, K-shell only)
        
        //double tau = energy0 / geant::kElectronMassC2;
        
        //gamma and beta: Lorentz factors of the photoelectron
        double gamma = tau + 1.0;
        double beta = std::sqrt(tau * (tau + 2.0)) / gamma;
        
        double z = 1 - beta * cosTheta;
        double z2 = z * z;
        double z4 = z2 * z2;
        double y = 1 - cosTheta * cosTheta; //sen^2(theta)
        
        double dsigma = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
        return dsigma;
        
    }
    
    void SauterGavrilaPhotoElectricModel::BuildOneLinAlias(int indx, double tau){
        
        fAliasData[indx]              = new LinAlias();
        fAliasData[indx]->fNumdata    = fNumSamplingAngles;
        fAliasData[indx]->fXdata      = new double[fNumSamplingAngles];
        fAliasData[indx]->fYdata      = new double[fNumSamplingAngles];
        fAliasData[indx]->fAliasW     = new double[fNumSamplingAngles];
        fAliasData[indx]->fAliasIndx  = new    int[fNumSamplingAngles];
        
        //fAliasData[indx]->fXdata[0] = std::log(1.e-12); //min
        //fAliasData[indx]->fXdata[1] = (std::log(2)+std::log(1.e-12))/2;
        //fAliasData[indx]->fXdata[2] = std::log(2);
        //fAliasData[indx]->fYdata[0] = CalculateDiffCrossSectionLog(tau, fAliasData[indx]->fXdata[0]);
        //fAliasData[indx]->fYdata[1] = CalculateDiffCrossSectionLog(tau, fAliasData[indx]->fXdata[1]);
        //fAliasData[indx]->fYdata[2] = CalculateDiffCrossSectionLog(tau, fAliasData[indx]->fXdata[2]);
        
        // note: the variable is a cos in [-1,1]
        // so fill 3 initial values of cos:
        //  -  xi_0 = x_min = -1.
        //  -  xi_1 = (x_max-x_min)/2 = 0.
        //  -  xi_2 = x_max = 1.
        // and the corresponding y(i.e.~PDF) values
        fAliasData[indx]->fXdata[0] = -1.;
        fAliasData[indx]->fXdata[1] = 0.;
        fAliasData[indx]->fXdata[2] = 1.0;
        fAliasData[indx]->fYdata[0] = CalculateDiffCrossSection(tau, -1.);
        fAliasData[indx]->fYdata[1] = CalculateDiffCrossSection(tau, 0.);
        fAliasData[indx]->fYdata[2] = CalculateDiffCrossSection(tau, 1.);
        
        
        int curNumData = 3;
        // expand the data up to numdata points
        while (curNumData<fAliasData[indx]->fNumdata) {
            // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
            double maxerr     = 0.0; // value of the current maximum error
            double thexval    = 0.0;
            double theyval    = 0.0;
            int    maxerrindx = 0;   // the lower index of the corresponding bin
            for (int i=0; i<curNumData-1; ++i) {
                double xx    = 0.5*(fAliasData[indx]->fXdata[i]+fAliasData[indx]->fXdata[i+1]);    // mid x point
                double yy    = 0.5*(fAliasData[indx]->fYdata[i]+fAliasData[indx]->fYdata[i+1]);    // lin. interpolated pdf value at the mid point
                //double val   = CalculateDiffCrossSectionLog(tau, xx); // real pdf value at the mid point
                double val   = CalculateDiffCrossSection(tau, xx); // real pdf value at the mid point
                double err   = std::abs(1.-(yy/val));
                if (err>maxerr) {
                    maxerr     = err;
                    maxerrindx = i;
                    thexval    = xx;
                    theyval    = val;
                }
            }
            // extend x,y data by puting a new real value at the mid point of the highest error bin
            // first shift all values to the right
            for (int j=curNumData; j>maxerrindx+1; --j) {
                fAliasData[indx]->fXdata[j] = fAliasData[indx]->fXdata[j-1];
                fAliasData[indx]->fYdata[j] = fAliasData[indx]->fYdata[j-1];
            }
            // fill x mid point
            fAliasData[indx]->fXdata[maxerrindx+1] = thexval;
            fAliasData[indx]->fYdata[maxerrindx+1] = theyval;
            // increase number of data
            ++curNumData;
        } // end while
        // prepare the alias data for this PDF(x,y)
        fAliasSampler->PreparLinearTable(fAliasData[indx]->fXdata, fAliasData[indx]->fYdata,
                                         fAliasData[indx]->fAliasW, fAliasData[indx]->fAliasIndx,
                                         fAliasData[indx]->fNumdata);
    }
    
}   // namespace geantphysics
