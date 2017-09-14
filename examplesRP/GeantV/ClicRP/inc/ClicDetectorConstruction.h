//===--- ClicDetectorConstruction.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ClicDetectorConstruction.h
 * @brief Implementation of geometry for GeantV calorimeter prototype
 * @author S. Vallecorsa midified from testEm5 and CaloRP
 * @date August 1, 2017
 */
//===----------------------------------------------------------------------===//

#ifndef CALO_DETECTOR_CONSTRUCTION
#define CALO_DETECTOR_CONSTRUCTION
#define NEW_DETECTOR
#include "GeantRunManager.h"
#include "GeantVDetectorConstruction.h"
#include "Material.h"
#include "Region.h"
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
class GeantRunManager;
class Material;
}
}
namespace userapplication {
class ClicDetectorConstruction : public Geant::GeantVDetectorConstruction {

static const int maxAbsorbers = 25;
  private:
	Geant::GeantRunManager *fRunMgr = nullptr;

	std::string fWorldMaterialName;
	std::string fCellMaterialName;
        std::string fAbsMaterialName[maxAbsorbers]; 
	
	geantphysics::Material *fAbsMaterial[maxAbsorbers];
	geantphysics::Material *fWorldMaterial;	
	geantphysics::Material *fCellMaterial;	
	bool userLayerNum=false;
        bool userAbsorberNum=false;
        bool userCaloXY=false;
	bool userThickness[maxAbsorbers];
	bool userMaterial[maxAbsorbers];

	bool prodCutByLength=true;

        int numAbsorbers;
        int numAbsorbers1;
        int numAbsorbers2;
        int numLayers;
        int numRows;
        int numCells;

	int fAbsLogicVolumeID[maxAbsorbers];	
	int fDetectorRegionIndex;
	double fGammaCut=0.1;
	double fElectronCut=0.1;
	double fPositronCut=0.1;

        double fAbsThickness[maxAbsorbers];
        double fCellThickness;
        double fTotThickness;
        double fCellX;
        double fCellY;

        double fCaloSizeXY;
        double fWorldSizeZ;
	double fWorldSizeXY;	

	void SetDetectorMaterials();

  public:
	ClicDetectorConstruction(Geant::GeantRunManager *runmgr) : GeantVDetectorConstruction(runmgr) {}
	~ClicDetectorConstruction();

  public:
        void SetAbsorberMaterialName(int,std::string);
        std::string GetAbsorberMaterialName(int);
        geantphysics::Material* GetAbsorberMaterial(int);

        void SetNumLayers(int);
        void SetNumAbsorbers(int);
        int GetNumLayers();
        int GetNumAbsorbers();

	void SetProductionCutsByEnergy(double);
	void SetProductionCutsByLength(double);
	void SetDetectorGammaProductionCut(double);
	void SetDetectorElectronProductionCut(double);
	void SetDetectorPositronProductionCut(double);

        void SetAbsorberThickness(int,double);
        double GetAbsorberThickness(int);

	int GetAbsorberLogicalVolumeID(int);
	int GetDetectorRegionIndex();

        void SetDetectorXY(double);
        double GetDetectorZ();
        double GetDetectorXY();
        double GetWorldZ();
        double GetWorldXY();

	void CreateMaterials();
	void CreateGeometry();
};
}
#endif
