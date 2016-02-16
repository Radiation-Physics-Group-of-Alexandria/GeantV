//
// Create Detector geometry, sensitive detectors and (magnetic) field
//
// Started from the following Geant4 class
//    examples/extended/persistency/gdml/G01/include/G01DetectorConstruction.hh
// and adjusted substanially to accommodate:
//    i) preconstructed geometry
//   ii) uniform or class-defined magnetic field
//
//    Authors: Mihaly Novak, CERN            March 2015
//             John Apostolakis, CERN        February 2016

#ifndef DETECTORCONSTRUCTION_H
#define DETECTORCONSTRUCTION_H 1

#include "G4VUserDetectorConstruction.hh"

class TGeoManager;
class G4MagneticField;
class G4UniformMagField;
class CMSMagneticFieldG4;

/// Detector construction allowing to use the geometry read from the GDML file

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction(G4VPhysicalVolume *setWorld = 0, bool useUniformField= true );
    ~DetectorConstruction();

    virtual G4VPhysicalVolume *Construct()
    {
      return fWorld;
    }

   virtual void ConstructSDandField();

   void SetUniformBzMagField(G4double);               // Set uniform magnetic field
   // static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 

   // void SetMagField(G4MagneticField* pFld);  // Set non-uniform magnetic field
   void UseInterpolatedField() { fUseUniformField = false; }

   // void SetDistanceConst( double dist ){ fDistanceConst= dist; }

  private:

    G4VPhysicalVolume  *fWorld;
    static TGeoManager *fgGeomMgrRoot;

    bool                fUseUniformField;
    double              fBzFieldValue;
    G4MagneticField    *fpMagField;           //pointer to the magnetic field

    std::string         fFieldFilename;
    // double              fDistanceConst;   // Field is assumed to be uniform over this distance
    static CMSMagneticFieldG4* fpMasterCmsField;  // Object in the master -- others are its 'clones'
};

#endif
