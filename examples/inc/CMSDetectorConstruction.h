#ifndef GEANT_CMS_Detector_Construction
#define GEANT_CMS_Detector_Construction

#include "UserDetectorConstruction.h"
#include "CMSDetectorConstruction.h"

#ifdef  USE_ROOT_TObject
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#endif

class CMSmagField;

class CMSDetectorConstruction : public UserDetectorConstruction 
#ifdef USE_ROOT_TObject
   , public TObject
#endif
{
  public:
    /** @brief Destructor */
    CMSDetectorConstruction();
    CMSDetectorConstruction(const char* fieldFilename);
    CMSDetectorConstruction(std::string fieldFilename);
    ~CMSDetectorConstruction();
    
    /** @brief Destructor */
    void SetFileForField(const char *filename){ fFieldFilename= filename; }
    void SetFileForField(std::string filename){ fFieldFilename= filename; }

    /** @brief Method to register a B-field, and create integrator for it. */
    bool CreateFieldAndSolver(bool useRungeKutta= true);  // override final;

  private:
    std::string   fFieldFilename;
    CMSmagField*  fCMSfield;
    TUniformMagField*  fUniformField; // Alternative - for debugging only
    /** Field is created and owned by this class */

    // ClassDef(CMSDetectorConstruction, 1) // User application    
};
#endif
