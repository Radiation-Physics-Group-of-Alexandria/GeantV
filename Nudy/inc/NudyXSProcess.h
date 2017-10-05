//-----------------------------------------------------------------
// @file NudyXSProcess.h
// @brief prototype to caller functions for Nudy called by Nudy interface for GV
// @author Abhijit Bhattacharyya
//----------------------------------------------------------------
#ifndef NUDY_XS_PROCESS_H
#define NUDY_XS_PROCESS_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>

#include "Material.h"
#include "Element.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#ifdef USE_ROOT
#include "TNudyDB.h"
#include "TNudyAlias.h"
#include "TNudyElementTable.h"
#include "TNudyENDF.h"
#include "TNudyEndfCont.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfINTG.h"
#include "TNudyEndfList.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfRecord.h"
#include "TNudyEndfSec.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyLibrary.h"
#include "TNudyManager.h"
#include "TNudySubLibrary.h"
#include "TVNudyModel.h"
#include "TNudyEndfEnergy.h"
#include "TNudyEndfSigma.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfAng.h"

#include "TSystem.h"
//class TFile;
#endif


using Nudy::TNudyAlias;
using Nudy::TNudyCore;
using Nudy::TNudyDB;
using Nudy::TNudyElementTable;
using Nudy::TNudyEndfFile;
using Nudy::TNudyENDF;
using Nudy::TNudyEndfMat;
using Nudy::TNudyEndfList;
using Nudy::TNudyEndfSec;
using Nudy::TNudyEndfCont;
using Nudy::TNudyEndfRecord;
using Nudy::TNudyEndfTape;
using Nudy::TNudyEndfTab1;
using Nudy::TNudyEndfTab2;
using Nudy::TNudyLibrary;
using Nudy::TNudyManager;
using Nudy::TVNudyModel;
using Nudy::TNudyEndfINTG;

using NudyPhysics::TNudyEndfDoppler;
using NudyPhysics::TNudyEndfEnergyAng;
using NudyPhysics::TNudyEndfEnergy;
using NudyPhysics::TNudyEndfFissionYield;
using NudyPhysics::TNudyEndfNuPh;
using NudyPhysics::TNudyEndfPhAng;
using NudyPhysics::TNudyEndfPhEnergy;
using NudyPhysics::TNudyEndfPhProd;
using NudyPhysics::TNudyEndfPhYield;
using NudyPhysics::TNudyEndfRecoPoint;
using NudyPhysics::TNudyEndfSigma;


namespace geantphysics {
  class LightTrack;
  class HadronicCrossSection;
  // class HadronicProcess;
  inline namespace GEANT_IMPL_NAMESPACE {
    class Isotope;
    class Material;
    class Element;
  }
}


namespace NudyPhysics{
  class NudyXSProcess {
  public:
    NudyXSProcess();
    NudyXSProcess(const int projCode, const double projKE, const double temp, const std::string isoName,
      const int tZ, const int tN, std::string reactType );
    virtual ~NudyXSProcess();

  public:
    inline std::string GetIsotopeName ();
    inline int GetProjectileCode ();
    inline int GetZ ();
    inline int GetN ();
    inline double GetProjectileKE ();
    inline double GetTemp ();
    inline double GetCrossSection ();
    inline std::string GetReactionType();

    inline void SetIsotopeName ( const std::string &isoName );
    inline void SetProjectileCode ( const int projCode );
    inline void SetZ ( const int tZValue );
    inline void SetN ( const int tNValue ) ;
    inline void SetA ( const int tNvalue, const int tZvalue );
    inline void SetProjectileKE ( const double projKE );
    inline void SetTemp (const double temp );
    inline void SetCrossSection ( const double XSvalue );
    inline void SetReactionType ( const std::string reactType );
    inline void SetEndfDataFileName ( const char * fileName );
    inline void SetEndfSubDataFileName ( const char * fileName );
    inline void SetRootFileName ( const char * fileName );

    double GetXS( int projCode, double projKE, double temp, std::string isoName, int tZ, int tN, std::string reactType ) ;
    std::string findENDFFileName( std::string projID, int projCode, std::string ele, int tZ, int tN ) ;
    std::string GetDataFileName( std::string str1, int projCode, std::string str2 ); // projID, isoName
    std::string FixRootDataFile( std::string str1 );                   // ENDF filename without path and extension
    std::string GetCWD();
    double ProcFission();

  private :
    std::string fIsoName;
    int fProjCode;
    int ftZ;
    int ftN;
    int ftA;
    double fProjKE;
    double fTemperature;
    double fXS;
    const char* fEndfDataFileName;
    const char* fEndfSubDataFileName;
    const char* fRootFileName;
    std::string fReactType;
    std::vector<int> fFissionFragmentsMass;
    std::vector<int> fFissionFragmentCharge;
  };

  //--------- GETTERS -------
  inline std::string NudyXSProcess::GetIsotopeName () { return fIsoName; }
  inline int NudyXSProcess::GetProjectileCode () { return fProjCode; }
  inline int NudyXSProcess::GetZ () { return ftZ; }
  inline int NudyXSProcess::GetN () { return ftN; }
  inline double NudyXSProcess::GetProjectileKE () { return fProjKE; }
  inline double NudyXSProcess::GetTemp () { return fTemperature; }
  inline double NudyXSProcess::GetCrossSection () { return fXS; }
  inline std::string NudyXSProcess::GetReactionType () { return fReactType; }

  //--------- SETTERS ---------
  inline void NudyXSProcess::SetIsotopeName ( const std::string &isoName ) { fIsoName = isoName; }
  inline void NudyXSProcess::SetProjectileCode ( const int projCode ) { fProjCode = projCode; }
  inline void NudyXSProcess::SetZ ( const int tZValue ) { ftZ = tZValue; }
  inline void NudyXSProcess::SetN ( const int tNValue ) { ftN = tNValue; }
  inline void NudyXSProcess::SetA ( const int tNvalue, const int tZvalue ) { ftA = ftZ + ftN; }
  inline void NudyXSProcess::SetProjectileKE ( const double projKE ) { fProjKE = projKE; }
  inline void NudyXSProcess::SetTemp ( const double temp ) { fTemperature = temp; }
  inline void NudyXSProcess::SetCrossSection ( const double XSvalue ) { fXS = XSvalue; }
  inline void NudyXSProcess::SetReactionType ( const std::string reactType ) { fReactType = reactType; }
  inline void NudyXSProcess::SetEndfDataFileName ( const char * fileName ) { fEndfDataFileName = fileName; }
  inline void NudyXSProcess::SetEndfSubDataFileName ( const char * fileName ) { fEndfSubDataFileName = fileName; }
  inline void NudyXSProcess::SetRootFileName ( const char * fileName ) { fRootFileName = fileName; }

} // namespace NudyPhysics

#endif
