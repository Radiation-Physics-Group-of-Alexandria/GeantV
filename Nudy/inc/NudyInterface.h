//-----------------------------------------------------------------
// @file NudyInterface.h
// @brief prototype Nudy interface for GV
// @author Abhijit Bhattacharyya
//----------------------------------------------------------------
#ifndef NUDY_INTERFACE_H
#define NUDY_INTERFACE_H

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
//#include "NudyXSProcess.h"

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
//using NudyPhysics::NudyXSProcess;


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
  class NudyInterface {

  public:
    NudyInterface () ;  
    NudyInterface( const int projCode, const double projKE, const double temp, const std::string isoName,
      const int tZ, const int tN );
    virtual ~NudyInterface ();

    std::vector<double> GetXS( int projCode, double projKE, double temp, std::string isoName, int tZ, int tN ) ;
    std::string  SetDataFileNameENDF( int projCode, std::string isoName, double projKE, int tZ, int tN );
    std::string findENDFFileName( std::string ele, int tZ, int tN ) ;
    std::string GetDataFileName( std::string str1, std::string str2 ); // projID, isoName
    std::string FixRootDataFile( std::string str1 );                   // ENDF filename without path and extension
    std::string SetDataFileNameROOT( std::string isoName, int tZ, int tN );
    std::string SetDataFileNameENDFSUB( std::string isoName, int tZ, int tN );
    std::string GetCWD();
    void SetProjIDFn( int projCode, double projKE, unsigned int channel);
    void ComputeCrossSection(  unsigned int channel );
    void InitializeXSChannel();
    bool GetFisCha(int inKey);

  public:
    inline std::string GetIsotopeName ();
    inline int GetProjectileCode ();
    inline int GetZ ();
    inline int GetN ();
    inline double GetProjectileKE ();
    inline double GetTemp();
    inline double GetCrossSection ();
    inline bool GetIsFissKey();
    inline std::string GetProjID();
    inline std::vector<double> GetXSTable();

    inline void SetIsotopeName (const std::string &isoName );
    inline void SetProjectileCode ( const int projCode );
    inline void SetZ ( const int tZValue );
    inline void SetN ( const int tNValue ) ;
    inline void SetA ( const int tNvalue, const int tZvalue );
    inline void SetProjectileKE ( const double projKE );
    inline void SetTemp ( const double temp );
    inline void SetCrossSection ( const double XSvalue );
    inline void SetEndfDataFileName ( const char * fileName );
    inline void SetEndfSubDataFileName ( const char * fileName );
    inline void SetRootFileName ( const char * fileName );
    inline void SetIsFissKey( const bool theKey);
    inline void SetProjID (const std::string &theID);
    inline void AppendXS ( const double xsvalue );

  private :
    unsigned int fNumberOfReactionChannels = 895;
    std::string fIsoName;
    int fProjCode;
    int ftZ;
    int ftN;
    int ftA;
    double fProjKE;
    double fTemperature;
    double fXS;
    std::string fProjID;
    const char* fEndfDataFileName;
    const char* fEndfSubDataFileName;
    const char* fRootFileName;
    std::vector<int> fFissionFragmentsMass;
    std::vector<int> fFissionFragmentCharge;
    std::vector<double>  fChannelXSArray;
    bool fIsFiss;
    std::vector<int> fChannelFiss {18, 19, 20, 21, 38, 151, 452, 454, 455, 456, 457, 458, 459};
  };

  //--------- GETTERS -------
  inline std::string NudyInterface::GetIsotopeName () { return fIsoName; }
  inline int NudyInterface::GetProjectileCode () { return fProjCode; }
  inline int NudyInterface::GetZ () { return ftZ; }
  inline int NudyInterface::GetN () { return ftN; }
  inline double NudyInterface::GetProjectileKE () { return fProjKE; }
  inline double NudyInterface::GetTemp () { return fTemperature; }
  inline bool NudyInterface::GetIsFissKey() { return fIsFiss; }
  inline std::string NudyInterface::GetProjID() { return fProjID; }
  inline std::vector<double> NudyInterface::GetXSTable() { return fChannelXSArray; }

  //--------- SETTERS ---------
  inline void NudyInterface::SetIsotopeName ( const std::string &isoName ) { fIsoName = isoName; }
  inline void NudyInterface::SetProjectileCode ( const int projCode ) { fProjCode = projCode; }
  inline void NudyInterface::SetZ ( const int tZValue ) { ftZ = tZValue; }
  inline void NudyInterface::SetN ( const int tNValue ) { ftN = tNValue; }
  inline void NudyInterface::SetA ( const int tNvalue, const int tZvalue ) { ftA = tNvalue + tZvalue; }
  inline void NudyInterface::SetProjectileKE ( const double projKE ) { fProjKE = projKE; }
  inline void NudyInterface::SetTemp (const double temp ) { fTemperature = temp; }
  inline void NudyInterface::SetCrossSection ( const double XSvalue ) { fXS = XSvalue; }
  inline void NudyInterface::SetEndfDataFileName ( const char * fileName ) { fEndfDataFileName = fileName; }
  inline void NudyInterface::SetEndfSubDataFileName ( const char * fileName ) { fEndfSubDataFileName = fileName; }
  inline void NudyInterface::SetRootFileName ( const char * fileName ) { fRootFileName = fileName; }
  inline void NudyInterface::SetIsFissKey( const bool theKey ) { fIsFiss = theKey; }
  inline void NudyInterface::SetProjID (const std::string &theID ) { fProjID = theID; }
  inline void NudyInterface::AppendXS ( const double xsvalue ) { fChannelFiss.push_back(xsvalue); }

 } // namespace NudyPhysics

#endif
