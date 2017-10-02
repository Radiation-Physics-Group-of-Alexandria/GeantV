/**********************************************
 * @file NudyInterface.cxx
 * @author Abhijit Bhattacharyya
 * @brief This interface links to the GV through NUDYCrossSection of realphysics/physics/HAD/CrossSection
 **********************************************/
#include <iostream>
#include "NudyInterface.h"

using namespace geantphysics;
using namespace Nudy;
using namespace NudyPhysics;



  NudyInterface::NudyInterface() : fProjCode(2112), fProjKE(1.0*geant::MeV), fIsoName("Fe"), ftZ(26), ftN(56) {};

  NudyInterface::NudyInterface(
			       const int projCode, const double projKE, const std::string isoName,
			       const int tZ, const int tN
			       ) : fProjCode(projCode), fProjKE(projKE), fIsoName(isoName), ftZ(tZ), ftN(tN) {};


  NudyInterface::~NudyInterface() {}


  double NudyInterface::GetXS( int projCode, double projKE, std::string isoName, int tZ, int tN) {
    double XSvalue = 0.0;
    std::cout << std::getenv("ENDFDATADIR") << std::endl;

    return XSvalue;
  }


  std::string NudyInterface::findENDFFileName(std::string elementName, int tZ, int tA) {
    std::stringstream ss;
    std::string fName="n-";
    ss << tZ;
    std::string stZ = ss.str();
    ss.str("");
    ss << tA;
    std::string stA = ss.str();
    ss.str("");

    switch(stZ.length()) {
    case 0: return "";
    case 1: stZ = "00" + stZ;   break;
    case 2: stZ = "0" + stZ;    break;
    case 3: stZ = stZ;
    }

    if (tZ == 12) {
      stA = "000";
    } else {
      switch(stA.length()) {
      case 0: return "";
      case 1: stA = "00" + stA;  break;
      case 2: stA = "0" + stA;   break;
      case 3: stA = stA;
      }
    }

    //  fName = "/home/vega/ENDFDATA/ENDF-B-VII.1/neutrons/" + fName + stZ + "_" + ele + "_" + stA;
    fName = fName + stZ + "_" + elementName + "_" + stA;
    return fName;
  }
