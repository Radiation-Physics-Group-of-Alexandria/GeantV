//-----------------------------------------------------------------
// testing the n elastic xsec produced by Nudy
// @file nelastic.h
// @brief tests neutron elastic cross section using nudy without GV
// @Abhijit Bhattacharyya
//----------------------------------------------------------------
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>

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
#include "Particle.h"
#include "G4ParticleTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "TSystem.h"


class TFile;
using namespace Nudy;
using namespace NudyPhysics;

namespace NudyPhysics {
    class TNudyEndfEnergy;
    class TNudyEndfSigma;
	class TNudyEndfEnergyAng;
    class TNudyEndfRecoPoint;
    class TNudyEndfAng;
	class TNudyEndfDoppler;
	class TNudyEndfNuPh;
	class TNudyEndfFissionYield;
	class TNudyEndfPhAng;
}

namespace Nudy {
	class TNudyDB; 
	class TNudyElement;
	class TNudyEndfTape;
    class TNudyEndfMat;
    class TNudyEndfFile;
    class TNudyEndfSec;
    class TNudyEndfCont;
    class TNudyEndfList;
    class TNudyEndfTab1;
    class TNudyEndfTab2;
    class TNudyEndfINTG;
	class TNudyEndfAlias;
	class TNudyAliasCont;
    class TNudyENDF;
}

	// inline namespace GEANT_IMPL_NAMESPACE {
	//	class Isotope;
	//	class Material;
	//	class Element;
	//}


string findENDFFileName(string ele, int tZ, int tA) {
	std::stringstream ss;
	string fName="n-";
	ss << tZ;
	string stZ = ss.str();
	ss.str("");
	ss << tA;
	string stA = ss.str();
	ss.str("");

	switch(stZ.length()) {
		case 0: return "";
		case 1: stZ = "00" + stZ;
			break;
		case 2: stZ = "0" + stZ;
			break;
		case 3: stZ = stZ;
	}

	if (tZ == 12) {
		stA = "000";
	} else {
		switch(stA.length()) {
			case 0: return "";
			case 1: stA = "00" + stA;
				break;
			case 2: stA = "0" + stA;
				break;
			case 3: stA = stA;
		}
	}

	fName = "/home/vega/ENDFDATA/ENDF-B-VII.1/neutrons/" + fName + stZ + "_" + ele + "_" + stA;
	return fName;
}






