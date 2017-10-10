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


NudyInterface::NudyInterface() //:
// fProjCode(2112), fProjKE(4.0e+6), fTemperature(293.60608), fIsoName("Fe"), ftZ(26), ftN(56),
// fEndfDataFileName(""), fEndfSubDataFileName(""), fRootFileName("")
{};


NudyInterface::NudyInterface(
  const int projCode, const double projKE, double temp, const std::string isoName, const int tZ, const int tN)
  // : fProjCode(projCode), fProjKE(projKE), fIsoName(isoName), ftZ(tZ), ftN(tN) {
  {
    std::string newReact;
    SetProjectileCode (projCode);
    SetProjectileKE (projKE);
    SetTemp (temp);
    SetIsotopeName (isoName);
    SetZ (tZ);
    SetN (tN);
    SetA (tN, tZ);
    NudyInterface::InitializeXSChannel();
  };

NudyInterface::~NudyInterface() {}

////////////////////////////////////
std::vector<double> NudyInterface::GetXS( int projCode, double projKE, double temp, std::string isoName, int tZ, int tN ) {
  NudyInterface::InitializeXSChannel();
  NudyInterface::SetA (tN, tZ);
  NudyInterface::SetProjectileKE (projKE);
  NudyInterface::SetTemp(temp);
  NudyInterface::SetIsotopeName(isoName);
  NudyInterface::SetIsFissKey(false); // initializing to false first

  Nudy::TNudyENDF *proc;

  //  Fix the name of the ENDF, ROOT and ENDFSUB filenames here
  std::string fileENDF1 = NudyInterface::SetDataFileNameENDF(projCode, isoName, projKE, tZ, tN);
  SetEndfDataFileName(fileENDF1.c_str());
  std::string fileENDF2 = NudyInterface::SetDataFileNameROOT(isoName, tZ, tN);
  SetRootFileName (fileENDF2.c_str());
  std::string fileENDF3 = NudyInterface::SetDataFileNameENDFSUB( isoName, tZ, tN);
  NudyInterface::SetEndfSubDataFileName (fileENDF3.c_str());

  // Create and process with NUDY with keywords
  proc = new Nudy::TNudyENDF (fEndfDataFileName, fRootFileName, "recreate");
  proc->Process();
  bool LFIval = proc->GetLFI();
  SetIsFissKey(LFIval);
  proc->SetEndfSub(fEndfSubDataFileName);
  proc->Process();

 // Processing starts for channels
  for (unsigned int channel = 2; channel <= NudyInterface::fNumberOfReactionChannels; channel++) {
    bool isCHFiss = NudyInterface::GetFisCha(channel);
    // std::cout << " CH # " << channel << " LFI :: " << LFIval << "  ChFiss : " << isCHFiss << std::endl;
    if (LFIval && !isCHFiss) continue;
    SetProjIDFn(projCode, projKE, channel);
    NudyInterface::ComputeCrossSection(channel);  // call routines from Nudy
  }
  return fChannelXSArray;
}

void NudyInterface::InitializeXSChannel() {
  fChannelXSArray.resize(fNumberOfReactionChannels, -999.99); // any value to indicate error
}

std::string NudyInterface::SetDataFileNameENDF( int projCode, std::string isoName, double projKE, int tZ, int tN ){
  NudyInterface::SetProjIDFn( projCode, projKE, 0);
  std::string DataFileNameString = NudyInterface::findENDFFileName( isoName, tZ, tN);
  std::string fileENDFName1 = NudyInterface::GetDataFileName(fProjID, DataFileNameString);
  return fileENDFName1;
}

std::string NudyInterface::SetDataFileNameROOT( std::string isoName, int tZ, int tN ){
  std::string DataFileNameString = NudyInterface::findENDFFileName(isoName, tZ, tN);
  std::string fileENDFName2 = NudyInterface::FixRootDataFile(DataFileNameString);
  return fileENDFName2;
}

void NudyInterface::SetProjIDFn(int projCode, double projKE, unsigned int channel) {
    bool isChFiss = NudyInterface::GetFisCha(channel);
    if ( projKE < 0.03*geant::eV && projCode == 2112 ) NudyInterface::SetProjID("thermal_scatt");
    if ( projKE > 0.03*geant::eV && projCode == 2112 ) NudyInterface::SetProjID("neutrons");
    if (isChFiss) NudyInterface::SetProjID("nfy");
  }

std::string NudyInterface::SetDataFileNameENDFSUB( std::string isoName, int tZ, int tN ){
  std::string fileENDFName3;
  std::string DataFileNameString;
  NudyInterface::SetProjID("nfy");
  std::string prjId = NudyInterface::GetProjID();
  DataFileNameString = NudyInterface::findENDFFileName(isoName, tZ, tN);
  fileENDFName3 = NudyInterface::GetDataFileName(prjId, DataFileNameString);
  return fileENDFName3;
}

// Actual Nudy CrossSection computation method
void NudyInterface::ComputeCrossSection(unsigned int channel) {
  int iElementID = 0;  //<------------- confusing testing by Abhijit 419 ?
  double xsvalue = 0.0;
  double iSigDiff = 0.001;   // trial value for test documentation reqd.
  unsigned int MT;
  //std::vector<double> crs;

  NudyPhysics::TNudyEndfSigma();
  NudyPhysics::TNudyEndfSigma xsec;
  xsec.GetData(fRootFileName, iSigDiff);

  NudyPhysics::TNudyEndfRecoPoint *recoPoint = new NudyPhysics::TNudyEndfRecoPoint(iElementID, fRootFileName);
  // This is  under testing to check Harphool code for interfacing to GV :: Abhijit

  for ( unsigned int crsp = 0; crsp < recoPoint->MtValues[iElementID].size(); crsp++) {
    MT = recoPoint->MtValues[iElementID][crsp]; //[isel];
    if ( MT == channel) {
      xsvalue = recoPoint->GetSigmaTotal(iElementID, fProjKE);
      break;
    }
  }
  NudyInterface::AppendXS(xsvalue);
}

  // selects the data file name for ENDF data
  std::string NudyInterface::GetDataFileName(std::string str1, std::string str2) {
    std::string EndfDataPath="";
    if (std::getenv("ENDFDATADIR")!=NULL){
      EndfDataPath = std::getenv("ENDFDATADIR");
      std::string ENDFDataString = EndfDataPath + "/" + str1 + "/" + str2 + ".endf";
      return ENDFDataString;
  } else {
    std::cout << " Please set environment ENDFDATADIR pointing to ENDF data directory." << std::endl;
    exit(99);
  }
  return EndfDataPath;
}

  // selects name of root data file in the current working directory
  std::string NudyInterface::FixRootDataFile(std::string str1){
    std::string cwdPath = NudyInterface::GetCWD();
    std::string rootENDFFile = cwdPath + "/" + str1 + ".root";
    return rootENDFFile;
  }

  // store root data file in current working directory
  // This is a bad technique. Using for the testing as quick solution.
  std::string NudyInterface::GetCWD() {
    char * tempDIRarray = new char[1024];
    std::string cwdPath= (getcwd(tempDIRarray, 1024)) ? std::string(tempDIRarray) : std::string("");
    delete [] tempDIRarray;
    return cwdPath;
  }

  std::string NudyInterface::findENDFFileName(std::string elementName, int tZ, int tN) {
    int tA = tZ + tN; // this is mass number to be derived from Z and N number
    std::stringstream ss;
    std::string fName = "";
    if (fProjID == "thermal_scatt") {
      fName = "tsl-";
    } else if (fProjID == "neutrons" ) {
      fName = "n-";
    } else {
      fName = fProjID + "-";
    }

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


    bool NudyInterface::GetFisCha(int inKey) {
      bool isIn = false;
      isIn = std::find(
        NudyInterface::fChannelFiss.begin(), NudyInterface::fChannelFiss.end(), inKey
      ) != NudyInterface::fChannelFiss.end();
      return isIn;
    }
