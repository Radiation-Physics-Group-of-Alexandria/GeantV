/****************************************************
 * @file NudyXSProcess.cxx
 * @author Abhijit Bhattacharyya
 * @brief Supplier of useful functions for NudyInterface
 ****************************************************/
 #include <iostream>
 #include "NudyXSProcess.h"

 using namespace geantphysics;
 using namespace Nudy;
 using namespace NudyPhysics;

NudyXSProcess::NudyXSProcess() :
 fProjCode(2112), fProjKE(1.0*geant::MeV), fTemperature(293.60608), fIsoName("Fe"), ftZ(26), ftN(56), fReactType("Total") {};

NudyXSProcess::NudyXSProcess(
  const int projCode, const double projKE, const double temp, const std::string isoName,
  const int tZ, const int tN, const std::string reactType )
   // :  fProjCode(projCode), fProjKE(projKE), fIsoName(isoName), ftZ(tZ), ftN(tN)
   {
     SetProjectileCode (projCode);
     SetProjectileKE (projKE);
     SetTemp (temp);
     SetIsotopeName (isoName);
     SetZ (tZ);
     SetN (tN);
     SetReactionType(reactType);
   }

NudyXSProcess::~NudyXSProcess() {}

double NudyXSProcess::GetXS( int projCode, double projKE, double temp, std::string isoName, int tZ, int tN, std::string reactType ) {
  double XSvalue = 0.0;
  std::string projID = "";

  NudyXSProcess::SetA (tZ, tN);
  NudyXSProcess::SetProjectileKE (projKE);

  if (projCode == 2112) projID = "neutrons";   // this is for the neutrons at the time being for testing
  if (reactType == "Fission") projID = "nfy";
  if (reactType == "Thermal") projID = "thermal_scatt";

  // This filename deduction is general
  std::string DataFileNameString = NudyXSProcess::findENDFFileName(projID, projCode, isoName, tZ, tN);
  std::string fileENDF1 = NudyXSProcess::GetDataFileName(projID, projCode, DataFileNameString);
  SetEndfDataFileName (fileENDF1.c_str());

  std::string fileENDF2 = NudyXSProcess::FixRootDataFile(DataFileNameString);
  SetRootFileName (fileENDF2.c_str());

  std::string fileENDF3;
  if (reactType == "Fission"){
    DataFileNameString = NudyXSProcess::findENDFFileName(projID, 0, isoName, tZ, tN);
    fileENDF3 = NudyXSProcess::GetDataFileName(projID, 0, DataFileNameString);
    SetEndfSubDataFileName (fileENDF3.c_str());
  }

  XSvalue = NudyXSProcess::ProcFission();  // call routines for fission from Nudy
  return XSvalue;
}

double NudyXSProcess::ProcFission() {
  int iElementID = 0;  //<------------- confusing testing by Abhijit

  double xsvalue = 0.0;
  double iSigDiff = 0.001;

  double sigmaTotal = 0.0;
  double sigmaPartial = 0.0;

  Nudy::TNudyENDF *proc = new Nudy::TNudyENDF(fEndfDataFileName, fRootFileName, "recreate");
  proc->Process();
  proc->SetEndfSub(fEndfSubDataFileName);
  proc->Process();
  NudyPhysics::TNudyEndfSigma();
  NudyPhysics::TNudyEndfSigma xsec;
  xsec.GetData(fRootFileName, iSigDiff);

  NudyPhysics::TNudyEndfRecoPoint *recoPoint = new NudyPhysics::TNudyEndfRecoPoint(iElementID, fRootFileName);

  // This is  under testing to check Harphool code for interfacing to GV :: Abhijit
  int isel = 0;
  int MT;
  double sumRnd = 0.0;
  std::vector<double> crs;
  TRandom3 *fRnd = new TRandom3(0);
  double randomSeries = fRnd->Uniform();
  std::cout << " The projectile KE:-> " << fProjKE << std::endl;
  for ( unsigned int crsp = 0; crsp < recoPoint->MtValues[iElementID].size(); crsp++) {
    MT = recoPoint->MtValues[iElementID][crsp]; //[isel];
    if (MT == 18) { // Fission
      xsvalue = recoPoint->GetNuTotal(iElementID, fProjKE);
      break;
    }
  }
  return xsvalue;
}


// selects the data file name for ENDF data
std::string NudyXSProcess::GetDataFileName(std::string str1, int Pcode, std::string str2) {
  std::string EndfDataPath = std::getenv("ENDFDATADIR");
  if (Pcode == 2112) str1 = "neutrons";
  std::string ENDFDataString = EndfDataPath + "/" + str1 + "/" + str2 + ".endf";
  return ENDFDataString;
}

// selects name of root data file in the current working directory
std::string NudyXSProcess::FixRootDataFile(std::string str1){
  std::string cwdPath = NudyXSProcess::GetCWD();
  std::string rootENDFFile = cwdPath + "/" + str1 + ".root";
  return rootENDFFile;
}

// store root data file in current working directory
// This is a bad technique. Using for the testing as quick solution.
std::string NudyXSProcess::GetCWD() {
  char * tempDIRarray = new char[1024];
  std::string cwdPath= (getcwd(tempDIRarray, 1024)) ? std::string(tempDIRarray) : std::string("");
  delete [] tempDIRarray;
  return cwdPath;
}

std::string NudyXSProcess::findENDFFileName(std::string projID, int projCode, std::string elementName, int tZ, int tN) {
  int tA = tZ + tN; // this is mass number to be derived from Z and N number
  std::stringstream ss;
  std::string fName = "";
  if (projID == "thermal_scatt") {
    fName = "tsl-";
  } else if (projID == "neutrons" || projCode == 2112) {
    fName = "n-";
  } else {
    fName = projID + "-";
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
