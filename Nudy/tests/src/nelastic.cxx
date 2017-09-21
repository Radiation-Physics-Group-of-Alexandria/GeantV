#include "nelastic.h"

using namespace Nudy;
using namespace NudyPhysics;

void setEnv() {
	gSystem->Load("libRIO");
	gSystem->Load("libGeom");		
	gSystem->Load("libNudy");
	gSystem->Load("libMathMore");
   gSystem->Load("libHist");
   gSystem->Load("libGraf");
   gSystem->Load("libGpad");
   gSystem->Load("libEG");

}

int main () {
	setEnv();

	int ielementID;
	double nuEnergy = 1.0e+6;
	double sigmaEn;
	std::vector<double> xsecE;

	std::string symbolT      = "Fe";
	int         tZ           = 26;
	int         tA           = 56;

	const char* irENDF;
	const char* rootENDF;
	const char* nENDFfileName;
	std::string ENDFNAMESTR;
	std::string ENDFNameTxt;
	std::string ENDFNameRoot;

	ENDFNAMESTR = findENDFFileName(symbolT, tZ, tA);
	ENDFNameTxt += ENDFNameRoot = ENDFNAMESTR;
	ENDFNameTxt+=".endf";
	ENDFNameRoot += ".root";
	nENDFfileName = ENDFNameTxt.c_str();
	rootENDF = ENDFNameRoot.c_str();
	
	TNudyENDF *proc = new TNudyENDF(nENDFfileName, rootENDF, "recreate");
	proc->Process();



	return 0;
}


