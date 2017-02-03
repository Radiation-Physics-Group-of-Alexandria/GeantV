#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "TNudyElement.h"
using namespace std;
#ifdef USE_ROOT
ClassImp(TNudyElement)
#endif

TNudyElement::TNudyElement(){}

TNudyElement::TNudyElement(std::string sym, int Z, int A)
{
  symbol = sym;
  charge = Z;
  mass   = A;
}
//______________________________________________________________________________
const char* TNudyElement::endfFileName()
{
  //Atomic number
  int targetZ = charge;
  stringstream ss;
  ss << targetZ;
  string strZ = ss.str();
  int len=strZ.length();
  std::string finalZ;
  if(len==1){
    finalZ="00" + strZ;
  } else{
      finalZ="0"+strZ;
    }
//Atomic mass
  int targetA = mass;
  stringstream ssa;
  ssa << targetA;
  string strA = ssa.str();
  int lenA=strA.length();
  std::string finalA;
  if(lenA==1){
    finalA="00" + strA;
  } else if(lenA==2){
     finalA="0"+strA;
    } else{
       finalA=strA;
      }
  if(targetA == 12)finalA = "000"; // for 12C, in the endf data, the file was named as n-006_C_000.endf
  std::string directoryName="/home/shiba/endffile/n-ddd_XX_fff.endf";
  std::string str1;
  std::string str2;
  std::string str3;
  std::string name;

  std::string s = "XX";
  std::string dummyZ = "ddd";
  std::string dummyA = "fff";
  std::string m = finalA;
  str1=replaceName(directoryName, s, symbol);
  str2=replaceName(str1, dummyZ, finalZ);
  str3=replaceName(str2, dummyA,finalA);
  const char* endfFile=str3.c_str();
  return endfFile;
}
//______________________________________________________________________________

std::string TNudyElement::replaceName(std::string &s, std::string &toReplace,
                       std::string &replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}
//______________________________________________________________________________
TNudyElement::~TNudyElement(){}
