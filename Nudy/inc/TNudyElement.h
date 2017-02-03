#include<iostream>
#ifndef TNudyElement_H
#define TNudyElement_H
#include <string>
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif


class TNudyElement
{

public:
TNudyElement();
TNudyElement(std::string,int,int);
std::string GetatomicSymbol(void){ return symbol ;}
int GetatomicNumber( void ){return charge;}
int GetatomicMass( void ){return mass;}
double GetatomicDensity( void );
std::string replaceName(std::string &, std::string &, std::string &);
const char* endfFileName();

virtual ~TNudyElement();
private:
int mass;
int charge;
double rho;
std::string symbol;

#ifdef USE_ROOT
  ClassDef(TNudyElement, 1) //
#endif
};
#endif


