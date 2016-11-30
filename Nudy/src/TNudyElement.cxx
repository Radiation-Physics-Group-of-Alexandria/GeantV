#include <iostream>
#include "TNudyElement.h"
using namespace std;
#ifdef USE_ROOT
ClassImp(TNudyElement)
#endif

TNudyElement::TNudyElement(){}

TNudyElement::TNudyElement(std::string sym, int Z, int A){
symbol = sym;
charge = Z;
mass   = A;
}
TNudyElement::~TNudyElement(){}






