

#include "GUVVecField.h"
#include <iostream>

 
GUVVecField& GUVVecField::operator = (const GUVVecField &p)
{
  
   if (&p != this){
   
     this->fChangesEnergy= p.fChangesEnergy;   
   }
   return *this;
}


