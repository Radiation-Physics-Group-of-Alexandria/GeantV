//===----------------------------------------------------------------------===//
/**
 * @file TemplateGUVField.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

//
//
// class TemplateGUVField
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows any kind of field (vector, scalar, tensor and any set of them)
// to be defined by implementing the inquiry function interface.
//
// The key method is  GetFieldValue( const  double Point[4],
//                    *************         double *fieldArr ) 
// Given an input position/time vector 'Point', 
// this method must return the value of the field in "fieldArr".
//
// A field must also specify whether it changes a track's energy:
//                    DoesFieldChangeEnergy() 
//                    *********************
// A field must co-work with a corresponding Equation of Motion, to
// enable the integration of a particle's position, momentum and, optionally, 
// spin.  For this a field and its equation of motion must follow the
// same convention for the order of field components in the array "fieldArr"
// -------------------------------------------------------------------

#ifndef VECGUVFIELD_HH
#define VECGUVFIELD_HH

#include <vector>
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include "backend/Backend.h"
#include "AlignedBase.h"


class VecGUVField : public AlignedBase
{
  public: 

      inline
      VecGUVField( int NumberOfComponents, bool changesEnergy );
      inline
      VecGUVField( const VecGUVField &);
      virtual ~VecGUVField();

      //Vector interface with specialization
      template <class Backend>
      void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position,
                                vecgeom::Vector3D<typename Backend::precision_v> &FieldValue );

      bool DoesFieldChangeEnergy() const { return fChangesEnergy; } 
      int  GetNumberOfComponents() const { return fNumberOfComponents; } 

      VecGUVField& operator = (const VecGUVField &p); // Useful ?
      
      virtual VecGUVField* Clone() const;

  private:
      const int  fNumberOfComponents; 
      bool       fChangesEnergy; 
};


inline VecGUVField::VecGUVField( int numberOfComponents, bool changesEnergy )
   : fNumberOfComponents(numberOfComponents),
     fChangesEnergy(changesEnergy)
     //GUVField(numberOfComponents, changesEnergy)
{
  // std::cout<<"-- entered VecGUVField  constructor--"<<std::endl;
}


inline VecGUVField::VecGUVField( const VecGUVField &field) 
  :  fNumberOfComponents(field.fNumberOfComponents)
    //GUVField(field)
{
  fChangesEnergy= field.fChangesEnergy;
}


VecGUVField::~VecGUVField()
{
}


VecGUVField& VecGUVField::operator = (const VecGUVField &p)
{
  
   if (&p != this){

    //line below if inheriting from GUVField. Comment 2nd line in that case
    // this->GUVField::operator=(p);
   
     this->fChangesEnergy= p.fChangesEnergy;   
   }
   return *this;
}

VecGUVField* VecGUVField::Clone() const
{
    std::cout << "Derived class does not implement cloning,\n"
              << "but Clone method called.\n"
              << "Cannot continue;" << std::endl;
    exit(1); 
    return NULL;
}


#endif /* VECGUVFVECTORIELD_HH */
