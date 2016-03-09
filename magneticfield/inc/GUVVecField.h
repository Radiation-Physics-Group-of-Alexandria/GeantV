

#ifndef GUVVECFIELD_HH
#define GUVVECFIELD_HH

#include <vector>
#include "base/Vector3D.h"
#include "base/Global.h"
#include "backend/Backend.h"

#include "GUVField.h"
#include "VcFloatBackend.h"

class GUVVecField 
{
  public: 

      inline
      GUVVecField( int NumberOfComponents, bool changesEnergy );
      inline
      GUVVecField( const GUVVecField &);
      virtual ~GUVVecField(){};

      //Vector interface with specialization
      virtual void GetFieldValue( const vecgeom::Vector3D<typename vecgeom::kVc::precision_v>      &Position,
                                        vecgeom::Vector3D<typename vecgeom::kVcFloat::precision_v> &FieldValue ) = 0;

      bool DoesFieldChangeEnergy() const { return fChangesEnergy; } 
      int  GetNumberOfComponents() const { return fNumberOfComponents; } 

      GUVVecField& operator = (const GUVVecField &p);

      // virtual GUVVecField* Clone() const;

  private:
      const int  fNumberOfComponents; 
      bool       fChangesEnergy; 
};

inline GUVVecField::GUVVecField( int numberOfComponents, bool changesEnergy )
   : fNumberOfComponents(numberOfComponents),
     fChangesEnergy(changesEnergy)
{
  std::cout<<"-- entered GUVVecField  constructor--"<<std::endl;
}


inline GUVVecField::GUVVecField( const GUVVecField &field) 
  :  fNumberOfComponents(field.fNumberOfComponents)
{
  fChangesEnergy= field.fChangesEnergy;
}
#endif /* GUVVECFIELD_HH */
