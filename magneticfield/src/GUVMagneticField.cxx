//
// First implementation class for GUVField
//
// Used to enable creation of vtable.
// -------------------------------------------------------------------

#include "GUVMagneticField.h"
#include <iostream>
 
GUVMagneticField::~GUVMagneticField()
{
}

void
GUVMagneticField::GetFieldValue( const double  Point[4],     // The old interface
                                 double* FieldArr )
{
   vecgeom::Vector3D<double> PositionV3D( Point[0], Point[1], Point[2]);
   vecgeom::Vector3D<float>  Field_v3f;
   this->GetFieldValue( PositionV3D, Field_v3f );
   FieldArr[0]= Field_v3f.x();
   FieldArr[1]= Field_v3f.y();
   FieldArr[2]= Field_v3f.z();
}
