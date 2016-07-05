#ifndef _CMS_MAGNETIC_FIELD_H_
#define _CMS_MAGNETIC_FIELD_H_

#include "G4MagneticField.hh"

#include <vector>
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4ThreeVector.hh"

typedef G4ThreeVector Vector3D;

class CMSMagneticFieldG4 : public G4MagneticField
{
  public:
    CMSMagneticFieldG4();
    CMSMagneticFieldG4(std::string inputMap);
   
    void GetFieldValue( const G4double  Point[4],
                              G4double *Bfield ) const override final;

  public:
    //New stuff 
    //Takes as input x,y,z; Gives output Bx,By,Bz
    void GetFieldValueXYZ(const Vector3D &position, Vector3D &xyzField) const;

    //Stores rz field as well for cross checking purpose 
    void GetFieldValueTest(const Vector3D &position, Vector3D &rzField);

    //Takes as input an SOA3D for position, gives field
    //   void GetFieldValues(const vecgeom::SOA3D<double> &position, vecgeom::SOA3D<double> &field); //not tested yet with given input

    //Reads data from given 2D magnetic field map. Can be easily modified to read a given 2D map, in case the file changes
    void ReadVectorData(std::string inputMap);

    ~CMSMagneticFieldG4();

  private:
    //  Invariants -- parameters of the field 
    const double kRMax=  9000. * CLHEP::millimeter;   //  Maximum value of R =  9.00 meters
    const double kZMax= 16000. * CLHEP::millimeter;   //  Max value of Z = 16.00 meters
    const int kNoZValues = 161;
    const int kNoRValues = 181;
    const int kHalfZValues = 80;

    // Derived values
    //kRDiff and kZDiff take care of mm because they come from kRMax and kZMax which have mm in them
    const double kRDiff = kRMax   / (kNoRValues-1);  //  Radius increment between lattice points
    const double kZDiff = 2*kZMax / (kNoZValues-1);  //  Z increment

    const double kZ0 = -kZMax;
    const double kRDiffInv = 1.0/kRDiff;
    const double kZDiffInv = 1.0/kZDiff;
    const double kAInverse = 1/(kRDiff*kZDiff);

    //For (R,Z) pairs : obtain estimated field in cylindrical coordinates - output in rzfield
    void GetFieldValueRZ(const double radius,const double z, Vector3D &rzField) const; 
    void GetFieldValueRZ(const std::vector<double>   & radius,
                         const std::vector<double>   & z,
                               std::vector<Vector3D> & rzField ) const;
   
    //Used to convert cartesian coordinates to cylindrical coordinates R-Z-phi
    //  Does not calculate phi
    inline void CartesianToCylindrical(const Vector3D &cart, double cyl[2]) const;

    //Converts cylindrical magnetic field to field in cartesian coordinates 
    inline void CylindricalToCartesian(const Vector3D &B_cyl, const double sinTheta, const double cosTheta, Vector3D &B_cart) const;

private:
    std::vector<double> fRadius, fPhi, fZ, fBr, fBz, fBphi;
    bool  fInitialised;
};

inline void CMSMagneticFieldG4::CartesianToCylindrical(const Vector3D &cart, double cyl[2]) const
{
    //cyl[3] =[r,z,phi]
    cyl[0] = sqrt(cart[0]*cart[0] + cart[1]*cart[1]); // r = sqrt(x^2 + y^2)
    cyl[1] = cart[2]; //z = z 
}

inline void CMSMagneticFieldG4::CylindricalToCartesian(const Vector3D &rzField,
                                                     const double    sinTheta,
                                                     const double    cosTheta,
                                                           Vector3D &xyzField) const
{
    //B_cyl[] has r, phi and z
    xyzField[0] = rzField[0]*cosTheta - rzField[1]*sinTheta; // Bx= Br cos(theta) - Bphi sin(theta)
    xyzField[1] = rzField[0]*sinTheta + rzField[1]*cosTheta; //By = Br sin(theta) + Bphi cos(theta)
    xyzField[2] = rzField[2];   //Bz = Bz 
}

inline void
CMSMagneticFieldG4::GetFieldValue( const G4double  Point[4],
                                       G4double *Bfield ) const
{
   G4ThreeVector Position3v( Point[0], Point[1], Point[2] );
   G4ThreeVector Bfield3v( 0.0, 0.0, 0.0 );
   GetFieldValueXYZ( Position3v, Bfield3v );
   Bfield[0]= Bfield3v.x();
   Bfield[1]= Bfield3v.y();
   Bfield[2]= Bfield3v.z();   
}
#endif
