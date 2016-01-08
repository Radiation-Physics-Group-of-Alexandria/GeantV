//===----------------------------------------------------------------------===//
/**
 * @file GUVField.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

//
//
// class GUVField
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

#ifndef GUVFIELD_HH
#define GUVFIELD_HH

// #include "GUVTypes.hh"
// #include "globals.hh"

/**
 * @brief Class GUVField
 */

class GUVField
{
  public:  // with description

      /**
       * @brief GeantTrack parametrized constructor
       *
       * @param Point    - position (0,1,2=x,y,z) and time (3) [Input]
       * @param fieldArr - output values of field. Usual convention:
       *                   0,1,2 = B_x, B_y, B_z
       *                   3,4,5 = E_x, E_y, E_z  (foreseen extension)
       *        Units are expected to be native GeantV units.
       */
      virtual void  GetFieldValue( const double Point[4],
                                         double *fieldArr ) const = 0;
      inline
      GUVField( int NumberOfComponents, bool changesEnergy );
      inline
      GUVField( const GUVField & );
      virtual ~GUVField();
      // inline GUVField& operator = (const GUVField &p); 

     // A field signature function that can be used to insure
     // that the Equation of motion object and the GUVField object
     // have the same "field signature"?

      bool   DoesFieldChangeEnergy() const { return fChangesEnergy; } 
      int    GetNumberOfComponents() const { return fNumberOfComponents; } 

      GUVField& operator = (const GUVField &p); // Useful ?
      
      virtual GUVField* Clone() const;
        // Implements cloning, likely needed for MT 

private:
      const int  fNumberOfComponents; // E.g.  B -> N=3 , ie x,y,z 
                                     //       E+B -> N=6 
      bool fChangesEnergy; 
       //  Each type/class of field set this accordingly:
       //    - an electric field     - "true"
       //    - a pure magnetic field - "false"
};

inline GUVField::GUVField( int NumberOfComponents, bool changesEnergy )
   : fNumberOfComponents(NumberOfComponents),
     fChangesEnergy(changesEnergy)
{
}

inline GUVField::GUVField( const GUVField &field )
   : fNumberOfComponents(field.fNumberOfComponents)
{
   // *this = field; 
   fChangesEnergy= field.fChangesEnergy;
}
#endif /* GUVFIELD_HH */
