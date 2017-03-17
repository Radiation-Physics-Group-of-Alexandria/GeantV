
#ifndef PHYSICSDATA_H
#define PHYSICSDATA_H

#include <vector>

namespace geantphysics {

class LightTrack;

class PhysicsData {
public:
  PhysicsData();
 ~PhysicsData(){}

  int  GetNumUsedSecondaries () const  { return fNumUsedSecondaries; }
  void SetNumUsedSecondaries (int val) { fNumUsedSecondaries = val;  }

  int  GetSizeListOfSecondaries() const  { return fListOfSecondaries.size(); }
  void SetSizeListOfSecondaries(int val) { fListOfSecondaries.resize(val); }

  //  Resize list if needed to accomodate the extra secondaries
  void PrepareForSecondaries( int numNewSecondaries );

  std::vector<LightTrack>& GetListOfSecondaries() { return fListOfSecondaries; }

  static void ClearAll();
  static std::vector<PhysicsData*> gThePhysicsDataTable;

private:
  int    fNumUsedSecondaries;   // number of secondary tracks currently used from fListOfSecondaries
  std::vector<LightTrack>  fListOfSecondaries;
};

inline
void PhysicsData::PrepareForSecondaries( int newSecondaries )
{
   int  currentSize= fListOfSecondaries.size();
   if(  fNumUsedSecondaries + newSecondaries > currentSize  ) {
      SetSizeListOfSecondaries( 2 * currentSize );
   }
}


} // namespace geantphysics

#endif
