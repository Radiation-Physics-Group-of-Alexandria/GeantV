#ifndef TNudyEndfFissionYield_H
#define TNudyEndfFissionYield_H

#include "TNudyEndfRecoPoint.h"

namespace Nudy {
	class TNudyEndfFile;
}

typedef std::vector<double> rowd;
typedef std::vector<rowd> matrixd2;
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

namespace NudyPhysics {
class TNudyEndfFissionYield : public NudyPhysics::TNudyEndfRecoPoint {

public:
  TNudyEndfFissionYield();
  TNudyEndfFissionYield(Nudy::TNudyEndfFile *file);
  virtual double GetFisYield(int elemid, double energyK);
  virtual ~TNudyEndfFissionYield();

private:
  double A, AWR, ABN, QX;                                      // standard ENDF parameters
  rowd ein, einc;                                              // incident energy
  matrixd2 zafp, fps, zafpc, fpsc, yi, cyi, dyi, yc, dyc;      // charge, mass, yield (independent and cummulative)
  rowd zafp1, fps1, zafpc1, fpsc1, yi1, cyi1, dyi1, yc1, dyc1; // charge, mass, yield (independent and cummulative)
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif
  ClassDef(TNudyEndfFissionYield, 1) // class for an ENDF fission yield reconstruction
};

} //namespace
#endif
