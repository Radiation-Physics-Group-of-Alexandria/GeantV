
#include "Region.h"

namespace vecgeom {
  inline namespace VECGEOM_IMPL_NAMESPACE {

std::vector<Region*> Region::gTheRegionTable;

Region::Region(const std::string &name, bool iscutinlength, double gcut, double emcut, double epcut) : fName(name),
  fIsProductionCutInLength(iscutinlength), fGammaCutValue(gcut), fElectronCutValue(emcut), fPositronCutValue(epcut) {
  fIndex = gTheRegionTable.size();
  gTheRegionTable.push_back(this);
}

Region::Region(const std::string &name, bool iscutinlength) : fName(name), fIsProductionCutInLength(iscutinlength),
  fGammaCutValue(-1.0), fElectronCutValue(-1.0), fPositronCutValue(-1.0) {
  fIndex = gTheRegionTable.size();
  gTheRegionTable.push_back(this);
}

} // namespace vecgeom
} // namespace VECGEOM_IMPL_NAMESPACE
