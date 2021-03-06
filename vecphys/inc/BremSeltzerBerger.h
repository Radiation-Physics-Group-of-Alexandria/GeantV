#ifndef BremSeltzerBerger_H
#define BremSeltzerBerger_H 1

#include "base/PhysicalConstants.h"
#include "base/VecPhys.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "Physics2DVector.h"
#include <fstream>

#include "EmModelBase.h"

namespace vecphys {

VECPHYS_DEVICE_FORWARD_DECLARE(class GUAliasSampler;)
inline namespace VECPHYS_IMPL_NAMESPACE {

class BremSeltzerBerger : public EmModelBase<BremSeltzerBerger> {
public:
  VECCORE_ATT_HOST
  BremSeltzerBerger(Random_t *states = 0, int threadId = -1);

  VECCORE_ATT_HOST_DEVICE
  BremSeltzerBerger(Random_t *states, int threadId, GUAliasSampler *sampler, Physics2DVector *sbData);

  VECCORE_ATT_HOST_DEVICE
  ~BremSeltzerBerger();

  VECCORE_ATT_HOST
  void Initialization();

  // interfaces for tables
  VECCORE_ATT_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_ATT_HOST
  void BuildPdfTable(int Z, double *p);

public:
  // Auxiliary methods
  VECCORE_ATT_HOST_DEVICE
  Physics2DVector *GetSBData() { return fDataSB; }

  VECCORE_ATT_HOST bool RetrieveSeltzerBergerData(std::ifstream &in, Physics2DVector *vec2D);

  // Implementation methods
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v CrossSectionKernel(typename Backend::Double_v energyIn,
                                                                        Index_v<typename Backend::Double_v> zElement);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernel(typename Backend::Double_v energyIn,
                                              Index_v<typename Backend::Double_v> zElement,
                                              typename Backend::Double_v &energyOut,
                                              typename Backend::Double_v &sinTheta);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernelCR(typename Backend::Double_v energyIn,
                                                Index_v<typename Backend::Double_v> zElement,
                                                typename Backend::Double_v &energyOut,
                                                typename Backend::Double_v &sinTheta);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernelUnpack(typename Backend::Double_v energyIn,
                                                    Index_v<typename Backend::Double_v> zElement,
                                                    typename Backend::Double_v &energyOut,
                                                    typename Backend::Double_v &sinTheta,
                                                    Mask_v<typename Backend::Double_v> &status);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v SampleSinTheta(typename Backend::Double_v energyIn);

  VECCORE_ATT_HOST_DEVICE
  void SampleByCompositionRejection(int elementZ, double energyIn, double &energyOut, double &sinTheta);

  VECCORE_ATT_HOST_DEVICE
  double CalculateDiffCrossSection(int Zelement, double Ein, double Eout) const;

  // the cross section calculation from Geant4

  VECCORE_ATT_HOST double GetG4CrossSection(int Z, double energyIn);

  VECCORE_ATT_HOST_DEVICE
  void SetCurrentElement(double Z);

  VECCORE_ATT_HOST_DEVICE double ComputeXSectionPerAtom(double cut, double energy);

  VECCORE_ATT_HOST_DEVICE
  double ComputeRelDXSectionPerAtom(double gammaEnergy);

  VECCORE_ATT_HOST_DEVICE
  double ComputeDXSectionPerAtom(double gammaEnergy);

  VECCORE_ATT_HOST_DEVICE
  void CalcLPMFunctions(double k);

  VECCORE_ATT_HOST_DEVICE
  double Phi1(double gg);

  VECCORE_ATT_HOST_DEVICE
  double Phi1M2(double gg);

  VECCORE_ATT_HOST_DEVICE
  double Psi1(double eps);

  VECCORE_ATT_HOST_DEVICE
  double Psi1M2(double eps);

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<BremSeltzerBerger>;

private:
  Physics2DVector *fDataSB;

  // Geant4 cross section
private:
  double totalEnergy;
  double currentZ;
  double densityFactor;
  double densityCorr;

  double lpmEnergy;
  double xiLPM;
  double phiLPM;
  double gLPM;

  double z13;
  double z23;
  double lnZ;
  double Fel;
  double Finel;
  double fCoulomb;
  double fMax;
};

// Implementation
template <class Backend>
VECCORE_ATT_HOST_DEVICE
    typename Backend::Double_v BremSeltzerBerger::CrossSectionKernel(typename Backend::Double_v /*energy*/,
                                                                     Index_v<typename Backend::Double_v> /*Z*/)
{
  return 1.0;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void BremSeltzerBerger::InteractKernel(typename Backend::Double_v energyIn,
                                                               Index_v<typename Backend::Double_v> zElement,
                                                               typename Backend::Double_v &energyOut,
                                                               typename Backend::Double_v &sinTheta)
{
  using Double_v = typename Backend::Double_v;

  Index_v<Double_v> irow;
  Index_v<Double_v> icol;
  Double_v fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn, irow, icol, fraction);

  Double_v probNA;
  Index_v<Double_v> aliasInd;

  // this did not used to work - Fixed SW
  Index_v<Double_v> ncol(fAliasSampler->GetSamplesPerEntry());
  Index_v<Double_v> index = ncol * irow + icol;
  fAliasSampler->GatherAlias<Backend>(index, zElement, probNA, aliasInd);

  // Seltzer-Berger specific

  // To-do: apply densityFactor (dummy for now) and calculate deltaY correctly
  // densityFactor = (Migdal constant)x(electron density of the material);
  Double_v densityFactor = 1.0;

  Double_v emin = math::Min(Double_v(fAliasSampler->GetIncomingMin()), energyIn);
  Double_v emax = math::Min(Double_v(fAliasSampler->GetIncomingMax()), energyIn);

  Double_v totalEnergy = energyIn + electron_mass_c2;
  Double_v densityCorr = densityFactor * totalEnergy * totalEnergy;
  Double_v minY        = math::Log(emin * emin + densityCorr);
  Double_v maxY        = math::Log(emax * emax + densityCorr);
  Double_v deltaY      = maxY - minY;

  Double_v yhat = fAliasSampler->SampleX<Backend>(deltaY, probNA, aliasInd, icol, fraction);

  energyOut = math::Sqrt(math::Max(math::Exp(minY + yhat) - densityCorr, Double_v(0.0)));
  sinTheta  = SampleSinTheta<Backend>(energyOut);
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v BremSeltzerBerger::SampleSinTheta(
    typename Backend::Double_v energyIn)
{
  using Double_v = typename Backend::Double_v;

  // angle of the radiated photon
  // based on G4DipBustGenerator::SampleDirection

  Double_v c     = Double_v(4.0) - Double_v(8.0) * UniformRandom<Double_v>(fRandomState, fThreadId);
  Double_v signc = math::Sign(c);
  Double_v a     = math::Abs(c);

  Double_v delta = 0.5 * (math::Sqrt(4.0 + a * a) + a);

  Double_v cofA = -signc * math::Pow(delta, Double_v(1.0 / 3.0));

  Double_v cosTheta = cofA - 1.0 / cofA;

  Double_v tau  = energyIn / electron_mass_c2;
  Double_v beta = math::Sqrt(tau * (tau + 2.0)) / (tau + 1.0);

  cosTheta = (cosTheta + beta) / (1.0 + cosTheta * beta);

  Double_v sinTheta = math::Sqrt((1.0 - cosTheta) * (1.0 + cosTheta));

  return sinTheta;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void BremSeltzerBerger::InteractKernelCR(typename Backend::Double_v energyIn,
                                                                 Index_v<typename Backend::Double_v> /*Z*/,
                                                                 typename Backend::Double_v &energyOut,
                                                                 typename Backend::Double_v &sinTheta)
{
  // dummy for now
  energyOut = energyIn;
  sinTheta  = 0.0;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void BremSeltzerBerger::InteractKernelUnpack(typename Backend::Double_v energyIn,
                                                                     Index_v<typename Backend::Double_v> /*Z*/,
                                                                     typename Backend::Double_v &energyOut,
                                                                     typename Backend::Double_v &sinTheta,
                                                                     Mask_v<typename Backend::Double_v> &status)
{
  // dummy for now
  energyOut = energyIn;
  sinTheta  = 0.0;
  status    = false;
}

} // end namespace impl
} // end namespace vecphys

#endif
