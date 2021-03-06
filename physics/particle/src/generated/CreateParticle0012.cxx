// This files was autogenerated by geant::ParticleOld::ReadFile

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "ParticleOld.h"
#ifdef VECCORE_CUDA
#include "base/Vector.h"
template <typename T>
using vector = vecgeom::Vector<T>;
#else
using std::vector;
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void CreateParticle0012() {
   vector<int> daughters;
   ParticleOld *part = nullptr;

   // Creating Xi_b+
   new ParticleOld("Xi_b+", -5132, 0, "B-Baryon", 100, 1, 5.7924, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-5132));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));

   // Creating Lambda_b0_bar
   new ParticleOld("Lambda_b0_bar", -5122, 0, "B-Baryon", 100, 0, 5.6202, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-5122));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-2101);
   part->AddDecay(ParticleOld::Decay(48, 0.4291,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4122);
   part->AddDecay(ParticleOld::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4122);
   part->AddDecay(ParticleOld::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-2101);
   part->AddDecay(ParticleOld::Decay(13, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-2101);
   part->AddDecay(ParticleOld::Decay(13, 0.07,  daughters));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(433);
   part->AddDecay(ParticleOld::Decay(0, 0.0435,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4122);
   part->AddDecay(ParticleOld::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(431);
   part->AddDecay(ParticleOld::Decay(0, 0.0285,  daughters));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(20213);
   part->AddDecay(ParticleOld::Decay(0, 0.0235,  daughters));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(213);
   part->AddDecay(ParticleOld::Decay(0, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-2101);
   part->AddDecay(ParticleOld::Decay(13, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-2101);
   part->AddDecay(ParticleOld::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.0077,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-2101);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-20443);
   daughters.push_back(-3122);
   part->AddDecay(ParticleOld::Decay(0, 0.0044,  daughters));
   daughters.clear();
   daughters.push_back(-443);
   daughters.push_back(-3122);
   part->AddDecay(ParticleOld::Decay(0, 0.0022,  daughters));
   daughters.clear();
   daughters.push_back(-441);
   daughters.push_back(-3122);
   part->AddDecay(ParticleOld::Decay(0, 0.0011,  daughters));

   // Creating Sigma*_b-_bar
   new ParticleOld("Sigma*_b-_bar", -5114, 0, "B-Baryon", 100, 1, 5.8364, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-5114));
   daughters.clear();
   daughters.push_back(-5122);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 1,  daughters));

   // Creating Sigma_b-_bar
   new ParticleOld("Sigma_b-_bar", -5112, 0, "B-Baryon", 100, 1, 5.8152, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-5112));
   daughters.clear();
   daughters.push_back(-5122);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 1,  daughters));

   // Creating bd_1_bar
   new ParticleOld("bd_1_bar", -5103, 0, "Unknown", 100, 0.666667, 5.40145, 0, 100, 100, 1, 100, 1);

   // Creating bd_0_bar
   new ParticleOld("bd_0_bar", -5101, 0, "Unknown", 100, 0.666667, 5.38897, 0, 100, 100, 1, 100, 1);

   // Creating Omega*_ccc--
   new ParticleOld("Omega*_ccc--", -4444, 0, "CharmedBaryon", 100, -2, 4.91594, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-4444));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.76,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(3);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.08,  daughters));

   // Creating Omega*_cc-
   new ParticleOld("Omega*_cc-", -4434, 0, "CharmedBaryon", 100, -1, 3.82466, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-4434));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.76,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(3);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.08,  daughters));

   // Creating Omega_cc-
   new ParticleOld("Omega_cc-", -4432, 0, "CharmedBaryon", 100, -1, 3.78663, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-4432));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.76,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(3);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.08,  daughters));

   // Creating Xi*_cc--
   new ParticleOld("Xi*_cc--", -4424, 0, "CharmedBaryon", 100, -2, 3.65648, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-4424));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.76,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(3);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(ParticleOld::Decay(11, 0.08,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
