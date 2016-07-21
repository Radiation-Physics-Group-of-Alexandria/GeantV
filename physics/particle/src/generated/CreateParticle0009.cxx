// This files was autogenerated by geant::Particle::ReadFile

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
#ifdef GEANT_NVCC
#include "base/Vector.h"
template <typename T>
using vector = vecgeom::Vector<T>;
#else
using std::vector;
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void CreateParticle0009() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating Xi_bb+
   new Particle("Xi_bb+", -5512, 0, "Unknown", 100, 1, 10.4227, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5512));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));

   // Creating bb_1_bar
   new Particle("bb_1_bar", -5503, 0, "Unknown", 100, 0.666667, 10.0735, 0, 100, 100, 1, 100, 1);

   // Creating Omega*_bcc-
   new Particle("Omega*_bcc-", -5444, 0, "B-Baryon", 100, -1, 8.31325, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5444));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));

   // Creating Omega_bcc+_bar
   new Particle("Omega_bcc+_bar", -5442, 0, "B-Baryon", 100, -1, 8.30945, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5442));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));

   // Creating Omega*_bc0_bar
   new Particle("Omega*_bc0_bar", -5434, 0, "B-Baryon", 100, 0, 7.219, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5434));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));

   // Creating Omega'_bc0_bar
   new Particle("Omega'_bc0_bar", -5432, 0, "B-Baryon", 100, 0, 7.21101, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5432));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif