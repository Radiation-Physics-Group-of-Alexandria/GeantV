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
void CreateParticle0025() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating s_bar
   new Particle("s_bar", -3, 0, "Quark", 100, 0.333333, 0.104, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3));
   daughters.clear();
   daughters.push_back(-21);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-22);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-23);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-25);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating u_bar
   new Particle("u_bar", -2, 0, "Quark", 100, -0.666667, 0.0024, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2));
   daughters.clear();
   daughters.push_back(-21);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-22);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-23);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-7);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-25);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating d_bar
   new Particle("d_bar", -1, 0, "Quark", 100, 0.333333, 0.0048, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1));
   daughters.clear();
   daughters.push_back(-21);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-22);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-23);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-25);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating Rootino
   new Particle("Rootino", 0, 0, "Unknown", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating d
   new Particle("d", 1, 1, "Quark", 100, -0.333333, 0.0048, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating u
   new Particle("u", 2, 1, "Quark", 100, 0.666667, 0.0024, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(7);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating s
   new Particle("s", 3, 1, "Quark", 100, -0.333333, 0.104, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating c
   new Particle("c", 4, 1, "Quark", 100, 0.666667, 1.27, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(7);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));

   // Creating b
   new Particle("b", 5, 1, "Quark", 100, -0.333333, 4.68, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(102, 0,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
