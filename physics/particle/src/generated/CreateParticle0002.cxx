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
void CreateParticle0002() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating ~u_R_bar
   new Particle("~u_R_bar", -2000002, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000002));
   daughters.clear();
   daughters.push_back(-1000039);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000037);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(-23);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(-25);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(-35);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(-36);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000001);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000001);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000021);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));

   // Creating ~d_R_bar
   new Particle("~d_R_bar", -2000001, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000001));
   daughters.clear();
   daughters.push_back(-1000039);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000024);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000037);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(-23);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(-25);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(-35);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(-36);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000002);
   daughters.push_back(24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000002);
   daughters.push_back(37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000021);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));

   // Creating ~chi_2-
   new Particle("~chi_2-", -1000037, 0, "Sparticle", 100, -1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000037));
   daughters.clear();
   daughters.push_back(-1000039);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000039);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-23);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-11);
   daughters.push_back(11);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-13);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-15);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-12);
   daughters.push_back(12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-14);
   daughters.push_back(14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-16);
   daughters.push_back(16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-1);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-3);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-5);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-2);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-4);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-25);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-35);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(-36);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000002);
   daughters.push_back(1);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000001);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000001);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000004);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000004);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000003);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000003);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000006);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000006);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000005);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000005);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000012);
   daughters.push_back(11);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000012);
   daughters.push_back(11);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000011);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000011);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000014);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000014);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000013);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000013);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000016);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000016);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000015);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000015);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000021);
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000021);
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
