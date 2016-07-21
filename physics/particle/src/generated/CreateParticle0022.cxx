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
void CreateParticle0022() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating R0_bar
   new Particle("R0_bar", -40, 0, "Unknown", 100, 0, 5000, 417.465, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-40));
   daughters.clear();
   daughters.push_back(-1);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(32, 0.215134,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(32, 0.215134,  daughters));
   daughters.clear();
   daughters.push_back(-3);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(32, 0.215133,  daughters));
   daughters.clear();
   daughters.push_back(-4);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(32, 0.214738,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(0, 0.0699301,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.0699301,  daughters));
   daughters.clear();
   daughters.push_back(-5);
   daughters.push_back(7);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-6);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-15);
   daughters.push_back(17);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating LQ_ue_bar
   new Particle("LQ_ue_bar", -39, 0, "Unknown", 100, 0.333333, 200, 0.39162, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-39));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating H-
   new Particle("H-", -37, 0, "GaugeBoson", 100, -1, 300, 5.75967, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-37));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-25);
   part->AddDecay(Particle::Decay(0, 0.929792,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.067484,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(0, 0.002701,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 1.3e-05,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(0, 1e-05,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(17);
   daughters.push_back(-18);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(-1000024);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000022);
   daughters.push_back(-1000037);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(-1000024);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000023);
   daughters.push_back(-1000037);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(-1000024);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000025);
   daughters.push_back(-1000037);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(-1000024);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000035);
   daughters.push_back(-1000037);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000006);
   daughters.push_back(1000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000006);
   daughters.push_back(1000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000006);
   daughters.push_back(2000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-2000006);
   daughters.push_back(2000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000001);
   daughters.push_back(-1000002);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000003);
   daughters.push_back(-1000004);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000011);
   daughters.push_back(-1000012);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000013);
   daughters.push_back(-1000014);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000015);
   daughters.push_back(-1000016);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000015);
   daughters.push_back(-1000016);
   part->AddDecay(Particle::Decay(53, 0,  daughters));

   // Creating W'-
   new Particle("W'-", -34, 0, "GaugeBoson", 100, -1, 500, 16.6708, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-34));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.251276,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.250816,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.215459,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(0, 0.085262,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(0, 0.085262,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(0, 0.08526,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.012903,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.012903,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.000465,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.00038,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 8e-06,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 6e-06,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(17);
   daughters.push_back(-18);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-23);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-25);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating W-
   new Particle("W-", -24, 0, "GaugeBoson", 100, -1, 80.398, 2.07002, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-24));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.321502,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.320778,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(0, 0.108062,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(0, 0.108062,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(0, 0.107983,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.016509,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.016502,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.000591001,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 1e-05,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(17);
   daughters.push_back(-18);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif