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
void CreateParticle0017() {
   vector<int> daughters;
   ParticleOld *part = nullptr;

   // Creating B0_bar
   new ParticleOld("B0_bar", -511, 0, "B-Meson", 100, 0, 5.27953, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-511));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(4);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(48, 0.4291,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(4);
   daughters.push_back(1);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(13, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(-4);
   daughters.push_back(3);
   daughters.push_back(4);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(13, 0.07,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(413);
   part->AddDecay(ParticleOld::Decay(42, 0.055,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(413);
   part->AddDecay(ParticleOld::Decay(42, 0.055,  daughters));
   daughters.clear();
   daughters.push_back(-16);
   daughters.push_back(15);
   daughters.push_back(413);
   part->AddDecay(ParticleOld::Decay(42, 0.03,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-433);
   part->AddDecay(ParticleOld::Decay(0, 0.025,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(411);
   part->AddDecay(ParticleOld::Decay(42, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(411);
   part->AddDecay(ParticleOld::Decay(42, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(-4);
   daughters.push_back(4);
   daughters.push_back(3);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(13, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(-433);
   part->AddDecay(ParticleOld::Decay(0, 0.0185,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-20213);
   part->AddDecay(ParticleOld::Decay(0, 0.018,  daughters));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(-431);
   part->AddDecay(ParticleOld::Decay(0, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(2);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-431);
   part->AddDecay(ParticleOld::Decay(0, 0.0135,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(415);
   part->AddDecay(ParticleOld::Decay(42, 0.012,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(415);
   part->AddDecay(ParticleOld::Decay(42, 0.012,  daughters));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(-213);
   part->AddDecay(ParticleOld::Decay(0, 0.011,  daughters));
   daughters.clear();
   daughters.push_back(-16);
   daughters.push_back(15);
   daughters.push_back(411);
   part->AddDecay(ParticleOld::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-213);
   part->AddDecay(ParticleOld::Decay(0, 0.009,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(20413);
   part->AddDecay(ParticleOld::Decay(42, 0.008,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(20413);
   part->AddDecay(ParticleOld::Decay(42, 0.008,  daughters));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(-20213);
   part->AddDecay(ParticleOld::Decay(0, 0.0055,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(10411);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(10413);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(10411);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(10413);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-4);
   daughters.push_back(3);
   daughters.push_back(2);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.0042,  daughters));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.0035,  daughters));
   daughters.clear();
   daughters.push_back(-20443);
   daughters.push_back(-313);
   part->AddDecay(ParticleOld::Decay(0, 0.0025,  daughters));
   daughters.clear();
   daughters.push_back(-20443);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.0019,  daughters));
   daughters.clear();
   daughters.push_back(-443);
   daughters.push_back(-313);
   part->AddDecay(ParticleOld::Decay(0, 0.0014,  daughters));
   daughters.clear();
   daughters.push_back(-443);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.0008,  daughters));
   daughters.clear();
   daughters.push_back(-441);
   daughters.push_back(-313);
   part->AddDecay(ParticleOld::Decay(0, 0.0007,  daughters));
   daughters.clear();
   daughters.push_back(-441);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.0004,  daughters));

   // Creating D*_2s-
   new ParticleOld("D*_2s-", -435, 0, "CharmedMeson", 100, -1, 2.5726, 0.015, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-435));
   daughters.clear();
   daughters.push_back(-421);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.4,  daughters));
   daughters.clear();
   daughters.push_back(-411);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.4,  daughters));
   daughters.clear();
   daughters.push_back(-423);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.1,  daughters));
   daughters.clear();
   daughters.push_back(-413);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.1,  daughters));

   // Creating D*_s-
   new ParticleOld("D*_s-", -433, 0, "CharmedMeson", 100, -1, 2.1123, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-433));
   daughters.clear();
   daughters.push_back(-431);
   daughters.push_back(-22);
   part->AddDecay(ParticleOld::Decay(0, 0.94,  daughters));
   daughters.clear();
   daughters.push_back(-431);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.06,  daughters));

   // Creating D_s-
   new ParticleOld("D_s-", -431, 0, "CharmedMeson", 100, -1, 1.9685, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-431));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(-3);
   daughters.push_back(3);
   part->AddDecay(ParticleOld::Decay(13, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   part->AddDecay(ParticleOld::Decay(13, 0.0952,  daughters));
   daughters.clear();
   daughters.push_back(-331);
   daughters.push_back(-213);
   part->AddDecay(ParticleOld::Decay(0, 0.095,  daughters));
   daughters.clear();
   daughters.push_back(-221);
   daughters.push_back(-213);
   part->AddDecay(ParticleOld::Decay(0, 0.079,  daughters));
   daughters.clear();
   daughters.push_back(-333);
   daughters.push_back(-213);
   part->AddDecay(ParticleOld::Decay(0, 0.052,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(313);
   part->AddDecay(ParticleOld::Decay(0, 0.05,  daughters));
   daughters.clear();
   daughters.push_back(-331);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.037,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(311);
   part->AddDecay(ParticleOld::Decay(0, 0.033,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-333);
   part->AddDecay(ParticleOld::Decay(42, 0.03,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-333);
   part->AddDecay(ParticleOld::Decay(42, 0.03,  daughters));
   daughters.clear();
   daughters.push_back(-333);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.028,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(311);
   part->AddDecay(ParticleOld::Decay(0, 0.028,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(313);
   part->AddDecay(ParticleOld::Decay(0, 0.026,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-331);
   part->AddDecay(ParticleOld::Decay(42, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-221);
   part->AddDecay(ParticleOld::Decay(42, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-221);
   part->AddDecay(ParticleOld::Decay(42, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-331);
   part->AddDecay(ParticleOld::Decay(42, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(-221);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(ParticleOld::Decay(0, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(2112);
   part->AddDecay(ParticleOld::Decay(0, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(-10221);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.0078,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-321);
   daughters.push_back(321);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-311);
   daughters.push_back(311);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-221);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-331);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-333);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-221);
   daughters.push_back(-323);
   part->AddDecay(ParticleOld::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-311);
   daughters.push_back(311);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-321);
   daughters.push_back(321);
   part->AddDecay(ParticleOld::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-213);
   daughters.push_back(-113);
   part->AddDecay(ParticleOld::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-213);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-113);
   part->AddDecay(ParticleOld::Decay(0, 0.001,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
