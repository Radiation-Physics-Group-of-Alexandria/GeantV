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
void CreateParticle0020() {
   vector<int> daughters;
   ParticleOld *part = nullptr;

   // Creating phi3(1850)_bar
   new ParticleOld("phi3(1850)_bar", -337, 0, "Unknown", 100, 0, 1.854, 0.087, 100, 100, 0, 100, 1);

   // Creating k3_star(1780)-_bar
   new ParticleOld("k3_star(1780)-_bar", -327, 0, "Unknown", 100, -1, 1.776, 0.159, 100, 100, 0, 100, 1);

   // Creating K*_2-
   new ParticleOld("K*_2-", -325, 0, "Meson", 100, -1, 1.4256, 0.098, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-325));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.332,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.168,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.166,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(-211);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.086,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.084,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-213);
   part->AddDecay(ParticleOld::Decay(0, 0.059,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(-211);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-113);
   part->AddDecay(ParticleOld::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-223);
   part->AddDecay(ParticleOld::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-221);
   part->AddDecay(ParticleOld::Decay(0, 0.002,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-22);
   part->AddDecay(ParticleOld::Decay(0, 0.002,  daughters));

   // Creating K*-
   new ParticleOld("K*-", -323, 0, "Meson", 100, -1, 0.89166, 0.0498, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-323));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(3, 0.666,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(3, 0.333,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-22);
   part->AddDecay(ParticleOld::Decay(0, 0.001,  daughters));

   // Creating K-
   new ParticleOld("K-", -321, 0, "Meson", 100, -1, 0.493677, 5.31674e-17, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-321));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(ParticleOld::Decay(0, 0.6352,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.2116,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-211);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.0559,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(42, 0.0482,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(42, 0.0318,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-111);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.0173,  daughters));

   // Creating k3_star(1780)0_bar
   new ParticleOld("k3_star(1780)0_bar", -317, 0, "Unknown", 100, 0, 1.776, 0.159, 100, 100, 0, 100, 1);

   // Creating K*_20_bar
   new ParticleOld("K*_20_bar", -315, 0, "Meson", 100, 0, 1.4324, 0.109, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-315));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.168,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.166,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(211);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.087,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.084,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(213);
   part->AddDecay(ParticleOld::Decay(0, 0.059,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(-211);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-113);
   part->AddDecay(ParticleOld::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-223);
   part->AddDecay(ParticleOld::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-221);
   part->AddDecay(ParticleOld::Decay(0, 0.002,  daughters));

   // Creating K*0_bar
   new ParticleOld("K*0_bar", -313, 0, "Meson", 100, 0, 0.896, 0.0505, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-313));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(3, 0.665,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(3, 0.333,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-22);
   part->AddDecay(ParticleOld::Decay(0, 0.002,  daughters));

   // Creating K0_bar
   new ParticleOld("K0_bar", -311, 0, "Meson", 100, 0, 0.497614, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-311));
   daughters.clear();
   daughters.push_back(-130);
   part->AddDecay(ParticleOld::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(-310);
   part->AddDecay(ParticleOld::Decay(0, 0.5,  daughters));

   // Creating omega3(1670)_bar
   new ParticleOld("omega3(1670)_bar", -227, 0, "Unknown", 100, 0, 1.667, 0.168, 100, 100, 0, 100, 1);

   // Creating rho3(1690)-_bar
   new ParticleOld("rho3(1690)-_bar", -217, 0, "Unknown", 100, -1, 1.6888, 0.161, 100, 100, 0, 100, 1);

   // Creating a_2-
   new ParticleOld("a_2-", -215, 0, "Meson", 100, -1, 1.3183, 0.107, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-215));
   daughters.clear();
   daughters.push_back(-213);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.34725,  daughters));
   daughters.clear();
   daughters.push_back(-113);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.34725,  daughters));
   daughters.clear();
   daughters.push_back(-221);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.144,  daughters));
   daughters.clear();
   daughters.push_back(-223);
   daughters.push_back(-211);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(0, 0.104,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(311);
   part->AddDecay(ParticleOld::Decay(0, 0.049,  daughters));
   daughters.clear();
   daughters.push_back(-331);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.0057,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-22);
   part->AddDecay(ParticleOld::Decay(0, 0.0028,  daughters));

   // Creating rho-
   new ParticleOld("rho-", -213, 0, "Meson", 100, -1, 0.77549, 0.149, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-213));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-111);
   part->AddDecay(ParticleOld::Decay(3, 0.99955,  daughters));
   daughters.clear();
   daughters.push_back(-211);
   daughters.push_back(-22);
   part->AddDecay(ParticleOld::Decay(0, 0.00045,  daughters));

   // Creating pi-
   new ParticleOld("pi-", -211, 0, "Meson", 100, -1, 0.13957, 2.52837e-17, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-211));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(ParticleOld::Decay(0, 0.999877,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(ParticleOld::Decay(0, 0.000123,  daughters));

   // Creating pi_diffr-
   new ParticleOld("pi_diffr-", -210, 0, "Meson", 100, -1, 0, 0, 100, 100, 1, 100, 1);

   // Creating rho3(1690)0_bar
   new ParticleOld("rho3(1690)0_bar", -117, 0, "Unknown", 100, 0, 1.6888, 0.161, 100, 100, 0, 100, 1);

   // Creating b-hadron_bar
   new ParticleOld("b-hadron_bar", -85, 0, "Generator", 100, 0.333333, 5, 0, 100, 100, 1, 100, 1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(-85));
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
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
