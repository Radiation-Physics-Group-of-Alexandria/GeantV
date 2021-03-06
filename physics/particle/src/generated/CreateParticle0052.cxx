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
void CreateParticle0052() {
   vector<int> daughters;
   ParticleOld *part = nullptr;

   // Creating ~e_L-
   new ParticleOld("~e_L-", 1000011, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(1000011));
   daughters.clear();
   daughters.push_back(1000039);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000037);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000012);
   daughters.push_back(-24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000012);
   daughters.push_back(-24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000012);
   daughters.push_back(-37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000012);
   daughters.push_back(-37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));

   // Creating ~nu_eL
   new ParticleOld("~nu_eL", 1000012, 1, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(1000012));
   daughters.clear();
   daughters.push_back(1000039);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000024);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000037);
   daughters.push_back(11);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(12);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000011);
   daughters.push_back(24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000011);
   daughters.push_back(24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000011);
   daughters.push_back(37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000011);
   daughters.push_back(37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));

   // Creating ~mu_L-
   new ParticleOld("~mu_L-", 1000013, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(1000013));
   daughters.clear();
   daughters.push_back(1000039);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000037);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000014);
   daughters.push_back(-24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000014);
   daughters.push_back(-24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000014);
   daughters.push_back(-37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000014);
   daughters.push_back(-37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));

   // Creating ~nu_muL
   new ParticleOld("~nu_muL", 1000014, 1, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(1000014));
   daughters.clear();
   daughters.push_back(1000039);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000024);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000037);
   daughters.push_back(13);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(14);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000013);
   daughters.push_back(24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000013);
   daughters.push_back(24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000013);
   daughters.push_back(37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000013);
   daughters.push_back(37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));

   // Creating ~tau_1-
   new ParticleOld("~tau_1-", 1000015, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(1000015));
   daughters.clear();
   daughters.push_back(1000039);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000024);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000037);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000016);
   daughters.push_back(-24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000016);
   daughters.push_back(-24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000016);
   daughters.push_back(-37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000016);
   daughters.push_back(-37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));

   // Creating ~nu_tauL
   new ParticleOld("~nu_tauL", 1000016, 1, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(1000016));
   daughters.clear();
   daughters.push_back(1000039);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000024);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000037);
   daughters.push_back(15);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(16);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000015);
   daughters.push_back(24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000015);
   daughters.push_back(24);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000015);
   daughters.push_back(37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000015);
   daughters.push_back(37);
   part->AddDecay(ParticleOld::Decay(53, 0,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
