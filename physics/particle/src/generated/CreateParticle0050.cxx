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
void CreateParticle0050() {
   vector<int> daughters;
   ParticleOld *part = nullptr;

   // Creating K*_10
   new ParticleOld("K*_10", 20313, 1, "Unknown", 100, 0, 1.403, 0.174, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20313));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));

   // Creating K*_1+
   new ParticleOld("K*_1+", 20323, 1, "Unknown", 100, 1, 1.403, 0.174, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20323));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));

   // Creating f'_1
   new ParticleOld("f'_1", 20333, 0, "Unknown", 100, 0, 1.4264, 0.053, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20333));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(311);
   part->AddDecay(ParticleOld::Decay(0, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(321);
   part->AddDecay(ParticleOld::Decay(0, 0.25,  daughters));

   // Creating D*_1+
   new ParticleOld("D*_1+", 20413, 1, "Unknown", 100, 1, 2.372, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20413));
   daughters.clear();
   daughters.push_back(423);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));

   // Creating D*_10
   new ParticleOld("D*_10", 20423, 1, "Unknown", 100, 0, 2.372, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20423));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(423);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));

   // Creating D*_1s+
   new ParticleOld("D*_1s+", 20433, 1, "Unknown", 100, 1, 2.4596, 0.0055, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20433));
   daughters.clear();
   daughters.push_back(433);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.8,  daughters));
   daughters.clear();
   daughters.push_back(433);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.2,  daughters));

   // Creating chi_1c
   new ParticleOld("chi_1c", 20443, 0, "Unknown", 100, 0, 3.51066, 0.0009, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20443));
   daughters.clear();
   daughters.push_back(82);
   daughters.push_back(-82);
   part->AddDecay(ParticleOld::Decay(12, 0.727,  daughters));
   daughters.clear();
   daughters.push_back(443);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.273,  daughters));

   // Creating B*_10
   new ParticleOld("B*_10", 20513, 1, "Unknown", 100, 0, 5.78, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20513));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));

   // Creating B*_1+
   new ParticleOld("B*_1+", 20523, 1, "Unknown", 100, 1, 5.78, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20523));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(211);
   part->AddDecay(ParticleOld::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.333,  daughters));

   // Creating B*_1s0
   new ParticleOld("B*_1s0", 20533, 1, "Unknown", 100, 0, 6.02, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20533));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(-321);
   part->AddDecay(ParticleOld::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(-311);
   part->AddDecay(ParticleOld::Decay(0, 0.5,  daughters));

   // Creating B*_1c+
   new ParticleOld("B*_1c+", 20543, 1, "Unknown", 100, 1, 7.3, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20543));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(411);
   part->AddDecay(ParticleOld::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(421);
   part->AddDecay(ParticleOld::Decay(0, 0.5,  daughters));

   // Creating chi_1b
   new ParticleOld("chi_1b", 20553, 0, "Unknown", 100, 0, 9.8928, 0, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(20553));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(21);
   part->AddDecay(ParticleOld::Decay(32, 0.65,  daughters));
   daughters.clear();
   daughters.push_back(553);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.35,  daughters));

   // Creating delta(1910)-
   new ParticleOld("delta(1910)-", 21112, 1, "Unknown", 100, -1, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating delta(1920)-
   new ParticleOld("delta(1920)-", 21114, 1, "Unknown", 100, -1, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1910)0
   new ParticleOld("delta(1910)0", 21212, 1, "Unknown", 100, 0, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating N(1700)0
   new ParticleOld("N(1700)0", 21214, 1, "Unknown", 100, 0, 1.7, 0.1, -100, 0, -100, -1, -1);

   // Creating N(1535)0
   new ParticleOld("N(1535)0", 22112, 1, "Unknown", 100, 0, 1.535, 0.15, -100, 0, -100, -1, -1);

   // Creating delta(1920)0
   new ParticleOld("delta(1920)0", 22114, 1, "Unknown", 100, 0, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1910)+
   new ParticleOld("delta(1910)+", 22122, 1, "Unknown", 100, 1, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating N(1700)+
   new ParticleOld("N(1700)+", 22124, 1, "Unknown", 100, 1, 1.7, 0.1, -100, 0, -100, -1, -1);

   // Creating N(1535)+
   new ParticleOld("N(1535)+", 22212, 1, "Unknown", 100, 1, 1.535, 0.15, -100, 0, -100, -1, -1);

   // Creating delta(1920)+
   new ParticleOld("delta(1920)+", 22214, 1, "Unknown", 100, 1, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1910)++
   new ParticleOld("delta(1910)++", 22222, 1, "Unknown", 100, 2, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating delta(1920)++
   new ParticleOld("delta(1920)++", 22224, 1, "Unknown", 100, 2, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating sigma(1750)-
   new ParticleOld("sigma(1750)-", 23112, 1, "Unknown", 100, -1, 1.75, 0.09, -100, 0, -100, -1, -1);

   // Creating sigma(1940)-
   new ParticleOld("sigma(1940)-", 23114, 1, "Unknown", 100, -1, 1.94, 0.22, -100, 0, -100, -1, -1);

   // Creating lambda(1600)
   new ParticleOld("lambda(1600)", 23122, 1, "Unknown", 100, 0, 1.6, 0.15, -100, 0, -100, -1, -1);

   // Creating lambda(1890)
   new ParticleOld("lambda(1890)", 23124, 1, "Unknown", 100, 0, 1.89, 0.1, -100, 0, -100, -1, -1);

   // Creating lambda(2110)
   new ParticleOld("lambda(2110)", 23126, 1, "Unknown", 100, 0, 2.11, 0.2, -100, 0, -100, -1, -1);

   // Creating sigma(1750)0
   new ParticleOld("sigma(1750)0", 23212, 1, "Unknown", 100, 0, 1.75, 0.09, -100, 0, -100, -1, -1);

   // Creating sigma(1940)0
   new ParticleOld("sigma(1940)0", 23214, 1, "Unknown", 100, 0, 1.94, 0.22, -100, 0, -100, -1, -1);

   // Creating sigma(1750)+
   new ParticleOld("sigma(1750)+", 23222, 1, "Unknown", 100, 1, 1.75, 0.09, -100, 0, -100, -1, -1);

   // Creating sigma(1940)+
   new ParticleOld("sigma(1940)+", 23224, 1, "Unknown", 100, 1, 1.94, 0.22, -100, 0, -100, -1, -1);

   // Creating xi(1690)-
   new ParticleOld("xi(1690)-", 23314, 1, "Unknown", 100, -1, 1.69, 0.05, -100, 0, -100, -1, -1);

   // Creating xi(1690)0
   new ParticleOld("xi(1690)0", 23324, 1, "Unknown", 100, 0, 1.69, 0.05, -100, 0, -100, -1, -1);

   // Creating rho(1700)0
   new ParticleOld("rho(1700)0", 30113, 1, "Unknown", 100, 0, 1.72, 0.25, -100, 0, -100, -1, -1);

   // Creating rho(1700)+
   new ParticleOld("rho(1700)+", 30213, 1, "Unknown", 100, 1, 1.72, 0.25, -100, 0, -100, -1, -1);

   // Creating omega(1650)
   new ParticleOld("omega(1650)", 30223, 1, "Unknown", 100, 0, 1.67, 0.315, -100, 0, -100, -1, -1);

   // Creating k_star(1680)0
   new ParticleOld("k_star(1680)0", 30313, 1, "Unknown", 100, 0, 1.717, 0.32, -100, 0, -100, -1, -1);

   // Creating k_star(1680)+
   new ParticleOld("k_star(1680)+", 30323, 1, "Unknown", 100, 1, 1.717, 0.32, -100, 0, -100, -1, -1);

   // Creating delta(1600)-
   new ParticleOld("delta(1600)-", 31114, 1, "Unknown", 100, -1, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating N(1720)0
   new ParticleOld("N(1720)0", 31214, 1, "Unknown", 100, 0, 1.72, 0.2, -100, 0, -100, -1, -1);

   // Creating N(1650)0
   new ParticleOld("N(1650)0", 32112, 1, "Unknown", 100, 0, 1.655, 0.165, -100, 0, -100, -1, -1);

   // Creating delta(1600)0
   new ParticleOld("delta(1600)0", 32114, 1, "Unknown", 100, 0, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating N(1720)+
   new ParticleOld("N(1720)+", 32124, 1, "Unknown", 100, 1, 1.72, 0.2, -100, 0, -100, -1, -1);

   // Creating N(1650)+
   new ParticleOld("N(1650)+", 32212, 1, "Unknown", 100, 1, 1.655, 0.165, -100, 0, -100, -1, -1);

   // Creating delta(1600)+
   new ParticleOld("delta(1600)+", 32214, 1, "Unknown", 100, 1, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating delta(1600)++
   new ParticleOld("delta(1600)++", 32224, 1, "Unknown", 100, 2, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating lambda(1670)
   new ParticleOld("lambda(1670)", 33122, 1, "Unknown", 100, 0, 1.67, 0.035, -100, 0, -100, -1, -1);

   // Creating xi(1950)-
   new ParticleOld("xi(1950)-", 33314, 1, "Unknown", 100, -1, 1.95, 0.06, -100, 0, -100, -1, -1);

   // Creating xi(1950)0
   new ParticleOld("xi(1950)0", 33324, 1, "Unknown", 100, 0, 1.95, 0.06, -100, 0, -100, -1, -1);

   // Creating N(1900)0
   new ParticleOld("N(1900)0", 41214, 1, "Unknown", 100, 0, 1.9, 0.5, -100, 0, -100, -1, -1);

   // Creating N(1710)0
   new ParticleOld("N(1710)0", 42112, 1, "Unknown", 100, 0, 1.71, 0.1, -100, 0, -100, -1, -1);

   // Creating N(1900)+
   new ParticleOld("N(1900)+", 42124, 1, "Unknown", 100, 1, 1.9, 0.5, -100, 0, -100, -1, -1);

   // Creating N(1710)+
   new ParticleOld("N(1710)+", 42212, 1, "Unknown", 100, 1, 1.71, 0.1, -100, 0, -100, -1, -1);

   // Creating lambda(1800)
   new ParticleOld("lambda(1800)", 43122, 1, "Unknown", 100, 0, 1.8, 0.3, -100, 0, -100, -1, -1);

   // Creating N(2090)0
   new ParticleOld("N(2090)0", 52114, 1, "Unknown", 100, 0, 2.08, 0.35, -100, 0, -100, -1, -1);

   // Creating N(2090)+
   new ParticleOld("N(2090)+", 52214, 1, "Unknown", 100, 1, 2.08, 0.35, -100, 0, -100, -1, -1);

   // Creating lambda(1810)
   new ParticleOld("lambda(1810)", 53122, 1, "Unknown", 100, 0, 1.81, 0.15, -100, 0, -100, -1, -1);

   // Creating pi(1300)0
   new ParticleOld("pi(1300)0", 100111, 1, "Unknown", 100, 0, 1.3, 0.4, -100, 0, -100, -1, -1);

   // Creating rho(1450)0
   new ParticleOld("rho(1450)0", 100113, 1, "Unknown", 100, 0, 1.465, 0.4, -100, 0, -100, -1, -1);

   // Creating pi(1300)+
   new ParticleOld("pi(1300)+", 100211, 1, "Unknown", 100, 1, 1.3, 0.4, -100, 0, -100, -1, -1);

   // Creating rho(1450)+
   new ParticleOld("rho(1450)+", 100213, 1, "Unknown", 100, 1, 1.465, 0.4, -100, 0, -100, -1, -1);

   // Creating eta(1295)
   new ParticleOld("eta(1295)", 100221, 1, "Unknown", 100, 0, 1.294, 0.055, -100, 0, -100, -1, -1);

   // Creating omega(1420)
   new ParticleOld("omega(1420)", 100223, 1, "Unknown", 100, 0, 1.425, 0.215, -100, 0, -100, -1, -1);

   // Creating k(1460)0
   new ParticleOld("k(1460)0", 100311, 1, "Unknown", 100, 0, 1.46, 0.26, -100, 0, -100, -1, -1);

   // Creating k_star(1410)0
   new ParticleOld("k_star(1410)0", 100313, 1, "Unknown", 100, 0, 1.414, 0.232, -100, 0, -100, -1, -1);

   // Creating k2_star(1980)0
   new ParticleOld("k2_star(1980)0", 100315, 1, "Unknown", 100, 0, 1.973, 0.373, -100, 0, -100, -1, -1);

   // Creating k(1460)+
   new ParticleOld("k(1460)+", 100321, 1, "Unknown", 100, 1, 1.46, 0.26, -100, 0, -100, -1, -1);

   // Creating k_star(1410)+
   new ParticleOld("k_star(1410)+", 100323, 1, "Unknown", 100, 1, 1.414, 0.232, -100, 0, -100, -1, -1);

   // Creating k2_star(1980)+
   new ParticleOld("k2_star(1980)+", 100325, 1, "Unknown", 100, 1, 1.973, 0.373, -100, 0, -100, -1, -1);

   // Creating eta(1475)
   new ParticleOld("eta(1475)", 100331, 1, "Unknown", 100, 0, 1.476, 0.085, -100, 0, -100, -1, -1);

   // Creating phi(1680)
   new ParticleOld("phi(1680)", 100333, 1, "Unknown", 100, 0, 1.68, 0.15, -100, 0, -100, -1, -1);

   // Creating psi'
   new ParticleOld("psi'", 100443, 0, "Unknown", 100, 0, 3.68609, 0, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(100443));
   daughters.clear();
   daughters.push_back(443);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.324,  daughters));
   daughters.clear();
   daughters.push_back(82);
   daughters.push_back(-82);
   part->AddDecay(ParticleOld::Decay(12, 0.1866,  daughters));
   daughters.clear();
   daughters.push_back(443);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.184,  daughters));
   daughters.clear();
   daughters.push_back(10441);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.093,  daughters));
   daughters.clear();
   daughters.push_back(20443);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.087,  daughters));
   daughters.clear();
   daughters.push_back(445);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.078,  daughters));
   daughters.clear();
   daughters.push_back(443);
   daughters.push_back(221);
   part->AddDecay(ParticleOld::Decay(0, 0.027,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-13);
   part->AddDecay(ParticleOld::Decay(0, 0.0083,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(ParticleOld::Decay(0, 0.0083,  daughters));
   daughters.clear();
   daughters.push_back(441);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.0028,  daughters));
   daughters.clear();
   daughters.push_back(443);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.001,  daughters));

   // Creating Upsilon'
   new ParticleOld("Upsilon'", 100553, 0, "Unknown", 100, 0, 10.0233, 0, -100, -1, -100, -1, -1);
   part = const_cast<ParticleOld*>(&ParticleOld::Particles().at(100553));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(21);
   daughters.push_back(21);
   part->AddDecay(ParticleOld::Decay(4, 0.425,  daughters));
   daughters.clear();
   daughters.push_back(553);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(ParticleOld::Decay(0, 0.185,  daughters));
   daughters.clear();
   daughters.push_back(553);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(ParticleOld::Decay(0, 0.088,  daughters));
   daughters.clear();
   daughters.push_back(20553);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.067,  daughters));
   daughters.clear();
   daughters.push_back(555);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.066,  daughters));
   daughters.clear();
   daughters.push_back(10551);
   daughters.push_back(22);
   part->AddDecay(ParticleOld::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-2);
   part->AddDecay(ParticleOld::Decay(32, 0.024,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   part->AddDecay(ParticleOld::Decay(32, 0.024,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(21);
   daughters.push_back(21);
   part->AddDecay(ParticleOld::Decay(4, 0.02,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(ParticleOld::Decay(0, 0.014,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-13);
   part->AddDecay(ParticleOld::Decay(0, 0.014,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-15);
   part->AddDecay(ParticleOld::Decay(0, 0.014,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-1);
   part->AddDecay(ParticleOld::Decay(32, 0.008,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-3);
   part->AddDecay(ParticleOld::Decay(32, 0.008,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
