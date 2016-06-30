
#include "MollerBhabhaIonizationModel.h"

#include <iostream>

#include <cmath>

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"
#include "AliasTable.h"

namespace geant {

MollerBhabhaIonizationModel::MollerBhabhaIonizationModel(bool iselectron, const std::string &modelname)
 : fIsElectron(iselectron), fModelName(modelname) {
   fNumSamplingPrimEnergies =  89;//276; // between the min/max e-/e+ kinetic energies //276->25 /decade
   fNumSamplingElecEnergies =  101; // at each energy grid point we have this number of scattered e- energies
   fMinPrimEnergy           =   1.0*geant::keV; // minimum kinetic energy of the interacting e-/e+
   fMaxPrimEnergy           = 100.0*geant::TeV; // maximum kinetic energy of the interacting e-/e+
   fPrimEnLMin              = 0.0;
   fPrimEnILDelta           = 1.0;
   fSamplingPrimEnergies    = nullptr;
   fLSamplingPrimEnergies   = nullptr;

   fNumMaterialCuts         = 0;
   fNumDifferentElecCuts    = 0;
   fGlobalMatCutIndxToLocal = nullptr;

   fAliasData                = nullptr; //alias data for each different e- cut
   fAliasSampler             = nullptr;

   Initialise();
}

MollerBhabhaIonizationModel::~MollerBhabhaIonizationModel() {
  if (fSamplingPrimEnergies)
    delete [] fSamplingPrimEnergies;
  if (fLSamplingPrimEnergies)
    delete [] fLSamplingPrimEnergies;

  if (fGlobalMatCutIndxToLocal)
    delete [] fGlobalMatCutIndxToLocal;

  if (fAliasData)
    for (int i=0; i<fNumDifferentElecCuts; ++i)
      for (int j=0; j<fNumSamplingPrimEnergies; ++j) {
        int indx = i*fNumSamplingPrimEnergies+j;
        if (fAliasData[indx]) {
          delete [] fAliasData[indx]->fXdata;
          delete [] fAliasData[indx]->fYdata;
          delete [] fAliasData[indx]->fAliasW;
          delete [] fAliasData[indx]->fAliasIndx;
          delete fAliasData[indx];
        }
      }

  if (fAliasSampler)
    delete fAliasSampler;
}

void MollerBhabhaIonizationModel::Initialise() {
  fAliasSampler = new AliasTable();
  InitSamplingTables();
}


/**
  *  The sampling is based on the sampling tables prepared at initialisation. Statistical interpolation is used to
  *  select one of the primary particle kinetic energy grid points out of \f$ E_i \leq E_{kin} < E_{i+1}\f$ (linear
  *  interpolation in log kinetic energy) at a given primary particle kinetic energy \f$E_{kin}\f$. Then the transformed
  *  variable \f$\xi\in[0,1]\f$ is sampled from the sampling table (prepared at initialisation) that belongs to the
  *  selected primary particle kinetic energy grid point. The kinetic energy transfered to the electron \f$T\f$ then is
  *  obtained by applying the following transformation:
  *  \f[
  *     T =
  *     \begin{cases}
  *      T_{cut}^{e-}e^{\xi\ln(0.5E_{kin}/T_{cut}^{e-})} & \textrm{in case of Moller scattering }[e^-+e^-\to e^-+e^-]\\
  *      T_{cut}^{e-}e^{\xi\ln(E_{kin}/T_{cut}^{e-})}    & \textrm{in case of Bhabha scattering }[e^++e^-\to e^++e^-]
  *     \end{cases}
  *  \f]
  *  where \f$E_{kin}\f$ is the current primary particle (i.e. e-/e+) kinetic energy and \f$T_{cut}^{e-}\f$ is the
  *  current electron kinetic energy production threshold.
  */
double MollerBhabhaIonizationModel::SampleEnergyTransfer(MaterialCuts *matcut, double primekin, double r1, double r2, double r3) {
  if (fMaxPrimEnergy<primekin){
    std::cerr<<" **** Primary energy = "<<primekin/geant::GeV<<" [GeV] > fMaxPrimEnergy = "<<fMaxPrimEnergy<<std::endl;
    exit(-1);
  }
  if (primekin<fSamplingPrimEnergies[0])
    primekin = fSamplingPrimEnergies[0];

  double elProdCut = matcut->GetProductionCutsInEnergy()[1]; // e- production cut
  double tmax      = primekin;
  if (fIsElectron)
    tmax *= 0.5;
  if (!(tmax>elProdCut))
    return 0.0;

  int mcindx       = matcut->GetIndex();
  int macindxlocal = fGlobalMatCutIndxToLocal[mcindx];
  // the location of the first-primary-energy lin-alias data of this mut -- e-cut
  int indxstart    = macindxlocal*fNumSamplingPrimEnergies;
  // determine primary energy lower grid point
  double leenergy  = std::log(primekin);
  int eenergyindx  = (int) ((leenergy-fPrimEnLMin)*fPrimEnILDelta);

  if (eenergyindx>=fNumSamplingPrimEnergies-1)
    eenergyindx = fNumSamplingPrimEnergies-2;

  double ploweener = (fLSamplingPrimEnergies[eenergyindx+1]-leenergy)*fPrimEnILDelta;

  indxstart += eenergyindx;
  if (r1>ploweener || !(fAliasData[indxstart]))
     ++indxstart;

  // sample the transformed variable xi=[kappa-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)]
  // where kappa = ln(eps) with eps = T/T_0
  // so xi= [ln(T/T_0)-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)] that is in [0,1]
  double  xi = fAliasSampler->SampleLinear(fAliasData[indxstart]->fXdata, fAliasData[indxstart]->fYdata,
                                     fAliasData[indxstart]->fAliasW, fAliasData[indxstart]->fAliasIndx,
                                     fAliasData[indxstart]->fNumdata,r2,r3);
  // transform it back
  // note: can be further optimised away one log call!
  static const double loghalf = std::log(0.5);
  double dum0 = primekin/elProdCut;
  // we save this log if we store all  log(1/elProdCut) values because log(primekin) is already computed above.
  double dum1 = std::log(dum0);
  if (fIsElectron)
    dum1 += loghalf;
  // return with the sampled kinetic energy transfered to the electron
  return std::exp(xi*dum1)*elProdCut;
}

void MollerBhabhaIonizationModel::InitSamplingTables() {
  // set up the common electron energy grid
  if (fSamplingPrimEnergies) {
    delete [] fSamplingPrimEnergies;
    fSamplingPrimEnergies = nullptr;
  }
  fSamplingPrimEnergies  = new double[fNumSamplingPrimEnergies];
  fLSamplingPrimEnergies = new double[fNumSamplingPrimEnergies];
  double lMin    = std::log(fMinPrimEnergy);
  fPrimEnLMin    = lMin;
  double lMax    = std::log(fMaxPrimEnergy);
  double delta   = (lMax-lMin)/(fNumSamplingPrimEnergies-1.0);
  fPrimEnILDelta = 1.0/delta;
  double dumle   = 0.0;
  for(int i = 0; i<fNumSamplingPrimEnergies; ++i){
    double ddum = lMin+dumle;
    fLSamplingPrimEnergies[i] = ddum;
    fSamplingPrimEnergies[i]  = std::exp(lMin+dumle);
//    std::cerr<<" E("<<i<<") = "<<fSamplingPrimEnergies[i]/geant::MeV<< " [MeV]"<<std::endl;
    dumle+=delta;
  }
  // - get number of different e- production cuts
  // - allocate space and fill the global to local material-cut index map
  const std::vector<MaterialCuts*> theMaterialCutsTable = MaterialCuts::GetTheMaterialCutsTable();
  int fNumMaterialCuts = theMaterialCutsTable.size();
  if (fGlobalMatCutIndxToLocal) {
    delete [] fGlobalMatCutIndxToLocal;
    fGlobalMatCutIndxToLocal = nullptr;
  }
  fGlobalMatCutIndxToLocal = new int[fNumMaterialCuts];
  //std::cerr<<" === Number of global Material-Cuts = "<<fNumMaterialCuts<<std::endl;

  // count diffenet e- production cuts and set the global to local mat-cut index map
  int oldnumDif = fNumDifferentElecCuts;
  int oldnumSPE = fNumSamplingPrimEnergies;
  fNumDifferentElecCuts = 0;
  for (int i=0; i<fNumMaterialCuts; ++i) {
    bool isnew = true;
    int j = 0;
    for (; j<fNumDifferentElecCuts; ++j) {
      if (theMaterialCutsTable[i]->GetProductionCutsInEnergy()[1]==theMaterialCutsTable[j]->GetProductionCutsInEnergy()[1]) {
        isnew = false;
        break;
      }
    }
    if (isnew) {
      fGlobalMatCutIndxToLocal[i] = fNumDifferentElecCuts;
      ++fNumDifferentElecCuts;
    } else {
      fGlobalMatCutIndxToLocal[i] = j;
    }
  }

  // std::cerr<<" === Number of local e- cuts = "<<fNumDifferentElecCuts<<std::endl;
  // allocate space for the different e- cut sampling tables and init these pointers to null
  if (fAliasData)
    for (int i=0; i<oldnumDif; ++i)
      for (int j=0; j<oldnumSPE; ++j) {
        int indx = i*oldnumSPE+j;
        if (fAliasData[indx]) {
          delete [] fAliasData[indx]->fXdata;
          delete [] fAliasData[indx]->fYdata;
          delete [] fAliasData[indx]->fAliasW;
          delete [] fAliasData[indx]->fAliasIndx;
          fAliasData[indx] = nullptr;
        }
      }

  int *isdone = new int[fNumDifferentElecCuts]();
  // as many POSSIBLE alias table as (number of primary enery grig point)x(number of different e- cuts)
  // some of them (those that are kinematically not allowed will remain nullptr!!!)
  int  idum  = fNumDifferentElecCuts*fNumSamplingPrimEnergies;
  fAliasData = new LinAlias*[idum];
  for (int i=0; i<idum; ++i)
    fAliasData[i] = nullptr;

  for (int i=0; i<fNumMaterialCuts; ++i) {
    //std::cerr<<"   See if Material +  e-cut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  e- cut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[1]<<std::endl;
    int localindx = fGlobalMatCutIndxToLocal[i];
    int ialias    = localindx*fNumSamplingPrimEnergies;
    if (!isdone[localindx]) { // init for this e- cut if it has not been done yet
      // std::cerr<<"   -> Will init for Material  - e-cut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  e- cut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[1]<<std::endl;
       BuildOneLinAlias(ialias, theMaterialCutsTable[i]->GetProductionCutsInEnergy()[1]);
       isdone[localindx] = 1;
    }
  }
  delete [] isdone;
  // test
  //for (int i=0; i<fNumMaterialCuts; ++i)
  //  std::cerr<<"     --> Global MatCut-indx = "<< i << " local indx = "<<fGlobalMatGCutIndxToLocal[i] <<std::endl;
}

// build alias tables over the fSamplingPrimEnergies primary particle kinetic energy grid (only at the kinematically
// allowed e- production cut -- primary kinetic energy combinations) for the given e- production cut energy
void MollerBhabhaIonizationModel::BuildOneLinAlias(int ialias, double elecprodcut) {
  // go for the predefined kinetic eneries
  for (int iprimener=0; iprimener<fNumSamplingPrimEnergies; ++iprimener) {
    double primener = fSamplingPrimEnergies[iprimener]; // current primary kinetic energy
    //double tmin     = elecprodcut; // minimum kinetic energy transferable to the scattered e-
    double tmax     = primener;    // minimum kinetic energy transferable to the scattered e-
    if (fIsElectron) {
      tmax *= 0.5;
    }
    // check if it is a kinematically allowed primary-energy & e- production cut combination
    // - build the alias table if yes
    // - leave this alias tabe array point to be nullptr otherwise
    if (tmax>elecprodcut) {
      // create the alias data struct
      fAliasData[ialias] = new LinAlias();
      fAliasData[ialias]->fNumdata   = fNumSamplingElecEnergies;
      fAliasData[ialias]->fXdata     = new double[fNumSamplingElecEnergies]();
      fAliasData[ialias]->fYdata     = new double[fNumSamplingElecEnergies]();
      fAliasData[ialias]->fAliasW    = new double[fNumSamplingElecEnergies]();
      fAliasData[ialias]->fAliasIndx = new    int[fNumSamplingElecEnergies]();
      // fill the x-scale i.e. the transformed xi variable values that are \in [0,1]
      // and the corresponding scattered e- distribution
      double adum = 1.0/(fNumSamplingElecEnergies-1.0);
      for (int i=0; i<fNumSamplingElecEnergies; ++i) {
        double xi = i*adum;
        if (i==0) {
          xi = 0.0;
        } else if (i==fNumSamplingElecEnergies-1) {
          xi = 1.0;
        }
        fAliasData[ialias]->fXdata[i]   = xi;
        if (fIsElectron)
          fAliasData[ialias]->fYdata[i] = ComputeMollerPDF(xi, elecprodcut, primener);
        else
          fAliasData[ialias]->fYdata[i] = ComputeBhabhaPDF(xi, elecprodcut, primener);
      }
      // init the alias data structure for this table
      fAliasSampler->PreparLinearTable(fAliasData[ialias]->fXdata, fAliasData[ialias]->fYdata,
                                       fAliasData[ialias]->fAliasW, fAliasData[ialias]->fAliasIndx,
                                       fNumSamplingElecEnergies);//fAliasData[ialias]->fNumdata);
    }
    ++ialias;
  }
}


/**
 * Ionization part of the (restricted) dE/dx computed based on the formula given by Berger and Seltzer
 * \cite berger1964tables \cite crawford1970electron.
 * \f[
 *   \frac{\mathrm{d}E}{\mathrm{d}x} = \frac{2\pi r_{e}^{2} m_0c^2 n_{el}}{\beta^2}
 *   \left[ \ln \frac{2(\tau+2)}{(Im_0c^2)^2} + G^{\pm}(\tau,\tau_{up}) - \delta
 *   \right]
 * \f]
 * where
 * \f[
 *  \begin{array}{lcl}
 *  r_e       &\to&  \textrm{classical electron radius} \\
 *  m_0c^2    &\to&  \textrm{electron rest mass energy} \\
 *  n_{el}    &\to&  \textrm{electron density of the material} \\
 *  I         &\to&  \textrm{mean excitation energy of the material}(^*) \\
 *  \gamma    &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident particle (e-/e+)
 *                    in rest mass energy unit} \\
 *  \beta     &\to&  \textrm{ratio of the speed of the incident particle to } c;\quad = P_t/E_t \textrm{ with } P_t=pc
 *                   \textrm{ i.e. total momentum of the incident particle (e-/e+) in rest mass energy unit} \\
 *  \tau      &\to&  \textrm{kinetic energy of the incident particle(e-/e+) in rest mass energy unit} \\
 *  \tau_{up} &\to&  \textrm{restricted maximum energy transfer in incident particle rest mass energy unit i.e. }
 *                    \min[\tau_{max},\tau_{c}] \textrm{ where} \\
 *            &&       \tau_{max} \to \textrm{maximum possible energy transfer in incident particle (e-/e+) rest mass energy unit i.e.}\\
 *            &&       \quad \tau_{max} = \begin{cases}
 *                                         \tau           & \quad \textrm{incident particle is e+} \\
 *                                         \frac{\tau}{2} & \quad \textrm{incident particle is e-}
 *                                        \end{cases}\\
 *            &&       \tau_{c} \to \textrm{e- production cut (kinetic) energy in rest mass energy unit}\\
 *  \delta    &\to&  \textrm{density effect correction}(^*) \\
 *  G^{\pm}(\tau,\tau_{up}) &\to& \textrm{represents the term that the integration for the Moller } [e^-+e^-\to e^-+e^-]
 *                                \textrm{ and Bhabha }[e^++e^-\to e^++e^-] \textrm{ differential cross sections yield
 *                                different results i.e.}\\
 *            &&      \textrm{for e-:}\\
 *            &&      \quad G^{-}(\tau,\tau_{up}) =-1-\beta^2+\ln[(\tau-\tau_{up})\tau_{up}]+\frac{\tau}{\tau-\tau_{up}}
 *                                                 +\frac{1}{\gamma^2}
 *                                     \left[
 *                                      \frac{\tau_{up}^2}{2} + (2\tau+1)\ln\left(1-\frac{\tau_{up}}{\tau}\right)
 *                                     \right] \\
 *            &&      \textrm{for e+:}\\
 *            &&      \quad G^{+}(\tau,\tau_{up}) = \ln(\tau\tau_{up})-\frac{\beta^2}{\tau}
 *                          \left[
 *                             \tau+2\tau_{up}-y\frac{3\tau_{up}^2}{2}-y^2\left( \tau_{up}-\frac{\tau_{up}^3}{3} \right)
 *                             -y^3\left( \frac{\tau_{up}^2}{2}-\tau\frac{\tau_{up}^3}{3}+\frac{\tau_{up}^4}{4} \right)
 *                          \right]\\
 *           &&       \quad \textrm{with } y=1/(\gamma+1)
 *  \end{array}
 * \f]
 *
 * (\f$^*\f$see more about \f$I\f$ and \f$\delta\f$ at MaterialProperties)
 */
double MollerBhabhaIonizationModel::ComputeDEDXPerVolume(Material *mat, double prodcutenergy, double particleekin) {
  constexpr double factor     = geant::kTwoPi*geant::kClassicElectronRadius*geant::kClassicElectronRadius*geant::kElectronMassC2;
  const double twolog10inv    = 1.0/(2.0*std::log(10.0));
  // get the material properties
  MaterialProperties *matProp = mat->GetMaterialProperties();
  // get the electron denisty of the material
  double elDensity     = matProp->GetTotalNumOfElectronsPerVol();
  // get the mean excitation energy
  double meanExcEnergy = matProp->GetMeanExcitationEnergy();
  // effective atomic number for the kinetic energy threshold computation
  double effZ          = matProp->GetTotalNumOfElectronsPerVol()/matProp->GetTotalNumOfAtomsPerVol();
  // compute the kinetic energy threshold
  double kineTh        = 0.25*std::sqrt(effZ)*geant::keV;
  // set kinetic energy
  double kineticEnergy = particleekin;
  if (kineticEnergy<kineTh) {
    kineticEnergy = kineTh;
  }
  // set other parameters
  double tau       = kineticEnergy/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
  double gamma     = tau + 1.0;       // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  double gamma2    = gamma*gamma;     // \gamma^2 = [E_t/(mc^2)]^2
  double betagama2 = tau*(tau+2.0);   // = (\beta * \gamma)^2 = (P_t/(mc^2))^2 i.e. [P_t in (mc^2) units]^2
  double beta2     = betagama2/gamma2;// \beta2 i.e. [P_t/E_t]^2
  meanExcEnergy   /= geant::kElectronMassC2;  // mean excitation energy in (mc^2) units
  meanExcEnergy   *= meanExcEnergy;
  // set maximum kinetic energy (in mc^2 units) that can be transformed to a free electron
  double taumax    = tau;  // for e+ : E_kin is the maximum kinetic energy transfer
  if (fIsElectron) {       // for e- : E_kin/2 is the maximum kinetic energy transfer
    taumax *= 0.5;
  }
  // set upper limit of tau: upper limit of the integral corresponding to the continuous part i.e.
  // min(e-/e+ production cut energy, maximum kinetic energy transfer) in mc^2 units
  double tauUpLim  = prodcutenergy/geant::kElectronMassC2;
  if (tauUpLim>taumax) {
    tauUpLim = taumax;
  }
  //
  // compute dE/dx
  double dedx = 0.0;
  // first compute G^{-/+} that is different for e- (Moller) and e+ (Bhabha)
  if (fIsElectron) {
    dedx = std::log((tau-tauUpLim)*tauUpLim) + tau/(tau-tauUpLim)
           + (0.5*tauUpLim*tauUpLim + (2.0*tau + 1.0)*std::log(1.0-tauUpLim/tau))/gamma2 - 1.0 - beta2;
  } else {
    double tauUpLim2 = tauUpLim*tauUpLim*0.5;   // \tau_up^2/2
    double tauUpLim3 = tauUpLim2*tauUpLim/1.5;  // \tau_up^3/3
    double tauUpLim4 = tauUpLim3*tauUpLim*0.75; // \tau_up^4/4
    double y         = 1.0/(1.0 + gamma);
    dedx = std::log(tau*tauUpLim) - beta2*(tau + 2.0*tauUpLim - y*(3.0*tauUpLim2 + y*(tauUpLim - tauUpLim3
           + y*(tauUpLim2 - tau*tauUpLim3 + tauUpLim4))))/tau;
  }
  // add the common term
  dedx        += std::log(2.0*(tau + 2.0)/meanExcEnergy);
  // get the density effect correction term
  double dumx  = std::log(betagama2)*twolog10inv;
  dedx        -= matProp->GetDensityEffectFunctionValue(dumx);
  // apply the multiplicative factor to get the final dE/dx
  dedx        *= factor*elDensity/beta2;
  // final check
  if (dedx<0.0) {
    dedx = 0.0;
  }
  // very low energy extrapolation
  if (particleekin<kineTh) {
    dumx = particleekin/kineTh;
    if (dumx>0.25) {
      dedx /= std::sqrt(dumx);
    } else {
      dedx *= 1.4*std::sqrt(dumx)/(dumx+0.1);
    }
  }
  return dedx;
}


double MollerBhabhaIonizationModel::ComputeXSectionPerAtom(Element *elem, Material *mat, double prodcutenergy,
                                                           double particleekin) {
  double xsec = elem->GetZ()*ComputeXSectionPerElectron(prodcutenergy, particleekin);
  return xsec;
}


double MollerBhabhaIonizationModel::ComputeXSectionPerVolume(Material *mat, double prodcutenergy, double particleekin) {
  double elDensity = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  double xsec      = elDensity*ComputeXSectionPerElectron(prodcutenergy, particleekin);
  return xsec;
}


/**
 *  The (restricted) ionization cross section per electron is computed for Moller \f$[e^-+e^-\to e^-+e^-]\f$ and Bhabha
 *  \f$[e^++e^-\to e^++e^-]\f$ scattering.
 *  Need to multiply by atomic number(Z) of the target atom to get atomic cross section (in internal [lenght^2] units)
 *  and by the electron density of the material to get macroscopic cross section (in units of [1/length]).
 *  Integration of the Moller and Bhabha cross sections between a minimum(\f$T_{kin}^{min}\f$) and
 *  maximum(\f$T_{kin}^{max}\f$) kinetic energy for incident particle kinetic energy
 *  \f$ E_{kin} \in \{T_{kin}^{min},T_{kin}^{max}\}\f$ result in (restricted cross section per \f$e^-\f$):
 *
 *  - for Moller scattering:
 *  \f[
 *   \int_{T_{kin}^{min}=T_{pcut}^{e^-}}^{T_{kin}^{max}=0.5E_{kin}} \frac{\mathrm{d}\sigma^{Moller}(E_{kin})}{\mathrm{d}T}\mathrm{d}T =
 *      \frac{2\pi r_e^2 m_0c^2}{\beta^2 E_{kin}} \left\{
 *          [\varepsilon_2-\varepsilon_1] \left(1-C + \frac{1}{\varepsilon_2\varepsilon_1}
 *               +\frac{1}{(1-\varepsilon_1)(1-\varepsilon_2)}
 *               -C\ln\frac{\varepsilon_2(1-\varepsilon_1)}{\varepsilon_1(1-\varepsilon_2)}
 *          \right)
 *       \right\}
 *  \f]
 *  - for Bhabha scattering:
 *  \f[
 *    \int_{T_{kin}^{min}=T_{pcut}^{e^-}}^{T_{kin}^{max}=E_{kin}} \frac{\mathrm{d}\sigma^{Bhabha}(E_{kin})}{\mathrm{d}T}\mathrm{d}T =
 *      \frac{2\pi r_e^2 m_0c^2}{E_{kin}} \left\{
 *          [\varepsilon_2-\varepsilon_1] \left[
 *            \frac{1}{\beta^2\varepsilon_2\varepsilon_1} + B_2 -\frac{B_3}{2}(\varepsilon_1+\varepsilon_2)
 *            +\frac{B_4}{3}\left[(\varepsilon_1+\varepsilon_2)^2-\varepsilon_1\varepsilon_2 \right]
              \right] - B_1\ln\frac{\varepsilon_2}{\varepsilon_1}
 *       \right\}
 *  \f]
 *  where
 *  \f[
 *  \begin{array}{lcl}
 *  T_{pcut}^{e^-} &\to& \textrm{electron production cut (kinetic) energy} \\
 *  E_{kin}   &\to&  \textrm{kinetic energy of the incident particle i.e. e-(Moller) or e+(Bhabha)} \\
 *  T         &\to&  \textrm{kinetic energy of the scattered e-}\\
 *  r_e       &\to&  \textrm{classical electron radius} \\
 *  m_0c^2    &\to&  \textrm{electron rest mass energy} \\
 *  \gamma    &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident particle (e-/e+)
 *                    in rest mass energy unit} \\
 *  \beta     &\to&  \textrm{ratio of the speed of the incident particle to } c;\quad = P_t/E_t \textrm{ with } P_t=pc
 *                   \textrm{ i.e. total momentum of the incident particle (e-/e+) in rest mass energy unit} \\
 *  C         &\to&  (2\gamma-1)/\gamma^2 \\
 *  \varepsilon_1 &\to& e^- \textrm{production cut (kinetic) energy in incident particle kinetic energy unit} \\
 *  \varepsilon_2 &\to& \textrm{maximum (kinetic)energy transfer in incident particle kinetic energy unit i.e.}\\
 *             &&    \quad = \begin{cases}
 *                     1           & \quad \textrm{Bhabha scattering i.e. incident particle is e+} \\
 *                     \frac{1}{2} & \quad \textrm{Moller scattering i.e. incident particle is e-}
 *                   \end{cases}\\
 *  y         &\to&  1/(1+\gamma)\\
 *  B_1       &\to&  2-y^2\\
 *  B_2       &\to&  (1-2y)(3+y^2)\\
 *  B_4       &\to&  (1-2y)^2\\
 *  B_3       &\to&  B_4+(1-2y)^2
 *  \end{array}
 *  \f]
 *
 * Note, that due to the indistinguishability of the incident and secondary electron in Moller scattering, the electron
 * with the higher energy after the interaction is considered to be the primary so the restricted cross section in
 * Moller scattering is obtained by the integration of the differential Moller scattering cross section between
 * \f$ T_{pcut}^{e^-} \f$ and \f$E_{kin}/2\f$. Therefore, the kinetic energy above a discrete Moller event can accour is
 * \f$ E_{kin}^{TH-Moller} = 2T_{pcut}^{e^-}\f$ (the post interaction primary electron kinetic energy, i.e. the one with
 * higher kinetic energy, will be alway above this threshold kinetic energy). In case of Bhabha scattering, the
 * restricted cross section is obtained by the integration of Bhabha differential cross section between
 * \f$ T_{pcut}^{e^-} \f$ and \f$E_{kin}\f$ i.e. the post interaction positron kinetic energy can be lower than
 * \f$ T_{pcut}^{e^-} \f$.
 */
double MollerBhabhaIonizationModel::ComputeXSectionPerElectron(double prodcutenergy, double particleekin) {
  // secondary e- produced only above production cut energy T_c:
  //   - for Moller scattering: kinetic energy of incident e- should be higher than 2T_c
  //   - for Bhabha scattering: kinetic energy of incident e+ should be higher than T_c
  // otherwise the discrete restricted cross section is zero.
  double maxEkin = particleekin;
  if (fIsElectron) {
    maxEkin *= 0.5;
  }
  double xsec = 0.0;
  if (maxEkin>prodcutenergy) {
    //set min/max energies in incoming particle kinetic energy unit
    double epsmin = prodcutenergy/particleekin;
    double epsmax = maxEkin/particleekin;
    // set other parameters
    double tau       = particleekin/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
    double gamma     = tau + 1.0;            // = E_t/(mc^2)  i.e. E_t in (mc^2) units
    double gamma2    = gamma*gamma;          // \gamma^2 = [E_t/(mc^2)]^2
    double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
    if (fIsElectron) { // Moller scattering i.e. e- + e- -> e- + e-
      double parC = (2.0*gamma-1.0)/gamma2;
      xsec        = (epsmax-epsmin)*(1.0-parC + 1.0/(epsmin*epsmax) + 1.0/((1.0-epsmin)*(1.0-epsmax)))
                    - parC*std::log((epsmax*(1.0-epsmin))/(epsmin*(1.0-epsmax)));
      xsec /= beta2;
    } else {           // Bhabha scattering i.e. e+ + e- -> e+ + e-
      double y     = 1.0/(1.0+gamma);
      double y2    = y*y;
      double ydum  = 1.0-2.0*y;
      double ydum2 = ydum*ydum;
      double b1    = 2.0-y2;
      double b2    = ydum*(3.0+y2);
      double b4    = ydum*ydum2;
      double b3    = b4+ydum2;
      double e1e2  = epsmin*epsmax;
      double e1pe2 = epsmin+epsmax;
      xsec         = (epsmax-epsmin)*(1.0/(beta2*e1e2) + b2 - 0.5*b3*e1pe2 + b4*(e1pe2*e1pe2-e1e2)/3.0)
                     - b1*std::log(epsmax/epsmin);
    }
  }
  xsec *= geant::kTwoPi*geant::kClassicElectronRadius*geant::kClassicElectronRadius*geant::kElectronMassC2/particleekin;
  return xsec;
}


/**
  *  The differential(in fractional energy transfer) atomic cross section for Moller scattering \cite moller1932theorie
  *  \f[
  *     \frac{\mathrm{d}\sigma}{\mathrm{d}\varepsilon} = \frac{2\pi r_e^2 Z}{\beta^2 (\gamma-1)}\left[
  *      C_1 + \frac{1}{\varepsilon}\left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{1}{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right]
  *  \f]
  *  where
  *  \f[
  *  \begin{array}{lcl}
  *  E_{kin}     &\to&  \textrm{kinetic energy of the incident e-} \\
  *  T           &\to&  \textrm{kinetic energy of the scattered e- } \in [T_{cut}^{e-},0.5E_{kin}] \textrm{ where }
  *                     T_{cut} \textrm{ is the electron kinetic energy production threshold}\\
  *  \varepsilon &\to&  \textrm{fractional energy transfer i.e. } \varepsilon=T/E_{kin} \in [T_{cut}^{e-}/E_{kin},0.5]\\
  *  \varepsilon' &\to&  \textrm{fraction of remaining kinetic energy i.e. } \varepsilon'=1-\varepsilon \\
  *  r_e         &\to&  \textrm{classical electron radius} \\
  *  Z           &\to&  \textrm{atomic number of the target atom} \\
  *  \gamma      &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident e- in rest mass
  *                    energy unit} \\
  *  \beta       &\to&  \textrm{ratio of the speed of the incident e- to } c;\quad = P_t/E_t \textrm{ with } P_t=pc \\
  *  C_1         &\to&  \left(\frac{\gamma-1}{\gamma}\right)^2 \\
  *  C_2         &\to&  \frac{2\gamma-1}{\gamma^2}
  *  \end{array}
  *  \f]
  *
  *  The following variable transformations are applied:
  *  - first \f$\varepsilon \to u = \ln(\varepsilon)\f$ so \f$\varepsilon=e^u,\; \mathrm{d}\varepsilon/\mathrm{d}u=e^u=
  *    \varepsilon \f$ which leads to \f$p(u) \propto \varepsilon \left[
  *      C_1 + \frac{1}{\varepsilon}\left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{1}{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right] \f$
  *  - then \f$ u \to \xi = [u-\ln(T_{cut}^{e-}/E_{kin})]/[\ln(0.5E_{kin}/T_{cut}^{e-})] \in [0,1]\f$ so
  *    \f$ u = \xi\ln(0.5E_{kin}/T_{cut}^{e-})+\ln(T_{cut}^{e-}/E_{kin}),\;
  *    \mathrm{d}u/\mathrm{d}\xi = \ln(0.5E_{kin}/T_{cut}^{e-})\f$ which leads to \f$ p(\xi) \propto \varepsilon \left[
  *      C_1 + \frac{1}{\varepsilon}\left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{1}{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right] \ln(0.5E_{kin}/T_{cut}^{e-}) \f$ where the last factor is just a constant.
  *
  *  So the transformed distribution
  *  \f[
  *      p(\xi) \propto \left[
  *      \varepsilon C_1 + \left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{\varepsilon }{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right]
  *  \f]
  *  where the transformed variable in terms of \f$T\f$ kinetic energy transfer is
  *  \f$\xi = \ln(T/T_{cut}^{e-})/\ln(0.5E_{kin}/T_{cut}^{e-}) \f$ so after the sampling of
  *  \f$\xi\f$ the kinetic energy transfer can be obtained as \f$T=T_{cut}^{e-}e^{\xi\ln(0.5E_{kin}/T_{cut}^{e-})}\f$.
  *
  */
double MollerBhabhaIonizationModel::ComputeMollerPDF(double xi, double prodcutenergy, double particleekin) {
  double tau       = particleekin/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
  double gamma     = tau + 1.0;            // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  double gamma2    = gamma*gamma;          // \gamma^2 = [E_t/(mc^2)]^2
  double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
  double C1        = (gamma-1.0)/gamma;
  C1 *=C1;
  double C2        = (2.0*gamma-1.0)/gamma2;

  double dum0      = prodcutenergy/particleekin;
  double dum1      = std::log(0.5/dum0);
  double a         = std::exp(xi*dum1)*dum0; // this is eps =  exp(xi*ln(0.5*T_0/T_cut))*T_cut/T_0
  double b         = 1.0-a;                  // eps'

  return ((1.0/a-C2)+a*C1+a/b*(1.0/b-C2)) *dum0; // xdum0 is just scaling; this is the shape
}


/**
  *  The differential(in fractional energy transfer) atomic cross section for Bhabha scattering
  *  \cite bhabha1936scattering
  *  \f[
  *     \frac{\mathrm{d}\sigma}{\mathrm{d}\varepsilon} = \frac{2\pi r_e^2 Z}{(\gamma-1)}\left[
  *      \frac{1}{\beta^2\varepsilon^2} - \frac{B_1}{\varepsilon} +  B_2 - \varepsilon B_3 + \varepsilon^2 B_4
  *     \right]
  *  \f]
  *  where
  *  \f[
  *  \begin{array}{lcl}
  *  E_{kin}     &\to&  \textrm{kinetic energy of the incident e+} \\
  *  T           &\to&  \textrm{kinetic energy of the scattered e- } \in [T_{cut}^{e-},E_{kin}] \textrm{ where }
  *                     T_{cut} \textrm{ is the electron kinetic energy production threshold}\\
  *  \varepsilon &\to&  \textrm{fractional energy transfer i.e. } \varepsilon=T/E_{kin} \in [T_{cut}^{e-}/E_{kin},1]\\
  *  r_e         &\to&  \textrm{classical electron radius} \\
  *  Z           &\to&  \textrm{atomic number of the target atom} \\
  *  \gamma      &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident e- in rest mass
  *                    energy unit} \\
  *  \beta       &\to&  \textrm{ratio of the speed of the incident e- to } c;\quad = P_t/E_t \textrm{ with } P_t=pc \\
  *  y           &\to&  1/(1+\gamma)\\
  *  B_1         &\to&  2-y^2\\
  *  B_2         &\to&  (1-2y)(3+y^2)\\
  *  B_4         &\to&  (1-2y)^2\\
  *  B_3         &\to&  B_4+(1-2y)^2
  *  \end{array}
  *  \f]
  *
  *  The following variable transformations are applied:
  *  - first \f$\varepsilon \to u = \ln(\varepsilon)\f$ so \f$\varepsilon=e^u,\; \mathrm{d}\varepsilon/\mathrm{d}u=e^u=
  *    \varepsilon \f$ which leads to \f$p(u) \propto \varepsilon \left[
  *      \frac{1}{\beta^2\varepsilon^2} - \frac{B_1}{\varepsilon} +  B_2 - \varepsilon B_3 + \varepsilon^2 B_4
  *     \right] \f$
  *  - then \f$ u \to \xi = [u-\ln(T_{cut}^{e-}/E_{kin})]/[\ln(E_{kin}/T_{cut}^{e-})] \in [0,1]\f$ so
  *    \f$ u = \xi\ln(E_{kin}/T_{cut}^{e-})+\ln(T_{cut}^{e-}/E_{kin}),\;
  *    \mathrm{d}u/\mathrm{d}\xi = \ln(E_{kin}/T_{cut}^{e-})\f$ which leads to \f$ p(\xi) \propto \varepsilon \left[
  *      \frac{1}{\beta^2\varepsilon^2} - \frac{B_1}{\varepsilon} +  B_2 - \varepsilon B_3 + \varepsilon^2 B_4
  *     \right] \ln(E_{kin}/T_{cut}^{e-}) \f$ where the last factor is just a constant.
  *
  *  So the transformed distribution
  *  \f[
  *      p(\xi) \propto \left[
  *      \frac{1}{\beta^2\varepsilon} - B_1 +  \varepsilon ( B_2 - \varepsilon ( B_3 + \varepsilon B_4))
  *     \right]
  *  \f]
  *  where the transformed variable in terms of \f$T\f$ kinetic energy transfer is
  *  \f$\xi = \ln(T/T_{cut}^{e-})/\ln(E_{kin}/T_{cut}^{e-}) \f$ so after the sampling of
  *  \f$\xi\f$ the kinetic energy transfer can be obtained as \f$T=T_{cut}^{e-}e^{\xi\ln(E_{kin}/T_{cut}^{e-})}\f$.
  *
  */
double MollerBhabhaIonizationModel::ComputeBhabhaPDF(double xi, double prodcutenergy, double particleekin) {
  double tau       = particleekin/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
  double gamma     = tau + 1.0;            // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  double gamma2    = gamma*gamma;          // \gamma^2 = [E_t/(mc^2)]^2
  double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
  double y         = 1.0/(1.0+gamma);
  double y2        = y*y;
  double ydum      = 1.0-2.0*y;
  double ydum2     = ydum*ydum;
  double b1        = 2.0-y2;
  double b2        = ydum*(3.0+y2);
  double b4        = ydum*ydum2;
  double b3        = b4+ydum2;

  double dum0      = prodcutenergy/particleekin;
  double dum1      = std::log(1.0/dum0);
  double a         = std::exp(xi*dum1)*dum0; // this is eps =  = exp(xi*ln(T_0/T_cut))*T_cut/T_0

  return ((1.0/(a*beta2)-b1) + a*(b2+a*(a*b4-b3)));//
}

} // namespace geant