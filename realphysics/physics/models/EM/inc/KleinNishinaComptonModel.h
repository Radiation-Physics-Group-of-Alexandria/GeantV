
#ifndef KLEINNISHINACOMPTONMODEL_H
#define KLEINNISHINACOMPTONMODEL_H

#include "EMModel.h"

#include <string>

namespace geantphysics {

class Material;
class Element;
class MaterialCuts;
class AliasTable;
class Particle;
class LightTrack;

/**
 * @brief   Compton model for gamma.
 * @class   KleinNishinaComptonModel
 * @author  Marilena Bandieramonte
 * @date    December 2016
 *
 * Compton model for gamma based on Klein-Nishina numerical differential cross sections
 * \cite
 */

class KleinNishinaComptonModel : public EMModel {
public:

  //virtual methods
  virtual void Initialize(); // from EMModel
  //virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle,bool istotal=false);
  //virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  //virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
    
  virtual int    SampleSecondaries(LightTrack &track, std::vector<LightTrack> &sectracks, Geant::GeantTaskData *td);
  virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;



/**
* @name Constructor, destructor:
*/
//@{
   /**
    * @brief Constructor to build a model based on the numerical differential cross sections stored in files.
    *
    * There are 3 different representation of the Seltzer-Berger differential cross sections:
    * - datafileindx = 0 indicates the interpolated NIST data base which is the same as in Geant4.
    * - datafileindx = 1 indicates the NIST data base
    * - datafileindx = 2 indicates that the NRC data base is used: the electron-nuclear part of the differential
    *                    cross section is taken from the NIST data base but the electron-electron part is from
    *                    \cite tessier2008calculation (unlike the NIST free electron-electron interaction this
    *                    later computation includes the effect of binding; differencies expected at lower < 1MeV
    *                    energies).
    *
    * @param[in] iselectron   Flag to indicate that the model is for electron(true) or for psitron(false).
    * @param[in] datafileindx Index to specify the data set.
    * @param[in] modelname    Name of the model.
    */
  KleinNishinaComptonModel(int datafileindx = 0, const std::string &modelname = "KleinNishinaCompton");

  /** @brief Destructor. */
 ~KleinNishinaComptonModel();
//@}


  /**
   * @brief Public method to initilise the model.
   *
   * During the initialisation, Seltzer-Berger numerical differential cross section for bremsstrahlung photon emission
   * will be loaded from files, internal sampling tables for run time sampling of the emitted photon energy will be
   * created. This method always need to be invoked between the construction and usage of the model. -> mb: mantenere?
   *
   */
  void Initialise();

  /**
   * @brief Public method to obtain (restricted) atomic cross sections.
   *
   * Since the model includes some material dependent corrections, for consistency reasons one also needs to provide
   * the material that the element, that the atomic cross section is requested, belongs to.
   *
   * @param[in] elem               Pointer to the element object that the atomic cross section is required.
   * @param[in] mat                Pointer to the material that the element blongs to.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] electronekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Bremsstrahlung atomic cross section for the specified configuration
   *                               in internal [lenght^2] units.
   */
  //double ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gammaprodcutenergy, double electronekin);

  /**
   * @brief Public method to obtain (restricted) macroscopic cross sections.
   *
   * @param[in] mat                Pointer to the material that the cross section is requested.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] electronekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Macroscopic bremsstrahlung cross section for the specified configuration
   *                               in internal [1/lenght] units.
   */
  //double ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin);

  /**
   * @brief Public method to obtain (restricted) radiative stopping power.
   *
   * @param[in] mat                Pointer to the material that the cross section is requested.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] electronekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Radiative stopping power for the specified configuration in internal
   *                               [energy/lenght] units.
   */
  //double ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin);

  /**
   * @brief Public method to sample the emitted (restricted) bremsstrahlung photon energy.
   *
   * The emitted bremsstrahlung photon energy is always higher than the gamma particle kinetic energy production
   * threshold in the specified material-cut pair and not higher than the incident particle kinetic energy.
   *
   * @param[in] matcut     Pointer to the material-cut pair in which the interaction takes place.
   * @param[in] eekin      Kinetic energy of the incident particle i.e. e-/e+.
   * @param[in] r1         Random number distributed uniformly on the 0 1 interval.
   * @param[in] r2         Random number distributed uniformly on the 0 1 interval.
   * @param[in] r3         Random number distributed uniformly on the 0 1 interval.
   * @return               Emitted bremsstrahlung photon energy sampled from the distribution specified by the
   *                       given configuration and the model in internal [energy] units if the particle (e-/e+)
   *                       kinetic energy is higher than the gamma particle production energy threshold. Zero otherwise.
   */
    double SamplePhotonEnergy(const MaterialCuts *matcut, double gammaekin, double r1, double r2, double r3);
    double SamplePhotonEnergyAndDirection(const MaterialCuts *matcut, double gammaekin, double &sinTheta, double &cosTheta, Geant::GeantTaskData *td);

  //
  void   SamplePhotonDirection(double elenergy, double sint2, double onecost, double &sinTheta, double &cosTheta, double rndm);


private:
  /**
   * @brief Internal method to load Seltzer-Berger atomic differential cross sections for bremsstrahlung photon emission
   *        from file.
   *
   *        Used at the initialisation of the model.
   */
  void   LoadDCSData();

  /** @brief Internal method to build emitted photon energy sampling tables under <em>linear approximation of
   *         the p.d.f.</em>.
   *
   *  Used at initialisation of the model to prepare emitted photon energy related transformed variable
   *  sampling tables for all different material-gamma production cut pairs over an e-/e+ kinetic energy grid.
   *  These tables are prepared for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>
   *  within the bins.
   */
  void   InitSamplingTables();

  /** @brief Internal method to build one emitted photon energy sampling tables under <em>linear approximation of
   *         the p.d.f.</em>.
   *
   *  This method is used by SeltzerBergerBremsModel::InitSamplingTables() to build emitted photon energy related
   *  sampling tables for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>. Grid points
   *  are adaptively inserted such that the error between the spline and linearly interpolated p.d.f. is minimised.
   *  Number of grid points will be SeltzerBergerBremsModel::LinAlias::fNumdata which is always the maximum i.e.
   *  SeltzerBergerBremsModel::fNumSamplingPhotEnergies .
   */
  void   BuildOneLinAlias(int indxlalias, const Material *mat, double gcut);

  /**
   * @brief Correction for accounting some differencies between positrons and electrons.
   */
  //double PositronCorrection(double ekinelectron, double ibeta2electron, double ephoton, double z);

  /**
   * @brief Correction for accounting some differencies between positrons and electrons.
   */
  //double PositronCorrection1(double ekinelectron, double ephoton, double gcutener, double z);

private:
 /** @brief Flag to indicate if the model is for e-(true) or for e+(false). Must be set before initialisation. */
 //bool     fIsElectron;                      // flag to indicate if the model is for electron(true) or for positron(flase)
 /** @brief Index of the data file where the numerical DCS are taken from. Must be set before initialisation.
   *        (see more at the description of the constructor).
   */
 int      fDataFileIndx;                    // flag to indicate which SB data should be used

 // secondary related data
 int      fSecondaryInternalCode; // electron GV code set at initialization

 //
 // these are for the DCS data loaded from the file
 /** @brief Maximum atomic number that numerical DCS are available in the file. */
 int      fDCSMaxZet;                       // max Z that we have DCS data in files
 /** @brief Number of electron energies that numerical DCS are available in the file. */
 int      fLoadDCSNumElectronEnergies;      // the elelectron energy grid dimension of the DCS loaded
 /** @brief Number of reduced photon energies at each electron energy that numerical DCS are available in the file. */
 int      fLoadDCSNumReducedPhotonEnergies; // the reduced photon energy grid dimension of the DCS loaded
 /** @brief The electron energy grid that numerical DCS data are given in the file. Loaded from file.
   *        Size of the array is fLoadDCSNumElectronEnergies.
   */
 double  *fLoadDCSElectronEnergyGrid;       // the electron energy grid of the DCS
 /** @brief The reduced photon energy grid that numerical DCS data are given at each electron energy. Loaded from file.
   *         Size of the array is fLoadDCSNumReducedPhotonEnergies.
   */
 double  *fLoadDCSReducedPhotonEnergyGrid;  // the reduced photon energy grid of the DCS
 /** @brief The numerical DCS data loaded from file for each element that used in the detector. Size of the array is
   *         fDCSMaxZet. Each element of the array is a pointer to a double array with size of
   *         fLoadDCSNumElectronEnergies times fLoadDCSNumReducedPhotonEnergies.
   */
double **fLoadDCSForElements;              // container to store DCS data for the used elements

 ///mb: changed
 // these are to describe and define the common gamma kinetic energy grid above we build sampling tables for run-time
 // samling of emitted photon energy.
    
 /** @brief Number of gamma kinetic energy grid points in [KleinNishinaComptonModel::fMinGammaEnergy,
   *        KleinNishinaComptonModel::fMaxGammaEnergy] (default 71). mb: ?? 71
   */
 int     fNumSamplingGammaEnergies;
 /** @brief Number of transformed emitted photon energy related variable in [0,1] (default 54). */
 int     fNumSamplingPhotEnergies;
 /** @brief Minimum of the gamma kinetic energy grid (default 1.0 [keV] that is the minimum available.) */
 double  fMinGammaEnergy;
 /** @brief Maximum of the gamma kinetic energy grid (default 10.0 [GeV] that is the maximum available.) */
 double  fMaxGammaEnergy;
 /** @brief Logarithm of SeltzerBergerBremsModel::fMinGammaEnergy i.e. ln(fMinGammaEnergy) . */
 double  fGammaEnLMin;                         // log min gamma energy
 
  /** @brief Inverse of the gamma kinetic energy grid delta i.e.
   *        ln[fMaxGammaEnergy/fMinGammaEnergy]/(fNumSamplingGammaEnergies-1)
   */
 double  fGammaEnILDelta;                      // 1 log delta gamma energy of the electron energy grid
 
    
  /** @brief The logarithmically spaced gamma kinetic energy grid.
   *
   *        Size of the array is SeltzerBergerBremsModel::fNumSamplingGammaEnergies points in the
   *        [KleinNishinaComptonModel::fMinGammaEnergy, KleinNishinaComptonModel::fMaxGammaEnergy] interval.
   */
 double *fSamplingGammaEnergies;   // the common gamma energy grid which we build sampling tables above
 
  /** @brief The logarithm of KleinNishinaComptonModel::fSamplingGammaEnergies grid.
   *
   *        Size of the array is KleinNishinaComptonModel::fNumSamplingGammaEnergies points in the
   *        [ln(KleinNishinaComptonModel::fMinGammaEnergy), ln(KleinNishinaComptonModel::fMaxGammaEnergy)] interval.
   */
 double *fLSamplingGammaEnergies;            // log of sampling gamma energies

 //
 // data to map all material-electron production cut pair indices to local indices including only the subset of all
 // material-electron production cut that are different. These data used only internally by the model.
 /** @brief Number of all material-electron production cut pairs. */
 int     fNumMaterialCuts;                  // number of different material-electroncut pairs
 /** @brief Number of different material-electron production cut pairs. */
 int     fNumDifferentMaterialECuts;        // number of different material-electroncut pairs
 /** @brief Map from global to local material-electron production cut indices. The size of the array is
   *        KleinNishinaComptonModel::fNumMaterialCuts.
   */
 int    *fGlobalMatGCutIndxToLocal;         // maps the global mat.-cut indices to local indices that are used here

    
 /** @brief Internal data structure to store data for sampling the emitted photon energy.
   *
   *  This data structure is set up at initialisation for each different material-gamma production cut pairs
   *  over the e-/e+ kinetic energy grid to sample the emitted photon energy distribution using a combination
   *  of Walker's alias sampling and liner approximation. At most we will have as many data structure as
   *  SeltzerBergerBremsModel::fNumDifferentMaterialGCuts times SeltzerBergerBremsModel::fNumSamplingElecEnergies
   *  and these data structure pointers are stored in the SeltzerBergerBremsModel::fAliasData linear array.
   *  However, data structures are created only for the possible e-/e+ kinetic energy - emitted photon energy
   *  combinations (i.e. for those e-/e+ kinetic energy grid points in SeltzerBergerBremsModel::fSamplingElecEnergies
   *  that are above the gamma production cut value) and the other elements of SeltzerBergerBremsModel::fAliasData
   *  linear array are left to be nullptr.
   */
 struct LinAlias{
   /** @brief Number of data points i.e. size of the arrays = SeltzerBergerBremsModel::fNumSamplingPhotEnergies. */
   int     fNumdata;                        // size of the arrays
   /** @brief Reduced photon energy related transformed variable values. */
   double *fXdata;                          // reduced photon energies
   /** @brief The probability density function values (not necessarily normalised) over the reduced photon energy
     *        related transformed variable values.
     */
   double *fYdata;                          // p.d.f (not necessarily norm.)
   /** @brief The alias probabilities (not necessarily normalised) over the reduced photon energy related transformed
     *        variable values.
     */
   double *fAliasW;                         // alias probs (not necessarily norm.)
   /** @brief The alias indices over the reduced photon energy related transformed variable values. */
   int    *fAliasIndx; // alias indices
 };
 /** @brief Linear array to store pointers to LinAlias data structures.
   *
   * The size is SeltzerBergerBremsModel::fNumDifferentMaterialGCuts times
   * SeltzerBergerBremsModel::fNumSamplingElecEnergies. Some of the stored pointers are left to be nullptr (that
   * correspond to e-/e+ kinetic energy grid points that are below the gamma production cut value).
   */
 LinAlias   **fAliasData;                   //alias data structure for all different matrial-gammacut pairs
 /** @brief An alias sampler used at run-time sampling of the emitted photon energy. */
 AliasTable  *fAliasSampler;

};

}      // namespace geantphysics

#endif // KLEINNISHINACOMPTON_H
