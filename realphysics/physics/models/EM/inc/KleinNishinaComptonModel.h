
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
        virtual int    SampleSecondaries(LightTrack &track, std::vector<LightTrack> &sectracks, Geant::GeantTaskData *td);
        virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;
        
        
        /**
         * @brief Constructor to build a model based on the numerical differential cross sections calculated with the KleinNishina.
         *
         * @param[in] modelname    Name of the model.
         */
        KleinNishinaComptonModel(const std::string &modelname = "KleinNishinaCompton");
        
        /** @brief Destructor. */
        
        ~KleinNishinaComptonModel();
        
        
        
        /**
         * @brief Public method to initilise the model.
         *
         * During the initialisation, internal sampling tables for run time sampling of the outgoing photon energy will be
         * created. This method always need to be invoked between the construction and usage of the model.
         *
         */
        void Initialise();
        
        
        /**
         * @brief Public method to calculate ...
         *
         * @param[in]
         * @param[in]
         * @param[in]
         * @return
         *
         */
        double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle * );
        
        
        /**
         * @brief Public method to calculate ...
         *
         * @param[in]
         * @param[in]
         * @param[in]
         * @return
         *
         */
        double ComputeXSectionPerVolume(const Material *mat, double prodcutenergy, double particleekin);
        
        /**
         * @brief Public method to calculate ...
         *
         * @param[in]
         * @param[in]
         * @param[in]
         * @return
         *
         */
        double ComputeXSectionPerAtom(const Element *elem, double kinenergy);
        
        /**
         * @brief Public method to calculate the KleinNishinaCompton differential cross section
         *
         * @param[in] energy0       Incomming photon energy.
         * @param[in] energy1       Scattered photon energy.
         * @return    dsigma        Differential cross section.
         *
         */
        double CalculateDiffCrossSection(double energy0, double energy1);
        
        /**
         * @brief Public method to sample the outgoing compton photon energy using the Alias Table
         *
         *
         * @param[in] matcut     Pointer to the material-cut pair in which the interaction takes place.
         * @param[in] gammaekin  Kinetic energy of the incident gamma particle.
         * @param[in] r1         Random number distributed uniformly in the 0 1 interval.
         * @param[in] r2         Random number distributed uniformly in the 0 1 interval.
         * @param[in] r3         Random number distributed uniformly in the 0 1 interval.
         * @return               outgoing compton photon energy sampled from the distribution specified by the
         *                       given configuration and the model in internal [energy] units if the gamma
         *                       kinetic energy is higher than the gamma particle production energy threshold.
         Zero otherwise.
         */
        
        double SamplePhotonEnergy(const MaterialCuts *matcut, double gammaekin, double r1, double r2, double r3);
        
        
        /**
         * @brief Public method to sample the outgoing compton photon energy and direction of the outgoing particle using the composition-rejection sampling
         *
         * @param[in] gammaekin  Kinetic energy of the incident particle gamma.
         * @param[in] sinTheta   Sin of the scattering angle (theta) of the gamma in the interaction/process frame coordinate system.
         * @param[in] cosTheta   Cos of the scattering angle (theta) of the gamma in the interaction/process frame coordinate system.
         * @param[in] td         GeantTaskData necessary for random numbers generation
         * @return               outgoing compton photon energy sampled from the distribution specified by the
         *                       given configuration and the model in internal [energy] units if the gamma
         *                       kinetic energy is higher than the gamma particle production energy threshold. Zero otherwise.
         */
        double SamplePhotonEnergyAndDirection(double gammaekin, double &sinTheta, double &cosTheta, Geant::GeantTaskData *td);
        
        
        /**
         * @brief Public method to sample the outgoing compton photon energy and direction of the outgoing particle using the composition-rejection sampling
         *
         * @param[in] gammaEkinIn   Kinetic energy of the incident gamma particle.
         * @param[in] gammaEkinOut  Kinetic energy of the outgoing gamma particle.
         * @param[in] sinTheta      Sin of the scattering angle (theta) of the gamma in the interaction/process frame coordinate system.
         * @param[in] cosTheta      Cos of the scattering angle (theta) of the gamma in the interaction/process frame coordinate system.
         * @param[in] phi           Phi angle of the scattered gamma in the interaction/process frame coordinate system.
         * @param[in] rnd           Random number distributed uniformly in the 0 1 interval.
         */
        void CalculatePhotonDirection(double gammaEkinIn, double gammaEkinOut, double &sinTheta, double &cosTheta, double &phi, double rnd);
        
    private:
        
        /** @brief Internal method to build emitted photon energy sampling tables under <em>linear approximation of
         *         the p.d.f.</em>.
         *
         *  Used during the initialisation of the model to prepare emitted photon energy
         *  sampling tables for all different material-gamma production cut pairs over a gamma kinetic energy grid.
         *  These tables are prepared for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>
         *  within the bins.
         */
        void   InitSamplingTables();
        
        /** @brief Internal method to build one emitted photon energy sampling tables under <em>linear approximation of
         *         the p.d.f.</em>.
         *
         *  This method is used by KleinNishinaComptonModel::InitSamplingTables() to build emitted photon energy related
         *  sampling tables for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>. Grid points
         *  are adaptively inserted such that the error between the differential cross section and linearly interpolated p.d.f. is minimised.
         *  Number of grid points will be KleinNishinaComptonModel::LinAlias::fNumdata which is always the maximum i.e.
         *  KleinNishinaComptonModel::fNumSamplingElectronEnergies.
         */
        void   BuildOneLinAlias(int indxlalias, double gcut);
        
        
        // secondary related data
        int      fSecondaryInternalCode; // electron GV code set at initialization
        
        
        // The following members are used to describe and define the common gamma kinetic energy grid above we build
        // sampling tables for run-time sampling of emitted photon energy.
        
        /** @brief Number of gamma kinetic energy grid points in [KleinNishinaComptonModel::fMinGammaEnergy,
         *        KleinNishinaComptonModel::fMaxGammaEnergy] (default 71).
         */
        int     fNumSamplingGammaEnergies;
        
        /** @brief Number of transformed emitted photon energy related variable in [0,1] (default 54). */
        int     fNumSamplingElectronEnergies;
        
        /** @brief Minimum of the gamma kinetic energy grid (default 1.0 [keV] that is the minimum available.) */
        double  fMinGammaEnergy;
        
        /** @brief Maximum of the gamma kinetic energy grid (default 10.0 [GeV] that is the maximum available.) */
        double  fMaxGammaEnergy;
        
        /** @brief Logarithm of KleinNishinaComptonModel::fMinGammaEnergy i.e. ln(fMinGammaEnergy) . */
        double  fGammaEnLMin;
        
        /** @brief Inverse of the gamma kinetic energy grid delta i.e.
         *        ln[fMaxGammaEnergy/fMinGammaEnergy]/(fNumSamplingGammaEnergies-1)
         */
        double  fGammaEnILDelta;
        
        /** @brief The logarithmically spaced gamma kinetic energy grid.
         *
         *        Size of the array is KleinNishinaComptonModel::fNumSamplingGammaEnergies points in the
         *        [KleinNishinaComptonModel::fMinGammaEnergy, KleinNishinaComptonModel::fMaxGammaEnergy] interval.
         */
        double *fSamplingGammaEnergies;   // the common gamma energy grid which we build sampling tables above
        
        /** @brief The logarithm of KleinNishinaComptonModel::fSamplingGammaEnergies grid.
         *
         *        Size of the array is KleinNishinaComptonModel::fNumSamplingGammaEnergies points in the
         *        [ln(KleinNishinaComptonModel::fMinGammaEnergy), ln(KleinNishinaComptonModel::fMaxGammaEnergy)] interval.
         */
        double *fLSamplingGammaEnergies;            // log of sampling gamma energies
        
        
        // The following are data to map all material-electron production cut pair indices to local indices including only
        // the subset of all material-electron production cut that are different. These data used only internally by the
        // model.
        
        /** @brief Number of all material-electron production cut pairs. */
        int     fNumMaterialCuts;                  // number of different material-electroncut pairs
        
        /** @brief Number of different material-electron production cut pairs. */
        int     fNumDifferentMaterialECuts;        // number of different material-electroncut pairs
        
        /** @brief Map from global to local material-electron production cut indices. The size of the array is
         *        KleinNishinaComptonModel::fNumMaterialCuts.
         */
        int    *fGlobalMatGCutIndxToLocal;         // maps the global mat.-cut indices to local indices that are used here
        
        
        /** @brief Internal data structure to store data for sampling the emitted photon energy. --> emitted electron energy
         *
         *  This data structure is set up at initialisation for each different material-gamma production cut pairs
         *  over the e-/e+ kinetic energy grid to sample the emitted photon energy distribution using a combination
         *  of Walker's alias sampling and liner approximation. At most we KleinNishinaComptonModel have as many data structure as
         *  KleinNishinaComptonModel::fNumDifferentMaterialGCuts times KleinNishinaComptonModel::fNumSamplingElecEnergies
         *  and these data structure pointers are stored in the KleinNishinaComptonModel::fAliasData linear array.
         *  However, data structures are created only for the possible e-/e+ kinetic energy - emitted photon energy
         *  combinations (i.e. for those e-/e+ kinetic energy grid points in KleinNishinaComptonModel::fSamplingElecEnergies
         *  that are above the gamma production cut value) and the other elements of KleinNishinaComptonModel::fAliasData
         *  linear array are left to be nullptr.
         */
        struct LinAlias{
            /** @brief Number of data points i.e. size of the arrays = KleinNishinaComptonModel::fNumSamplingPhotEnergies. */
            int     fNumdata;                        // size of the arrays
            
            /** @brief This must be the gamma energies - maybe beginning of the bin? (transformed or not..it depends). */
            double *fXdata;                          // reduced photon energies
            
            /** @brief The probability density function values (not necessarily normalised) over the photon energy
             *        variable values.
             */
            double *fYdata;                          // p.d.f (not necessarily normalised)
            
            /** @brief The alias probabilities (not necessarily normalised) over the reduced photon energy
             *        variable values.
             */
            double *fAliasW;                         // alias probs (not necessarily normalised)
            
            /** @brief The alias indices over the photon energy variable values. */
            int    *fAliasIndx; // alias indices
        };
        
        /** @brief Linear array to store pointers to LinAlias data structures.
         *
         * The size is KleinNishinaComptonModel::fNumDifferentMaterialGCuts times
         * KleinNishinaComptonModel::fNumSamplingElecEnergies. Some of the stored pointers are left to be nullptr (that
         * correspond to e-/e+ kinetic energy grid points that are below the gamma production cut value).
         */
        LinAlias   **fAliasData;                   //alias data structure for all different material-gammacut pairs
        
        /** @brief An alias sampler used at run-time sampling of the emitted photon energy. */
        AliasTable  *fAliasSampler;
        
    };
    
}      // namespace geantphysics

#endif // KLEINNISHINACOMPTON_H
