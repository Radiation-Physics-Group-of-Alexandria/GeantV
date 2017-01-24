
#ifndef SAUTERGAVRILAPHOTOELECTRICMODEL_H
#define SAUTERGAVRILAPHOTOELECTRICMODEL_H

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
     * @brief   Photoelectric model for gamma.
     * @class   SauterGavrilaPhotoElectricModel
     * @author  Marilena Bandieramonte
     * @date    December 2016
     *
     * PhotoElectric model for gamma based on SauterGavrila differential cross section
     * \cite
     */
    
    class SauterGavrilaPhotoElectricModel : public EMModel {
    public:
        
        //virtual methods
        virtual void Initialize(); // from EMModel
        
        virtual int SampleSecondaries(LightTrack &track, std::vector<LightTrack> &sectracks, Geant::GeantTaskData *td);
        
        virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;
        
        //Constructor
        SauterGavrilaPhotoElectricModel(int datafileindx = 0, const std::string &modelname = "SauterGavrilaPhotoElectric");
        //Deconstructor
        ~SauterGavrilaPhotoElectricModel();
    
        //Initialization
        void Initialise();
        
        //Differential cross section based on SauterGavrila distribution for k-shell
        double CalculateDiffCrossSection(int /*Zelement*/, double energy0, double cosTheta);
        
        //Cross-section per atom
        double ComputeCrossSectionPerAtom(LightTrack &track,
                                          double Z,
                                          double A,
                                          double,
                                          double);
        //Cross-section per volume
        virtual double CrossSectionPerVolume(const MaterialCuts *matcut,
                                             LightTrack &track,
                                             double kineticEnergy,
                                             double cutEnergy,
                                             double maxEnergy);
    
        //Sample PhotoElectron direction with Alias sampling
        double SamplePhotoElectronDirection_Alias(const MaterialCuts *matcut,
                                                  double gammaEnIn,
                                                  double r1,
                                                  double r2,
                                                  double r3);

        //Sample PhotoElectron direction with Alias sampling
        void SamplePhotoElectronDirection_Rejection(double gammaEkinIn,
                                                double &sinTheta,
                                                double &cosTheta,
                                                double &phi,
                                                Geant::GeantTaskData *td);
        
    private:
        
        void   InitSamplingTables();
        void   BuildOneLinAlias(int indxlalias, double gcut);
        
        int      fSecondaryInternalCode; // electron GV code set at initialization
        
        /** @brief Number of gamma kinetic energy grid points in [SauterGavrilaPhotoElectricModel::fMinGammaEnergy,
         *        SauterGavrilaPhotoElectricModel::fMaxGammaEnergy] (default 71).
         */
        int     fNumSamplingGammaEnergies;
        
        /** @brief Number of transformed emitted photon energy related variable in [0,1] (default 54). */
        int     fNumSamplingElectronEnergies;
        
        /** @brief Minimum of the gamma kinetic energy grid (default 1.0 [keV] that is the minimum available.) */
        double  fMinGammaEnergy;
        
        /** @brief Maximum of the gamma kinetic energy grid (default 10.0 [GeV] that is the maximum available.) */
        double  fMaxGammaEnergy;
        
        /** @brief Logarithm of SauterGavrilaPhotoElectricModel::fMinGammaEnergy i.e. ln(fMinGammaEnergy) . */
        double  fGammaEnLMin;                         // log min gamma energy
        
        /** @brief Inverse of the gamma kinetic energy grid delta i.e.
         *        ln[fMaxGammaEnergy/fMinGammaEnergy]/(fNumSamplingGammaEnergies-1)
         */
        double  fGammaEnILDelta;                      // 1 log delta gamma energy of the electron energy grid
        
        
        /** @brief The logarithmically spaced gamma kinetic energy grid.
         *
         *        Size of the array is SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies points in the
         *        [SauterGavrilaPhotoElectricModel::fMinGammaEnergy, SauterGavrilaPhotoElectricModel::fMaxGammaEnergy] interval.
         */
        double *fSamplingGammaEnergies;   // the common gamma energy grid which we build sampling tables above
        
        
        /** @brief The logarithm of SauterGavrilaPhotoElectricModel::fSamplingGammaEnergies grid.
         *
         *  Size of the array is SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies points in the
         *  [ln(SauterGavrilaPhotoElectricModel::fMinGammaEnergy), ln(SauterGavrilaPhotoElectricModel::fMaxGammaEnergy)]
         *  interval.
         */
        double *fLSamplingGammaEnergies;            // log of sampling gamma energies
        
        // data to map all material-electron production cut pair indices to local indices including only the subset of all
        // material-electron production cut that are different. These data used only internally by the model.
        /** @brief Number of all material-electron production cut pairs. */
        int     fNumMaterialCuts;                  // number of different material-electroncut pairs
        
        /** @brief Number of different material-electron production cut pairs. */
        int     fNumDifferentMaterialECuts;        // number of different material-electroncut pairs
        
        /** @brief Map from global to local material-electron production cut indices. The size of the array is
         *        SauterGavrilaPhotoElectricModel::fNumMaterialCuts.
         */
        int    *fGlobalMatGCutIndxToLocal;         // maps the global mat.-cut indices to local indices that are used here
        
        
        /** @brief Internal data structure to store data for sampling the emitted photoelectron direction
         *
         *  This data structure is set up at initialisation for each different material-gamma production cut pairs
         *  over the e-/e+ kinetic energy grid to sample the emitted photon energy distribution using a combination
         *  of Walker's alias sampling and liner approximation. At most we SauterGavrilaPhotoElectricModel have as many data structure as
         *  SauterGavrilaPhotoElectricModel::fNumDifferentMaterialGCuts times SauterGavrilaPhotoElectricModel::fNumSamplingElecEnergies
         *  and these data structure pointers are stored in the SauterGavrilaPhotoElectricModel::fAliasData linear array.
         *  However, data structures are created only for the possible e-/e+ kinetic energy - emitted photon energy
         *  combinations (i.e. for those e-/e+ kinetic energy grid points in SauterGavrilaPhotoElectricModel::fSamplingElecEnergies
         *  that are above the gamma production cut value) and the other elements of SauterGavrilaPhotoElectricModel::fAliasData
         *  linear array are left to be nullptr.
         */
        struct LinAlias{
            /** @brief Number of data points i.e. size of the arrays = SauterGavrilaPhotoElectricModel::fNumSamplingPhotEnergies. */
            int     fNumdata;                        // size of the arrays
            
            /** @brief This must be the gamma energies - maybe beginning of the bin? (transformed or not..it depends). */
            double *fXdata;                          // reduced photon energies
            /** @brief The probability density function values (not necessarily normalised) over the photon energy
             *        variable values.
             */
            double *fYdata;                          // p.d.f (not necessarily norm.)
            /** @brief The alias probabilities (not necessarily normalised) over the reduced photon energy
             *        variable values.
             */
            double *fAliasW;                         // alias probs (not necessarily norm.)
            /** @brief The alias indices over the photon energy variable values. */
            int    *fAliasIndx; // alias indices
        };
        
        /** @brief Linear array to store pointers to LinAlias data structures.
         *
         * The size is SauterGavrilaPhotoElectricModel::fNumDifferentMaterialGCuts times
         * SauterGavrilaPhotoElectricModel::fNumSamplingElecEnergies. Some of the stored pointers are left to be nullptr (that
         * correspond to e-/e+ kinetic energy grid points that are below the gamma production cut value).
         */
        LinAlias   **fAliasData;                   //alias data structure for all different material-cut pairs
        /** @brief An alias sampler used at run-time sampling of the emitted photon energy. */
        AliasTable  *fAliasSampler;
        
    };
    
}      // namespace geantphysics

#endif // SAUTERGAVRILAPHOTOELECTRICMODEL_H
