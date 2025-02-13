#ifndef HADRONSAMPLING_H
#define HADRONSAMPLING_H

#include "../../../Utilities/Include/ChargeConfiguration.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../HadronizationHadronSet/Include/HadronData.h"

#include <vector>
#include <map>

using namespace std;

// TODO: Introduce all possible simplification in sampling weight calculation (e.g., 2J+1)
// TODO: Check usage or return in inlined function calling deeper function (remind error in CSEM getDomain)
// TODO: add equation for weight calculation and sampling strategy algorithm description
// TODO: describe adopted assumption: single heavy + light, etc, with link to methods below
// TODO: in all doxy long comment use . and ,!

/**
* @brief Hadron sampling class
*
* This class is used for hadron sampling during cluster hadronization
* channel generation. The sampling is performed by the run method,
* which generates a hadron set compatible with the provided input cluster
* charge configuration (and mass value). 
*
* @author Christopher Bignamini
*/
class HadronSampling
{
    
    public:
    
        /**
        * @brief Constructor
        *
        * @param i_hadronSamplingGroups Hadron sampling group handler used for channel composition sampling
        * @param i_samplingTemperature Hadron sampling temperature parameter
        * @param i_samplingEnergyDensity Hadron sampling cluster energy density 
        * @param io_randomGenerator Random number generator
        * @throw HadronizationException in case of non positive sampling temperature or energy density 
        *        or errors during preliminary sampling data calculation
        */
        HadronSampling(const HadronSamplingGroups& i_hadronSamplingGroups,
                       double i_samplingTemperature,
                       double i_samplingEnergyDensity,
                       RandomNumberGenerator& io_randomGenerator);
    
        /**
        * @brief Destructor
        *
        */
        ~HadronSampling(void);
    
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Hadron set sampling method
        *
        * This is the hadron set sampling method. After a (successful) execution of this method 
        * a hadronization hadron set for the provided input cluster is available as well as
        * a sampling weight.
        *
        * WARNING: at present only light flavored and single heavy flavored cluster are considered 
        * (e.g., ccbar and cbbar cluster are not accepted).
        *
        * @param i_chargeConfiguration Cluster charge configuration
        * @param i_mass Cluster mass
        * @param i_energyDensity Cluster energy density
        * @return 0 in case of sampling correctly performed
        *         1 in case of maximum number of sampling attempt exceeded
        *         2 in case of heavy flavored hadron unavailability for the hadronization 
        *           of a given heavy flavored cluster (due to mass values)
        *         3 in case of double heavy flavored cluster provided
        * @throw HadronizationException if class has not been correctly initialized or in case of errors during the sampling procedure
        */
        unsigned int run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                         double i_mass,
                         double i_energyDensity);

        /**
        * @brief Hadron sampling temperature parameter set method
        *
        * @param i_hadronSamplingTemperature Hadron sampling temperature parameter
        * @throw HadronizationException if non positive sampling temperature is provided or 
        *        in case of errors during preliminary sampling data calculation
        */
        void setHadronSamplingTemperature(double i_hadronSamplingTemperature);
    
        /**
        * @brief Hadron sampling cluster energy density parameter set method
        *
        * @param i_hadronSamplingEnergyDensity Hadron sampling cluster energy density parameter
        * @throw HadronizationException if non positive sampling energy density is provided or 
        *        in case of errors during preliminary sampling data calculation
        */
        void setHadronSamplingEnergyDensity(double i_hadronSamplingEnergyDensity);

        /**
        * @brief Hadron set get method
        *
        * @return Generated hadron set
        * @throw HadronizationException if hadron sampling has not been performed
        */
    // TODO: debug, after reordering removal this method must be const
        const vector<HadronData>& getHadronSet(void);
    
        /**
        * @brief Number of hadron set get method
        *
        * @return Generated hadron set multiplicity
        * @throw HadronizationException if hadron sampling has not been performed
        */
        unsigned int getNumberOfHadrons(void) const;
    
        /**
        * @brief Hadron sampling weight get method
        *
        * @return Hadron set sampling weight
        * @throw HadronizationException if hadron sampling has not been performed
        */
        double getSamplingWeight(void) const;
    
        // TEST
        inline void updateRapidities(void) { m_runNewCalculation=true; }
        // TEST
    
    private:
    
        /**
        * @brief Hadron group sampling probability calculation method
        *
        * @throw HadronizationException in case of error during hadron sampling data calculation 
        *        (e.g., zero mean multiplicity) or not identified hadron sampling group (stored 
        *        in m_hadronSamplingGroups).
        */
        void computeSamplingData(void);

        /**
        * @brief Heavy hadron sampling probability calculation method
        *
        * This method is internally used during each heavy cluster hadronization procedure
        * for the calculation of the heavy hadron sampling probability for a given cluster 
        * heavy charge configuration: only flavor compatible hadrons with mass smaller than
        * the cluster mass are considered. 
        *
        * @param i_clusterCharmCharge Cluster charm charge
        * @param i_clusterBottomCharge Cluster bottom charge
        * @return true if if no heavy hadron exist with mass smaller than the hadronizing cluster mass, 
        *         false otherwise
        */
        bool setHeavyHadronSamplingProbabilities(int i_clusterCharmCharge,int i_clusterBottomCharge);

        /**
        * @brief Single hadron group sampling data calculation method
        *
        * @param i_hadronGroup Hadron group whose sampling probabilities are required
        * @param i_isLightGroup Group type flag: true if light flavored, false otherwise
        * @throw HadronizationException in case of error during hadron sampling data calculation (e.g., zero mean multiplicity)
        */
        void computeHadronSamplingWeights(const HadronGroup i_hadronGroup,bool i_isLightGroup);
    
        /**
        * @brief Not normalized light flavored hadron sampling probability calculation method
        *
        * @param i_hadronMass Hadron mass
        * @param i_hadronSpin Hadron spin multiplicity (2J+1)
        * @return Not normalized hadron sampling probability
        */
        double computeLightHadronMeanMultiplicity(double i_hadronMass,unsigned int i_hadronSpinMultiplicity) const;

        /**
        * @brief Not normalized heavy flavored hadron sampling probability calculation method
        *
        * @param i_hadronMass Hadron mass
        * @return Not normalized hadron sampling probability
        */
        double computeHeavyHadronSamplingWeight(double i_hadronMass) const;
    
        /**
        * @brief Light flavored hadron sampling method
        *
        * This is the high level light hadron group sampling method: using this method
        * the number of hadrons, belonging to each light flavor group, to be extracted
        * is generated. The single hadrons are then sampled making use of the 
        * runHadronSampling method
        */
        void runLightHadronSampling(void);
    
        /**
        * @brief Single group hadron sampling method
        *
        * Given a specific sampling group and a number of hadrons to be sampled from it, 
        * this method performs the sampling itself making use of the sampled data previously
        * computed and stored in the class data members
        *
        * @param i_hadronGroup Hadron sampling group
        * @param i_isHeavyGroup Group type flag: true if group is heavy flavored, false otherwise
        * @param i_numberOfHadrons Number of hadrons to be generated 
        */
        void runHadronSampling(HadronGroup i_hadronGroup,bool i_isHeavyGroup,unsigned int i_numberOfHadrons = 1);

        /**
        * @brief Channel sampling weight calculation method
        *
        * This method performs the calculation of the microcanonical weight corresponding 
        * to the sampled hadronization channel (except for phase space weight contribution 
        * and for the weight normalization externally obtained including the cluster partition 
        * function value). This calculation requires the knowledge of the cluster physical volume, 
        * which appears in the factor V/(2pi)^3 of the microcanonical weight formula, computed from 
        * the cluster mass and energy density provided as input parameter of the present method. 
        * It must be noted that the energy density is in this case the physical one, unlike the 
        * sampling energy density internally used by the present class for channel sampling only.
        *
        * @param i_clusterEnergyDensity Cluster (physical) energy density
        */
        void computeSampledChannelWeight(double i_clusterEnergyDensity);

        /**
        * Hadron sampling setup status flag
        */
        bool m_isHadronSamplingReady;
    
        /**
        * Sampled hadron set availability flag
        */
        bool m_isHadronSetAvailable;
    
        /**
        * Hadron sampling weight
        */
        double m_weight;
        
        /**
        * Hadron sampling groups
        */
        const HadronSamplingGroups* m_hadronSamplingGroups;
        
        /**
        * Random number generator
        */
        RandomNumberGenerator* m_randomGenerator;
        
        /**
        * Hadron sampling temperature parameter
        */
        double m_samplingTemperature;
        
        /**
        * Hadron sampling energy density parameter
        */
        double m_samplingEnergyDensity;
        
        /**
        * Hadron sampling volume parameter
        */
        double m_clusterSamplingVolume;
        
        /**
        * Input cluster mass
        */
        double m_clusterMass;
    
        /**
        * Input cluster residual strange charge
        */
        int m_residualStrangeCharge;
        
        /**
        * Input cluster residual electric charge
        */
        double m_residualElectricCharge;
        
        /**
        * Input cluster residual baryonic charge
        */
        double m_residualBaryonicCharge;
   
        /**
        * Single light hadron mean multiplicities weight factors (light flavored hadron only):
        * this data member stores the sampling weight contribution factor (2*J+1)/\hat{nu}_j. 
        * The storage of these quantities is convenient due to their constancy
        * wrt cluster properties: they are reused during each sampling weight calculation. 
        * The \hat{nu}_j mean grandcanonical multiplicity (excluding cluster volume factor) is stored 
        * in m_hadronMeanMultiplicities.
        */
        map<HadronGroup,vector<double> > m_lightHadronMultiplicityWeightFactors;
    
        /**
        * Single light hadron mean multiplicities weight factors (light flavored hadron only):
        * this data member stores the sampling weight contribution factor exp(sum(hat{nu}_j)),
        * where the sum is performed over light hadrons only. As for the data member above,
        * the storage of these quantities is convenient due to their constancy wrt cluster properties, 
        * so that they are reused during each sampling weight calculation.
        * The \hat{nu}_j mean grandcanonical multiplicity (excluding cluster volume factor) is stored
        * in m_hadronMeanMultiplicities.
        */
        double m_lightHadronExponentialWeightFactor;
    
        /**
        * Light hadron groups mean multiplicities (light flavored hadron only)
        */
        map<HadronGroup,double> m_hadronGroupMeanMultiplicities;
    
        /**
        * Single light hadron cumulative sampling probabilities (light and heavy flavored)
        */
        map<HadronGroup,vector<double> > m_samplingCumulativeProbabilities;
    
        /**
        * Heavy flavored hadron masses (Stored for faster probability calculation)
        */
        map<HadronGroup,vector<double> > m_heavyHadronMasses;

        /**
        * Single heavy flavored hadron sampling weight (not normalized)
        */
        map<HadronGroup,vector<double> > m_heavyHadronSamplingWeights;
    
        /**
        * Number of heavy flavored hadron accepted for the sampling procedure
        */
        unsigned int m_acceptedHeavyHadronNumber;

        /**
        * Single heavy flavored hadron sampling probabilities
        */
        vector<double> m_heavyHadronSamplingProbabilities;

        /**
        * Heavy flavored sampling group required for hadron sampling 
        * (If any, this is actively used only in case of heavy flavored cluster)
        */
        HadronGroup m_heavyHadronSamplingGroup;
        
        /**
        * Positions (wrt the corresponding hadron mass vector) of heavy flavored hadrons accepted for the sampling procedure
        */
        vector<unsigned int> m_availableHeavyHadronPositions;
    
        /**
        * Sum of sampled hadron mass
        */
        double m_hadronMassSum;

        /**
        * Sampled hadrons
        */
        vector<HadronData> m_hadronSet;
    
        /**
        * Sampled hadrons map, used for final sampling weight calculation
        */
        map<HadronGroup,vector<unsigned int> > m_hadronSetPositionMap;
    
    
        // TODO: test
        // hat{nu}_j
        void computeFugacities(double i_electricCharge,
                               int i_strangeCharge,
                               double i_baryonicCharge,
                               double& io_lambdaElectricCharge,
                               double& io_lambdaBaryonicCharge,
                               double& io_lambdaStrangeCharge,
                               double& io_samplingT);
        bool m_runNewCalculation;
        vector<double> m_lightHadronSpinMultiplicities;
        vector<double> m_lightHadronMasses;
        vector<double> m_lightHadronElectricCharges;
        vector<double> m_lightHadronBaryonicCharges;
        vector<double> m_lightHadronStrangeCharges;
        double m_residualMass;
        double m_electricCharge;
        double m_strangeCharge;
        double m_baryonicCharge;
        double m_lambdaElectricCharge;
        double m_lambdaStrangeCharge;
        double m_lambdaBaryonicCharge;
        double m_samplingT;
    
};

#endif