#ifndef HADRONSAMPLINGNEW_H
#define HADRONSAMPLINGNEW_H

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
class HadronSamplingNew
{
    
    public:
    
        /**
        * @brief Constructor
        *
        * @param i_hadronSamplingGroups Hadron sampling group handler used for channel composition sampling
        * @param io_randomGenerator Random number generator
        * @throw HadronizationException in case of non positive sampling temperature or energy density 
        *        or errors during preliminary sampling data calculation
        */
        HadronSamplingNew(const HadronSamplingGroups& i_hadronSamplingGroups,
                          RandomNumberGenerator& io_randomGenerator);

        /**
        * @brief Destructor
        *
        */
        ~HadronSamplingNew(void);
    
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
    
        inline double getSamplingElectricChargeFugacity(void) const { return m_lambdaElectricCharge; }// TODO: why not using GC parameter structure
        inline double getSamplingStrangeChargeFugacity(void) const { return m_lambdaStrangeCharge; }// TODO: add check on parameter availability
        inline double getSamplingBaryonicChargeFugacity(void) const { return m_lambdaBaryonicCharge; }// TODO: cosa fare quando la stima non converge?
        inline double getSamplingCharmChargeFugacity(void) const { return m_lambdaCharmCharge; }
        inline double getSamplingBottomChargeFugacity(void) const { return m_lambdaBottomCharge; }
        inline double getSamplingTemperature(void) const { return m_samplingTemperature; }
        bool m_verboseLog;

    private:
    
    
        void retrieveHadronData(void);
        void computeFugacities(void);
        void computeSamplingData(void);
        void computeHadronSamplingWeights(const HadronGroup i_hadronGroup);
        double computeHadronMeanMultiplicity(const double i_hadronMass,const unsigned int i_hadronSpinMultiplicity) const;
        void runHadronSampling(void);
        void runHadronSampling(const HadronGroup i_hadronGroup,const unsigned int i_numberOfHadrons);
        void computeSampledChannelWeight(void);
    
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
        double m_energyDensity;
        
        /**
        * Input cluster mass
        */
        double m_clusterMass;

        /**
        * Hadron sampling volume parameter
        */
        double m_clusterSamplingVolume;
        

        map<HadronGroup,double> m_hadronGroupMeanMultiplicities;
        
        map<HadronGroup,vector<double> > m_samplingCumulativeProbabilities;
    
        map<HadronGroup,vector<double> > m_hadronMultiplicityWeightFactors;
    
        double m_hadronExponentialWeightFactor;
    
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
                               int i_charmCharge,
                               int i_bottomCharge,
                               double& io_lambdaElectricCharge,
                               double& io_lambdaBaryonicCharge,
                               double& io_lambdaStrangeCharge,
                               double& io_lambdaCharmCharge,
                               double& io_lambdaBottomCharge,
                               double& io_samplingT);
        vector<double> m_hadronSpinMultiplicities;
        vector<double> m_hadronMasses;
        vector<double> m_hadronElectricCharges;
        vector<double> m_hadronBaryonicCharges;
        vector<double> m_hadronStrangeCharges;
        vector<double> m_hadronCharmCharges;
        vector<double> m_hadronBottomCharges;
        double m_electricCharge;
        double m_strangeCharge;
        double m_baryonicCharge;
        double m_charmCharge;
        double m_bottomCharge;
        double m_residualElectricCharge;
        double m_residualStrangeCharge;
        double m_residualBaryonicCharge;
        double m_residualCharmCharge;
        double m_residualBottomCharge;
        double m_lambdaElectricCharge;
        double m_lambdaStrangeCharge;
        double m_lambdaBaryonicCharge;
        double m_lambdaCharmCharge;
        double m_lambdaBottomCharge;
};

#endif