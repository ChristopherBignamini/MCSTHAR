#ifndef PARTITIONFUNCTIONGENERATOROLD_H
#define PARTITIONFUNCTIONGENERATOROLD_H

// TODO: the not homogeneous include list which follows is a sign of bad design...
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include "../../HadronizationSetup/Include/MicrocanonicalParameterGrid.h"
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../HadronizationChannelGenerator/Include/HadronSampling.h"
#include "PartitionFunctionGenerationData.h"
#include "PhaseSpaceIntegrator.h"
#include "PhaseSpaceIntegrationData.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include <sstream>
#include <set>
#include <vector>

// TODO: check availability for code factorization with ClusterHadronization class

/**
* @brief Partition function channel ordering operator (based on particle Id list)
*
* @author Christopher Bignamini
*/
class OrderedHadronIdsComparison
{
    public:
    
        /**
        * @brief Partition function channel ordering operator (based on particle Id list)
        *
        * @param i_firstChannelOrderedHadronIds Increasingly ordered list of hadron Ids of first channel
        * @param i_secondChannelOrderedHadronIds Increasingly ordered list of hadron Ids of second channel
        * @return True if first channel is minor (check implementation) than the second, false otherwise 
        */
        bool operator()(const vector<int>& i_firstChannelOrderedHadronIds,
                        const vector<int>& i_secondChannelOrderedHadronIds) const;
};
    
/**
* @brief Microcanonical single mass partition function calculation class
*
* This object is used to perform the calculation of the partition function
* data set corresponding to a single mass value (for each call to the run
* method) and to the microcanonical parameter grid points. See (// TODO: add ref)
* for details concerning the adopted calculation strategy.
* 
* @author Christopher Bignamini
*/
// TODO: find better name for this class
class PartitionFunctionGeneratorOld
{
    public:
    
        // TODO: the sampling parameters are strongly linked with partition function ones, add compatibility check and
        //       describe difference between phenomenological and sampling parameters
        /**
        * @brief Constructor
        *
        * @param i_microcanonicalParameterGrid Microcanonical parameter grid structure
        * @param i_hadronSamplingTemperature Hadron sampling probailities temperature parameter
        * @param i_hadronSamplingGroups Hadronization channel sampling group handler
        * @param i_maxChannelSamplingNumber Maximum number of allowed channel samplings for partition function calculation
        * @param i_integrationErrorThreshold Integration stability threshold (relative error) for partition function calculation
        * @param i_maxPhaseSpaceSamplingNumber Maximum number of allowed phase space point samplings for phase space integration
        * @param i_phaseSpaceIntegrationErrorThreshold Integration stability threshold (relative error) for phase space integration
        * @param io_randomGenerator Random number generator (I/O parameter)
        * @throw HadronizationException in case of error in hadron sampling or phase space integration objects initialization
        */
        // TODO: move i_hadronSamplingTemperature somewhere else
        PartitionFunctionGeneratorOld(const MCSTHAR::Setup::MicrocanonicalParameterGrid& i_microcanonicalParameterGrid,
                                      double i_hadronSamplingTemperature,
                                      const HadronSamplingGroups& i_hadronSamplingGroups,
                                      unsigned long int i_maxChannelSamplingNumber,
                                      double i_integrationErrorThreshold,
                                      unsigned long int i_maxPhaseSpaceSamplingNumber,
                                      double i_phaseSpaceIntegrationErrorThreshold,
                                      RandomNumberGenerator& io_randomGenerator);

        /**
        * @brief Destructor
        */
        ~PartitionFunctionGeneratorOld(void);
        
        // Default copy constructor and assignement operator are being used //TODO: check
    
        /**
        * @brief Single mass value partition function data set calculation method
        *
        * A call to this method starts the calculation procedure
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_mass Mass value
        * @throw HadronizationException in case of error during calculation procedure, generator initialization or non positive mass provided
        */
        void run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                 double i_mass);
        
        /**
        * @brief Single mass partition function data get method
        *
        * @return Single mass partition function data
        * @throw HadronizationException in case of partition function data unavailability
        */
        const PartitionFunctionData& getPartitionFunctionData(void) const;
    
    
    private:

        /**
        * @brief Single mass and energy density value partition function data set calculation method
        *
        * Internally used method for the calculation of the partition function data set 
        * corresponding to the full set of strangeness suppression parameter values and to a 
        * single mass/energy density value pair 
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_mass Mass value
        * @param i_energyDensityIndex Energy density element index in energy density grid
        * @param i_energyDensity Energy density value corresponding to the above index
        * @throw HadronizationException in case of error during calculation procedure, channel sampling 
        *   or phase space integration (see data members for exceptions)
        */
        void runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                             double i_mass,
                                             unsigned int i_energyDensityIndex,
                                             double i_energyDensity);

        // TODO: remove after debugging/testing
        void computeStrangenessSuppression(vector<double>& io_strangenessSuppressionFactors);
    
        /**
        * @brief Partition function data reset method
        *
        * Internally used method for single mass partition function data reset, used in case
        * of hadronization channel unavailability only 
        */
        void resetPartitionFunctionData(void);
    
        /**
        * Partition function generator initialization status flag
        */
        bool m_isPartitionFunctionGeneratorReady;
        
        /**
        * Partition function data availability flag
        */
        bool m_isPartitionFunctionAvailable;
        
        /**
        * Single mass partition function data
        */
        PartitionFunctionData m_partitionFunctionData;

        /**
        * Hadron sampling handler
        */
        HadronSampling m_hadronSampling;
            
        /**
        * Phase space integrator
        */
        PhaseSpaceIntegrator m_phaseSpaceIntegrator;
    
        /**
        * Sampled channel phase space integration data set:
        * map key is the increasingly ordered list of hadron Ids for a given channel
        * and the value the corresponding phase space integration data. The usage of
        * this data member allows to avoid the phase space integral recalculation
        * for those channel sampled more than once during the calculation. 
        */
        map<vector<int>,PhaseSpaceIntegrationData,OrderedHadronIdsComparison> m_channelPhaseSpaceData;
    
        /**
        * Maximum number of allowed hadronization channel samplings
        */
        unsigned long int m_maxChannelSamplingNumber;
        
        /**
        * Global (phase space integration and sum over hadronization channels) integration stability squared (relative) error threshold
        */
        const double m_integrationSquaredErrorThreshold;
    
        /**
        * Channel availability status flag
        */
        bool m_noChannelAvailability;
    
        /**
        * Previous run last energy density element index (used for hadron sampling object setting)
        */
        unsigned int m_previousRunLastEnergyDensityIndex;
    
        /**
        * Number of energy density values (string)
        */
        stringstream m_numberOfEnergyDensityValueStr;
    
        /**
        * Minimum channel sampling to be performed before the first stability check (to avoid spurious stability)
        */
        // TODO: static const or not?  Add input parameter?
        unsigned long int m_minimumStabilityCheckSamplingNumber;// TODO: use static or not?
};

#endif