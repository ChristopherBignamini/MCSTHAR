#ifndef PARTITIONFUNCTIONGENERATORMASTER_H
#define PARTITIONFUNCTIONGENERATORMASTER_H

// TODO: the not homogeneous include list which follows is a sign of bad design...
#include "PartitionFunctionGenerator.h"
#include "PartitionFunctionGenerationData.h"
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include <vector>
#include <mpi.h>

using namespace std;

// TODO: check availability for code factorization with ClusterHadronization class

/**
* @brief Microcanonical single mass partition function calculation master class
*
* This object is used to perform the calculation of the partition function
* data set corresponding to a single mass/energy density value (for each call to the run
* method) and for a set of strangeness suppression parameter values . See (// TODO: add ref)
* for details concerning the adopted calculation strategy.
*
* @author Christopher Bignamini
*/
// TODO: find better name for this class
class PartitionFunctionGeneratorMaster : public PartitionFunctionGenerator
{
    public:
        
        // TODO: the sampling parameters are strongly linked with partition function ones, add compatibility check and
        //       describe difference between phenomenological and sampling parameters
        /**
        * @brief Constructor
        *
        * @param i_minGammaS Minimum strangeness suppression parameter value
        * @param i_maxGammaS Maximum strangeness suppression parameter value
        * @param i_numberOfGammaSValues Number of strangeness suppression parameter values
        * @param i_hadronSamplingGroups Hadronization channel sampling group handler
        * @param i_maxChannelSamplingNumber Maximum number of allowed channel samplings for partition function calculation
        * @param i_integrationErrorThreshold Integration stability threshold (relative error) for partition function calculation
        * @param i_randomNumberGeneratorSeed Random number generator seed
        * @param i_mpiCommunicator MPI communicator
        * @throw HadronizationException in case of error in hadron sampling or phase space integration objects initialization and
        *   input parameter values
        */
        PartitionFunctionGeneratorMaster(double i_minGammaS,
                                         double i_maxGammaS,
                                         unsigned int i_numberOfGammaSValues,
                                         const HadronSamplingGroups& i_hadronSamplingGroups,
                                         unsigned long long int i_maxChannelSamplingNumber,
                                         double i_integrationErrorThreshold,
                                         unsigned int i_randomNumberGeneratorSeed,
                                         const MPI_Comm& i_mpiCommunicator);
    
        /**
        * @brief Destructor
        */
        ~PartitionFunctionGeneratorMaster(void);
    
        /**
        * @brief Single mass and energy density partition function data get method
        *
        * @return Single mass and energy density partition function data
        * @throw HadronizationException in case of partition function data unavailability
        */
        const SingleMassEnergyDensityPairPartitionFunctionData& getPartitionFunctionData(void) const;
        
    private:
        
        /**
        * @brief Private copy constructor
        */
        PartitionFunctionGeneratorMaster(const PartitionFunctionGeneratorMaster&);
        
        /**
         * @brief Private assignement operator
         */
        PartitionFunctionGeneratorMaster& operator=(const PartitionFunctionGeneratorMaster&);
        
        /**
        * @brief Single mass and energy density value partition function data set calculation method
        *
        * Internally used method for the calculation of the partition function data set
        * corresponding to the full set of strangeness suppression parameter values and to a
        * single mass/energy density value pair
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_mass Mass value
        * @param i_energyDensity Energy density value
        * @throw HadronizationException in case of error during calculation procedure, channel sampling
        *   or phase space integration (see data members for exceptions)
        */
        void runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                             double i_mass,
                                             double i_energyDensity);
                
        /**
        * @brief Partition function data reset method
        *
        * Internally used method for single mass partition function data reset, used in case
        * of hadronization channel unavailability only
        */
        void resetPartitionFunctionData(void);
    
        /**
        * @brief Local calculation data merging method
        *
        * Internally used method for the communication to partition function generator
        * master class of the local calculation data and for the subsequent data merging.
        *
        * @throw HadronizationException in case of error during data communication/reduction
        */
        // TODO: move to base class using task id + if???
        void mergeLocalSums(void);
    
        /**
        * Partition function data availability flag
        */
        bool m_isPartitionFunctionAvailable;
        
        /**
        * Single mass and energy density partition function data
        */
        SingleMassEnergyDensityPairPartitionFunctionData m_partitionFunctionData;// TODO: this data member leads to data duplication
    
        /**
        * Partition function calculation global data
        *
        * This data member store the (cumulative) global calculationdata of all single.
        * task. Data is stored according to the structure described for m_weightSumData 
        * in base class.
        * 
        */
        double* m_cumulativeWeightSumData;// TODO: switch to MPI derived data types

        /**
        * Maximum number of allowed hadronization channel samplings
        */
        unsigned long long int m_totalMaxSamplingNumber;

        /**
        * Single task maximum number of (full) channel samplings
        */
        unsigned long long int m_localMaxSamplingNumber;

        /**
        * Minimum channel sampling to be performed before the first stability check (to avoid spurious stability)
        */
        // TODO: static const or not?  Add input parameter?
        unsigned long long int m_minimumStabilityCheckSamplingNumber;// TODO: use static or not?
    
        /**
        * Global (phase space integration and sum over hadronization channels) integration stability squared (relative) error threshold
        */
        const double m_integrationSquaredErrorThreshold;
    
        /**
        * Number of MPI tasks
        */
        int m_numberOfTasks;
};

#endif