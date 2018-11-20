#ifndef PARTITIONFUNCTIONGENERATOR_H
#define PARTITIONFUNCTIONGENERATOR_H

// TODO: the not homogeneous include list which follows is a sign of bad design...
#include "PartitionFunctionPhaseSpaceSampling.h"
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../HadronizationChannelGenerator/Include/HadronSampling.h"
#include "../../HadronizationChannelGenerator/Include/HadronSamplingNew.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include <vector>
#include <mpi.h>

using namespace std;

// TODO: check availability for code factorization with ClusterHadronization class
// TODO: add description of parallelization strategy
/**
* @brief Microcanonical single mass partition function calculation base class
*
* This object is used to perform the calculation of the partition function
* data set corresponding to a single mass/energy density value (for each call to the run
* method) and for a set of strangeness suppression parameter values . See (// TODO: add ref)
* for details concerning the adopted calculation strategy. This class provides the base
* interface and data members for the derived (MPI related) classes PartitionFunctionGeneratorMaster
* and PartitionFunctionGeneratorSlave
*
* @author Christopher Bignamini
*/
// TODO: find better name for this class
class PartitionFunctionGenerator
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
        * @param i_randomNumberGeneratorSeed Random number generator seed
        * @param i_mpiCommunicator MPI communicator
        * @throw HadronizationException in case of error in hadron sampling or phase space integration objects initialization and
        *   input parameter values
        */
        PartitionFunctionGenerator(double i_minGammaS,
                                   double i_maxGammaS,
                                   unsigned int i_numberOfGammaSValues,
                                   const HadronSamplingGroups& i_hadronSamplingGroups,
                                   unsigned int i_randomNumberGeneratorSeed,
                                   const MPI_Comm& i_mpiCommunicator);
        
        /**
        * @brief Destructor
        */
        virtual ~PartitionFunctionGenerator(void);
        
        /**
        * @brief Single mass and energy density value partition function data set calculation method
        *
        * A call to this method starts the calculation procedure of the partition function data set
        * corresponding to the full set of strangeness suppression parameter values and to the provided
        * mass/energy density value pair
        *
        * @param i_chargeConfiguration Abelian charge configuration
        * @param i_mass Mass value
        * @param i_energyDensity Energy density value
        * @throw HadronizationException in case of error during calculation procedure, generator initialization,
        *   channel sampling or phase space integration (see data members for exceptions) and in case of non positive
        *   mass/energy density provided
        */
        void run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                 double i_mass,
                 double i_energyDensity);
    
    protected:
        
        /**
        * @brief Private copy constructor
        */
        PartitionFunctionGenerator(const PartitionFunctionGenerator&);
        
        /**
        * @brief Private assignement operator
        */
        PartitionFunctionGenerator& operator=(const PartitionFunctionGenerator&);
    
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
        virtual void runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                     double i_mass,
                                                     double i_energyDensity) = 0;

        // TODO: check everywhere the exception case list
        /**
        * @brief Interation block calculation method
        *
        * Internally used method for the calculation result update through the execution of a new
        * set of calculation iterations. After each new iteration block the calculation 
        * stability status is being checked by the PartitionFunctionGeneratorMaster class
        *  
        * @param i_chargeConfiguration Charge configuration
        * @param i_mass Mass value
        * @param i_energyDensity Energy density value
        * @throw HadronizationException in case of error during calculation procedure, channel sampling
        *   or phase space integration (see data members for exceptions)
        */
        void runNewCalculationIterations(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                         double i_mass,
                                         double i_energyDensity);
    
        /**
        * @brief Local calculation data merging method
        *
        * Internally used method for the communication to partition function generator 
        * master class of the local calculation data and for the subsequent data merging.
        *
        * @throw HadronizationException in case of error during data communication/reduction
        */
        // TODO: move to base class using task id + if???
        virtual void mergeLocalSums(void) = 0;
    
        /**
        * @brief Strangeness suppression factor calculation method
        *
        * Internally used method for the calculation of the strangeness suppression factor
        * corresponding to the full set of strangeness suppression parameter values and for
        * the hadronization channel whose strange composition is given by the method input
        * parameters 
        * 
        * @param i_ssbarComponentValues ssbar wave function components of the channel hadrons
        * @param i_oneMinusSSBARComponentValues 1-ssbar wave function components of the channel hadrons
        * @param i_numbeOfStrangeQuarks Number of strange quarks of the channel hadrons
        */
        void computeStrangenessSuppression(const vector<double>& i_ssbarComponentValues,
                                           const vector<double>& i_oneMinusSSBARComponentValues,
                                           const vector<unsigned int>& i_numbeOfStrangeQuarks);
                
        /**
        * Partition function generator initialization status flag
        */
        bool m_isPartitionFunctionGeneratorReady;
        
        /**
        * Strangeness suppression parameter values
        */
        vector<double> m_gammaSValues;
        
        /**
        * Number of strangeness suppression parameter values
        */
        unsigned int m_numberOfGammaSValues;
    
        /**
        * Single task stability check number of iterations
        */
        unsigned long long int m_localStabilityCheckSamplingNumber;

        /**
        * Single task stability check number of iterations 
        * (this sampling number step is used until max number of sampling is reached)
        */
        unsigned long long int m_localStabilityCheckSamplingNumberStart;

        /**
        * Random number generator
        */
        RandomNumberGenerator* m_randomGenerator;
        
        /**
        * Hadron sampling handler
        */
//        HadronSampling m_hadronSampling;
        HadronSamplingNew m_hadronSamplingNew;
    
        /**
        * Phase space sampling handler
        */
        PartitionFunctionPhaseSpaceSampling m_phaseSpaceSampling;
        
        /**
        * Partition function calculation data
        *
        * Weigth data vector contains serialzed partition function calculation data as
        * sequence of weight and weight^2 sum pairs, each pair corresponding to a single 
        * strangeness suppression parameter value, followed by a last status element set 
        * to 0 if no error is found during calculation (and 1 otherwise). This particular 
        * data structure is used to allow data communication with MPI without additional 
        * memory copying.
        *  
        */
        double* m_weightSumData;// TODO: switch to MPI derived data types

        /**
        * Partition function calculation data size (number of elements in m_weightSumData array)
        */
        const unsigned int m_weightSumDataSize;
    
        /**
        * Strangeness suppression factor value storage vector
        */
        vector<double> m_strangenessSuppressionFactors;
    
        /**
        * MPI communicator
        */
        const MPI_Comm& m_mpiCommunicator;        
};

#endif