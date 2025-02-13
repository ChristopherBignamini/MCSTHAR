#ifndef PARTITIONFUNCTIONGENERATORSLAVE_H
#define PARTITIONFUNCTIONGENERATORSLAVE_H

#include "PartitionFunctionGenerator.h"
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include <vector>
#include <mpi.h>

using namespace std;

// TODO: check availability for code factorization with ClusterHadronization class

/**
* @brief Microcanonical single mass partition function calculation slave class
*
* This object is used to perform the calculation of the partition function
* data set corresponding to a single mass/energy density value (for each call to the run
* method) and for a set of strangeness suppression parameter values . See (// TODO: add ref)
* for details concerning the adopted calculation strategy.
*
* @author Christopher Bignamini
*/
// TODO: find better name for this class
class PartitionFunctionGeneratorSlave : public PartitionFunctionGenerator
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
        PartitionFunctionGeneratorSlave(double i_minGammaS,
                                        double i_maxGammaS,
                                        unsigned int i_numberOfGammaSValues,
                                        const HadronSamplingGroups& i_hadronSamplingGroups,
                                        unsigned int i_randomNumberGeneratorSeed,
                                        const MPI_Comm& i_mpiCommunicator);
        
        /**
        * @brief Destructor
        */
        ~PartitionFunctionGeneratorSlave(void);
            
    private:
        
        /**
        * @brief Private copy constructor
        */
        PartitionFunctionGeneratorSlave(const PartitionFunctionGeneratorSlave&);
        
        /**
        * @brief Private assignement operator
        */
        PartitionFunctionGeneratorSlave& operator=(const PartitionFunctionGeneratorSlave&);

        /**
        * @brief Single mass and energy density value partition function data set calculation method
        *
        * Internally used method for the calculation of the partition function data set
        * corresponding to the full set of strangeness suppression parameter values and to a
        * single mass/energy density value pair
        *
        * @param i_chargeConfiguration Charge configuration
        * @param i_mass Mass value
        * @param i_energyDensity Energy density value corresponding to the above index
        * @throw HadronizationException in case of error during calculation procedure, channel sampling
        *   or phase space integration (see data members for exceptions)
        */
        void runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
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
        void mergeLocalSums(void);
};

#endif