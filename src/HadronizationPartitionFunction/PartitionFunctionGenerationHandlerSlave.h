#ifndef PARTITIONFUNCTIONGENERATIONHANDLERSLAVE_H// TODO: check this practice
#define PARTITIONFUNCTIONGENERATIONHANDLERSLAVE_H

#include "PartitionFunctionGenerationHandler.h"
#include "PartitionFunctionGeneratorSlave.h"
#include <vector>
#include <string>
#include <mpi.h>

using namespace std;

// TODO: add descritpion of parallelization strategy
/**
* @brief Microcanonical partition function calculation handling slave class
*
* This class is used to handle the partition function calculation: according 
* to the described parallelization strategy, this class..// TODO: add missing info
*
* @author Christopher Bignamini
*/
class PartitionFunctionGenerationHandlerSlave : public PartitionFunctionGenerationHandler
{
    public:
    
        /**
        * @brief Constructor
        *
        * @param i_hadronDataSetFileName Hadron list file
        * @param i_lightHadronMaxMass Maximum light flavored hadron mass
        * @param i_minEnergyDensity Minimum energy density parameter value
        * @param i_maxEnergyDensity Maximum energy density parameter value
        * @param i_numberOfEnergyDensityValues Number of values in energy density grid
        * @param i_minGammaS Minimum strangeness suppression parameter value
        * @param i_maxGammaS Maximum strangeness suppression parameter value
        * @param i_numberOfGammaSValues Number of values in strangeness suppression parameter grid
        * @param i_massValueList Mass values list
        * @param i_chargeConfiguration Abelian charge configuration
        * @param i_randomNumberGeneratorSeed Random number generator seed
        * @param i_mpiCommunicator MPI communicator
        * @throw HadronizationException in case of error in provided setup (see also data member constructor exceptions)
        */
        PartitionFunctionGenerationHandlerSlave(const string& i_hadronDataSetFileName,
                                                double i_lightHadronMaxMass,
                                                double i_minEnergyDensity,
                                                double i_maxEnergyDensity,
                                                double i_numberOfEnergyDensityValues,
                                                double i_minGammaS,
                                                double i_maxGammaS,
                                                unsigned int i_numberOfGammaSValues,
                                                const vector<double>& i_massValueList,
                                                const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                unsigned int i_randomNumberGeneratorSeed,
                                                const MPI_Comm& i_mpiCommunicator);
    
        /**
        * @brief Destructor
        */
        ~PartitionFunctionGenerationHandlerSlave(void);// TODO: remove empty destructor
        
        /**
        * @brief Generation run method
        *
        * A call to this method starts the generation procedure
        *
        * @throw HadronizationException in case of error during generation procedure or handler initialization
        */
        void run(void);
        
    private:
        
        /**
        * @brief Private copy constructor (uncopiable class) // TODO: check this!
        */
        PartitionFunctionGenerationHandlerSlave(const PartitionFunctionGenerationHandlerSlave&);
        
        /**
        * @brief Private assignement (uncopiable class) // TODO: check this!
        */
        PartitionFunctionGenerationHandlerSlave& operator=(const PartitionFunctionGenerationHandlerSlave&);
    
        /**
        * Single mass value partition function data set generator (slave)
        */
        PartitionFunctionGeneratorSlave m_partitionFunctionGenerator;
};

#endif
