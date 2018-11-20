#ifndef PARTITIONFUNCTIONGENERATIONHANDLERMASTER_H// TODO: check this practice
#define PARTITIONFUNCTIONGENERATIONHANDLERMASTER_H

#include "PartitionFunctionGenerationHandler.h"
#include "PartitionFunctionGeneratorMaster.h"
#include "../../HadronizationPartitionFunctionStorage/Include/PartitionFunctionArchiveWriter.h"
#include "../../HadronizationSetup/Include/PartitionFunctionCalculationSetup.h"
#include <mpi.h>

// TODO: add descritpion of parallelization strategy
/**
* @brief Microcanonical partition function calculation handling master class
*
* This class is used to handle the partition function calculation and
* output, in particular through the PartitionFunctionGenerator and
* PartitionFunctionArchiveWriter classes. According to the described 
* parallelization stratefy, this class..// TODO: add missing info 
*
* @author Christopher Bignamini
*/
class PartitionFunctionGenerationHandlerMaster : public PartitionFunctionGenerationHandler
{
    public:
        
        /**
        * @brief Constructor
        *
        * @param i_calculationSetup Partition funciton data set generation setup
        * @param i_mpiCommunicator MPI communicator
        * @throw HadronizationException in case of error in provided setup (see also data member constructor exceptions)
        */
        PartitionFunctionGenerationHandlerMaster(const PartitionFunctionCalculationSetup& i_calculationSetup,
                                                 const MPI_Comm& i_mpiCommunicator);
        
        /**
        * @brief Destructor
        */
        ~PartitionFunctionGenerationHandlerMaster(void);// TODO: remove empty destructor
        
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
        PartitionFunctionGenerationHandlerMaster(const PartitionFunctionGenerationHandlerMaster&);
        
        /**
        * @brief Private assignement (uncopiable class) // TODO: check this!
        */
        PartitionFunctionGenerationHandlerMaster& operator=(const PartitionFunctionGenerationHandlerMaster&);
                
        /**
        * Single mass value partition function data set generator (master)
        */
        PartitionFunctionGeneratorMaster m_partitionFunctionGenerator;
        
        /**
        * Partition function data set writer
        */
        PartitionFunctionArchiveWriter m_partitionFunctionWriter;
};

#endif
