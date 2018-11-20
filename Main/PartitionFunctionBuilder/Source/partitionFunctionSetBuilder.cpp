#include "../../../Hadronization/HadronizationSetup/Include/PartitionFunctionCalculationSetupLoader.h"
#include "../../../Hadronization/HadronizationPartitionFunction/Include/PartitionFunctionGenerationHandlerMaster.h"
#include "../../../Hadronization/HadronizationPartitionFunction/Include/PartitionFunctionGenerationHandlerSlave.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/logMessage.h"
#include <iostream>
#include <sstream>
#include <mpi.h>

// TODO: add MPI exception/signal handling
/**
* @brief Partition function set calculation main function
*
* This is the main function used for partition function set calculation.
* The calculation configuration must be provided through an input configuration
* xml file (see the corresponding MCSTHARPartitionFunctionCalculationSetup.xsd 
* schema for details). The calculation results are stored in the output folder
* specified within the configuration file itself according to the structure and
* procedures detailed in the code documentation.
* The only available run configuration for the executable corresponding to the 
* present main is:
*
* mpirun -np <number_of_MPI_tasks> <path_to_partitionFunctionSetBuilder> <configuration_file_name> 
*
* In case of succesful execution the executable returns 0 while in case of error
* the exit code 1 is returned and the encountered error data is print to standard
* output. 
*
* @author Christopher Bignamini
*
* @param Partition function calculation setup file
* @return 0 in case of succesfull calculation run, 1 otherwise.
*/
int main(int argc, char* argv[])
{
    // TODO: check everywhere std::sqrt and std::abs
    try
    {
        // Initialize MPI environment
        int mpiExitStatus(MPI_Init(&argc,&argv));
        if(mpiExitStatus != 0)
        {
            // TODO: add exception with MPI error message
            return 1;
        }
        
        // Retrieve total task number and current task Id
        int taskId;
        int numberOfTasks;
        mpiExitStatus = MPI_Comm_rank(MPI_COMM_WORLD,
                                      &taskId);
        mpiExitStatus = MPI_Comm_size(MPI_COMM_WORLD,
                                      &numberOfTasks);
        if(mpiExitStatus != 0)
        {
            // TODO: add exception with MPI error message
            return 1;
        }
        
        // Log start run message
        if(taskId==0)
        {
            logMessage("Starting partition function set calculation");
        }
        
        // Check number of input arguments
        if(argc!=2)
        {
            // TODO: add exception handling for MPI
            throw HadronizationException("Error during partition function calculation run, wrong number of input parameters. Run as mpirun -np <number_of_MPI_tasks> <path_to_partitionFunctionSetBuilder> <configuration_file_name>",
                                         __FUNCTION__,801);
        }
        
        // Load partition function calculation setup file
        if(taskId==0)
        {
            logMessage("Loading calculation setup");
        }
        PartitionFunctionCalculationSetupLoader partitionFunctionCalculationSetup(argv[1]);
        
        // Build microcanonical partition function generation handler
        // (The single handler are being build in a MPI serial procedure in order to avoid same parallel file reading)
        // TODO: check problem above and find better procedure
        PartitionFunctionGenerationHandler* partitionFunctionGenerationHandler(NULL);
        const PartitionFunctionCalculationSetup& calculationSetup(
                                                                  partitionFunctionCalculationSetup.getPartitionFunctionCalculationSetup());
        for(int id=0;id<numberOfTasks;++id)
        {
            if(id==taskId)
            {
                if(id==0)
                {
                    logMessage("Building calculation handler");
                    
                    partitionFunctionGenerationHandler =
                    new PartitionFunctionGenerationHandlerMaster(calculationSetup,
                                                                 MPI_COMM_WORLD);                
                }
                else
                {
                    partitionFunctionGenerationHandler =
                    new PartitionFunctionGenerationHandlerSlave(calculationSetup.hadronDataSetFileName,
                                                                calculationSetup.lightHadronMaxMass,
                                                                calculationSetup.microcanonicalParameterGrid.minEnergyDensity,
                                                                calculationSetup.microcanonicalParameterGrid.maxEnergyDensity,
                                                                calculationSetup.microcanonicalParameterGrid.numberOfEnergyDensityValues,
                                                                calculationSetup.microcanonicalParameterGrid.minGammaS,
                                                                calculationSetup.microcanonicalParameterGrid.maxGammaS,
                                                                calculationSetup.microcanonicalParameterGrid.numberOfGammaSValues,
                                                                calculationSetup.massValueList,
                                                                calculationSetup.chargeConfiguration,
                                                                calculationSetup.randomNumberGeneratorSeed + taskId,// TODO: is this right?
                                                                MPI_COMM_WORLD);
                
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
                
        // Run generation
        if(taskId==0)
        {
            logMessage("Run calculation");
        }
        partitionFunctionGenerationHandler->run();
                
        // Free memory
        delete partitionFunctionGenerationHandler;
        
        // Log end run message
        // Run generation
        if(taskId==0)
        {
            logMessage("Partition function set calculation succesfully executed");
        }
        
        // Finalize MPI environment
        MPI_Finalize();


    }
    catch(HadronizationException& ex)
    {
        // TODO: MPI exception handling
        stringstream errCodeStr;
        errCodeStr<<ex.getReturnValue();
        const string errMessage("Error during partition function set calculation. " +
                                ex.getErrorMessage() +
                                "\nError code: " +
                                errCodeStr.str());
        logMessage(errMessage,ERROR);
        
        // TODO: find final solution
        MPI_Abort(MPI_COMM_WORLD,1);
                
        return 1;
    }
        
	return 0;
}

