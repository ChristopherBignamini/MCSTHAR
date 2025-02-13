#include "../Include/PartitionFunctionGeneratorSlave.h"
#include "../../../Utilities/Include/HadronizationException.h"

// TODO: check everywhere exit code
// TODO: check everywhere input parameter constness

PartitionFunctionGeneratorSlave::PartitionFunctionGeneratorSlave(const double i_minGammaS,
                                                                 const double i_maxGammaS,
                                                                 const unsigned int i_numberOfGammaSValues,
                                                                 const HadronSamplingGroups& i_hadronSamplingGroups,
                                                                 const unsigned int i_randomNumberGeneratorSeed,
                                                                 const MPI_Comm& i_mpiCommunicator)
                                                                :PartitionFunctionGenerator(i_minGammaS,
                                                                                            i_maxGammaS,
                                                                                            i_numberOfGammaSValues,
                                                                                            i_hadronSamplingGroups,
                                                                                            i_randomNumberGeneratorSeed,
                                                                                            i_mpiCommunicator)
{
    // Receive from master task stability check sampling number data
    int exitStatus(MPI_Bcast(&m_localStabilityCheckSamplingNumberStart,
                             1,
                             MPI_UNSIGNED_LONG_LONG,
                             0,
                             m_mpiCommunicator));
    
    // Set data members
    if(exitStatus==0)
    {
        m_isPartitionFunctionGeneratorReady = true;
    }
    else
    {
        // TODO: add MPI error status
        throw HadronizationException("Error during stability check data broadcasting",
                                     __FUNCTION__,
                                     851);
    }
}

// TODO: remove empty destructors
PartitionFunctionGeneratorSlave::~PartitionFunctionGeneratorSlave(void)
{
}

// TODO: remove the energy density input parameter with the new m_hadronSampling interface
// TODO: refactoring required
void PartitionFunctionGeneratorSlave::runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                                      const double i_mass,
                                                                      const double i_energyDensity)
{
    // TODO: unsigned long int or unsigned long long int? Check counters everywhere
    // TODO: Is there a better calculation method to avoid huge numbers?
    // TODO: check alternative strategy (no communication)
    bool newIterationsRequired(true);
    while(newIterationsRequired)
    {
        try
        {
            // Run calculation for m_localStabilityCheckSamplingNumberStep iterations
            runNewCalculationIterations(i_chargeConfiguration,
                                        i_mass,
                                        i_energyDensity);
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
        
        // Slave tasks receive calculation status check flag
        const int exitStatus(MPI_Bcast(&m_localStabilityCheckSamplingNumber,// TODO: switch to data member?
                                       1,
                                       MPI_UNSIGNED_LONG_LONG,
                                       0,
                                       m_mpiCommunicator));
        
        if(exitStatus!=0)
        {
            // TODO: add MPI error status
            throw HadronizationException("Error during data data broadcasting",
                                         __FUNCTION__,
                                         852);
        }
        if(m_localStabilityCheckSamplingNumber!=0)
        {
            newIterationsRequired = true;
        }
        else
        {
            newIterationsRequired = false;
        }
        
    }
    
    return;
}

void PartitionFunctionGeneratorSlave::mergeLocalSums(void)
{
    
    int exitStatus(MPI_Reduce(m_weightSumData,
                              m_weightSumData,
                              m_weightSumDataSize,
                              MPI_DOUBLE,
                              MPI_SUM,
                              0,
                              m_mpiCommunicator));
    if(exitStatus!=0)
    {
        // TODO: add MPI error status
        throw HadronizationException("Error during data calculation reduction",
                                     __FUNCTION__,
                                     853);
    }
    
    return;
}