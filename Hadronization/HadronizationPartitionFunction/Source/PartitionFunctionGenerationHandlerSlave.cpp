#include "../Include/PartitionFunctionGenerationHandlerSlave.h"
#include "../../../Utilities/Include/HadronizationException.h"

PartitionFunctionGenerationHandlerSlave::PartitionFunctionGenerationHandlerSlave(const string& i_hadronDataSetFileName,
                                                                                 const double i_lightHadronMaxMass,
                                                                                 const double i_minEnergyDensity,
                                                                                 const double i_maxEnergyDensity,
                                                                                 const double i_numberOfEnergyDensityValues,
                                                                                 double i_minGammaS,
                                                                                 double i_maxGammaS,
                                                                                 unsigned int i_numberOfGammaSValues,
                                                                                 const vector<double>& i_massValueList,
                                                                                 const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                                                 const unsigned int i_randomNumberGeneratorSeed,
                                                                                 const MPI_Comm& i_mpiCommunicator)
                                                                        :PartitionFunctionGenerationHandler(i_hadronDataSetFileName,
                                                                                                            i_lightHadronMaxMass,
                                                                                                            i_minEnergyDensity,
                                                                                                            i_maxEnergyDensity,
                                                                                                            i_numberOfEnergyDensityValues,
                                                                                                            i_massValueList,
                                                                                                            i_chargeConfiguration)
                                                                        ,m_partitionFunctionGenerator(i_minGammaS,
                                                                                                      i_maxGammaS,
                                                                                                      i_numberOfGammaSValues,
                                                                                                      m_hadronSamplingGroups,
                                                                                                      i_randomNumberGeneratorSeed,
                                                                                                      i_mpiCommunicator)
{
    m_isHandlerAvailable = true;
}

// TODO: all these empty destructors can be avoided!
PartitionFunctionGenerationHandlerSlave::~PartitionFunctionGenerationHandlerSlave(void)
{
}

void PartitionFunctionGenerationHandlerSlave::run(void)
{
    if(m_isHandlerAvailable==false)
    {
        throw HadronizationException("Error during partition function calculation run, calculation handler not available",
                                     __FUNCTION__,871);
    }
    
    // Loop over mass and energy density values
    for(unsigned int massIndex=0;massIndex<m_massValues.size();++massIndex)
    {        
        // Loop over energy density grid points
        for(unsigned int energyDensityIndex=0;energyDensityIndex<m_energyDensityValues.size();
            ++energyDensityIndex)
        {
            // TODO: check memory usage, it could be useful to print data during run (asincronously)
            try
            {
                // Run partition function set generation
                m_partitionFunctionGenerator.run(m_chargeConfiguration,
                                                 m_massValues.at(massIndex),
                                                 m_energyDensityValues.at(energyDensityIndex));
            }
            catch(HadronizationException& ex)
            {
                throw ex;
            }
        }
    }
    return;
}
