#include "../Include/PartitionFunctionGenerationHandlerMaster.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/logMessage.h"
#include <sstream>
#include <vector>

using namespace std;

PartitionFunctionGenerationHandlerMaster::PartitionFunctionGenerationHandlerMaster(const PartitionFunctionCalculationSetup& i_calculationSetup,
                                                                                   const MPI_Comm& i_mpiCommunicator)
    :PartitionFunctionGenerationHandler(i_calculationSetup.hadronDataSetFileName,
                                        i_calculationSetup.lightHadronMaxMass,
                                        i_calculationSetup.microcanonicalParameterGrid.minEnergyDensity,
                                        i_calculationSetup.microcanonicalParameterGrid.maxEnergyDensity,
                                        i_calculationSetup.microcanonicalParameterGrid.numberOfEnergyDensityValues,
                                        i_calculationSetup.massValueList,
                                        i_calculationSetup.chargeConfiguration)
    ,m_partitionFunctionGenerator(i_calculationSetup.microcanonicalParameterGrid.minGammaS,
                                  i_calculationSetup.microcanonicalParameterGrid.maxGammaS,
                                  i_calculationSetup.microcanonicalParameterGrid.numberOfGammaSValues,
                                  m_hadronSamplingGroups,
                                  i_calculationSetup.maxChannelSamplingNumber,
                                  i_calculationSetup.integrationErrorThreshold,
                                  i_calculationSetup.randomNumberGeneratorSeed,
                                  i_mpiCommunicator)
    ,m_partitionFunctionWriter(i_calculationSetup.outputFolder,
                               MicrocanonicalParameterGridStructure(i_calculationSetup.microcanonicalParameterGrid.minEnergyDensity,
                                                                    i_calculationSetup.microcanonicalParameterGrid.maxEnergyDensity,
                                                                    i_calculationSetup.microcanonicalParameterGrid.numberOfEnergyDensityValues,
                                                                    i_calculationSetup.microcanonicalParameterGrid.minGammaS,
                                                                    i_calculationSetup.microcanonicalParameterGrid.maxGammaS,
                                                                    i_calculationSetup.microcanonicalParameterGrid.numberOfGammaSValues),
                               m_chargeConfiguration,
                               m_massValues)
{
    m_isHandlerAvailable = true;
}

// TODO: all these empty destructors can be avoided!
PartitionFunctionGenerationHandlerMaster::~PartitionFunctionGenerationHandlerMaster(void)
{
}

void PartitionFunctionGenerationHandlerMaster::run(void)
{
    if(m_isHandlerAvailable==false)
    {
        throw HadronizationException("Error during partition function calculation run, calculation handler not available",
                                     __FUNCTION__,881);
    }
    
    // Loop over mass and energy density values
    stringstream numberOfMassesStr;
    numberOfMassesStr<<m_massValues.size();
    stringstream numberOfEnergyDensitiesStr;
    numberOfEnergyDensitiesStr<<m_energyDensityValues.size();
    vector<double>::const_iterator massValueIt(m_massValues.begin());
    for(unsigned int massIndex=0;massIndex<m_massValues.size();++massIndex,++massValueIt)
    {
        // Log status message
        // Create log message
        // TODO: store log string as const static (?);
        stringstream massIndexStr;
        massIndexStr<<(massIndex+1);
        // TODO: use function for message pattern creation
        const string message("Running calculation for mass value " +
                             massIndexStr.str() +
                             " of " +
                             numberOfMassesStr.str());
        // Log message
        logMessage(message);
        
        // Loop over energy density grid points
        vector<double>::const_iterator energyValueIt(m_energyDensityValues.begin());
        for(unsigned int energyDensityIndex=0;energyDensityIndex<m_energyDensityValues.size();
            ++energyDensityIndex,++energyValueIt)
        {
            // Create log message
            // TODO: store log string as const static (?);
            stringstream energyDensityIndexStr;
            energyDensityIndexStr<<(energyDensityIndex+1);
            // TODO: use function for message pattern creation
            const string message("Running calculation for energy density value number " +
                                 energyDensityIndexStr.str() +
                                 " of " +
                                 numberOfEnergyDensitiesStr.str());
            // Log message
            logMessage(message);
            
            // TODO: check memory usage, it could be useful to print data during run (asincronously)
            try
            {
                // Run partition function set generation
                m_partitionFunctionGenerator.run(m_chargeConfiguration,
                                                 *massValueIt,
                                                 *energyValueIt);
            }
            catch(HadronizationException& ex)
            {
                throw ex;
            }
            
            // Log status message and write partition function data
            // Log message
            logMessage("Writing single mass/energy density value partition function set");
            try
            {
                // Write partition fuction data for current mass/energy density value
                m_partitionFunctionWriter.updatePartitionFunctionData(m_partitionFunctionGenerator.getPartitionFunctionData(),
                                                                      massIndex);
            }
            catch(HadronizationException& ex)
            {
                throw ex;
            }
        }
    }
    return;
}
