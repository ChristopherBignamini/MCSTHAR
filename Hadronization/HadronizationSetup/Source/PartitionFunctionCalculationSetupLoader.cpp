#include "../Include/MCSTHARPartitionFunctionCalculationSetup.h"
#include "../Include/PartitionFunctionCalculationSetupLoader.h"
#include "../../../Utilities/Include/HadronizationException.h"


// TODO: move check for input parameter in other classes (also for hadronization setup),
//       so that the check is performed even in case of different class usage
PartitionFunctionCalculationSetupLoader::PartitionFunctionCalculationSetupLoader(const string& i_partitionFunctionCalculationSetupFileName)
                                                                                :m_partitionFunctionCalculationSetupAvailability(false)
{
    auto_ptr< ::MCSTHARPartitionFunctionCalculationSetup> partitionFunctionCalculationSetupFile;
    try
    {
        // Parse setup file
        partitionFunctionCalculationSetupFile =
            MCSTHARPartitionFunctionCalculationSetup_(i_partitionFunctionCalculationSetupFileName,
                                                      xml_schema::flags::dont_validate);
    }
    catch(const xml_schema::exception& exception)
    {
        // Error during setup file parsing, throw exception
        string errorMessage("Unable to parse the specified setup file. ");
        errorMessage += exception.what();
        throw HadronizationException(errorMessage,__FUNCTION__,411);
    }

    // Check if provided string is empty
    if(partitionFunctionCalculationSetupFile->hadronDataFileName().empty())
    {
        // Empty data file name provided, throw exception
        throw HadronizationException("Error during parition function calculation setup loading, empty hadron data file name provided",
                                     __FUNCTION__,412);
    }
    m_partitionFunctionCalculationSetup.hadronDataSetFileName = partitionFunctionCalculationSetupFile->hadronDataFileName();

    
    // Set charge configuration
    m_partitionFunctionCalculationSetup.chargeConfiguration.strangeCharge =
        partitionFunctionCalculationSetupFile->chargeConfiguration().strangeCharge();
    m_partitionFunctionCalculationSetup.chargeConfiguration.charmCharge =
        partitionFunctionCalculationSetupFile->chargeConfiguration().charmCharge();
    m_partitionFunctionCalculationSetup.chargeConfiguration.bottomCharge =
        partitionFunctionCalculationSetupFile->chargeConfiguration().bottomCharge();
    m_partitionFunctionCalculationSetup.chargeConfiguration.electricCharge =
        static_cast<double>(partitionFunctionCalculationSetupFile->chargeConfiguration().electricCharge());
    m_partitionFunctionCalculationSetup.chargeConfiguration.baryonicCharge =
        static_cast<double>(partitionFunctionCalculationSetupFile->chargeConfiguration().baryonicCharge());
        
    // Set microcanonical parameter grid
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.minEnergyDensity =
        partitionFunctionCalculationSetupFile->microcanonicalParameterGridStructure().minEnergyDensity();
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.maxEnergyDensity =
        partitionFunctionCalculationSetupFile->microcanonicalParameterGridStructure().maxEnergyDensity();
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.numberOfEnergyDensityValues =
        partitionFunctionCalculationSetupFile->microcanonicalParameterGridStructure().numberOfEnergyDensityValues();
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.minGammaS =
        partitionFunctionCalculationSetupFile->microcanonicalParameterGridStructure().minStrangenessSuppressionParameter();
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.maxGammaS =
        partitionFunctionCalculationSetupFile->microcanonicalParameterGridStructure().maxStrangenessSuppressionParameter();
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.numberOfGammaSValues =
        partitionFunctionCalculationSetupFile->microcanonicalParameterGridStructure().numberOfStrangenessSuppressionParameterValues();
    m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.totalNumberOfGridPoints =
        m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.numberOfEnergyDensityValues*
        m_partitionFunctionCalculationSetup.microcanonicalParameterGrid.numberOfGammaSValues;
    
    // Set mass value list
    // TODO: add cross-check for sequence lenght (also in other loaders)
    const unsigned int numberOfMassValues(partitionFunctionCalculationSetupFile->massValueList().length());
    m_partitionFunctionCalculationSetup.massValueList.resize(numberOfMassValues);
    for(unsigned int massValueIndex=0;massValueIndex<numberOfMassValues;++massValueIndex)
    {
        m_partitionFunctionCalculationSetup.massValueList[massValueIndex] =
            partitionFunctionCalculationSetupFile->massValueList().massValue().at(massValueIndex);
    }

    // Set integrationd data
    m_partitionFunctionCalculationSetup.maxChannelSamplingNumber =
        partitionFunctionCalculationSetupFile->integrationParameters().maxChannelSamplingNumber();
    m_partitionFunctionCalculationSetup.integrationErrorThreshold =
        partitionFunctionCalculationSetupFile->integrationParameters().integrationErrorThreshold();
        

    // Retrieve random number generator seed
    if(partitionFunctionCalculationSetupFile->randomNumberGeneratorSeed().present())
    {
        m_partitionFunctionCalculationSetup.randomNumberGeneratorSeed =
            partitionFunctionCalculationSetupFile->randomNumberGeneratorSeed().get();
    }

    // Set output folder
    m_partitionFunctionCalculationSetup.outputFolder = partitionFunctionCalculationSetupFile->outputFolder();
    
    // Update setup loading status
    m_partitionFunctionCalculationSetupAvailability = true;
    
    return;
}

PartitionFunctionCalculationSetupLoader::~PartitionFunctionCalculationSetupLoader(void)
{
}

const PartitionFunctionCalculationSetup& PartitionFunctionCalculationSetupLoader::getPartitionFunctionCalculationSetup(void) const
{
    if(!m_partitionFunctionCalculationSetupAvailability)
    {
        // Partition function calculation parameters not correctly set, throw exception
        throw HadronizationException("Error during partition function calculation parameter retrieveal, parameters not correctly set",
                                     __FUNCTION__,413);
    }
    else
    {
        return m_partitionFunctionCalculationSetup;
    }
}
