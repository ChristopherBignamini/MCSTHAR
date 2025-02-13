#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/fileUtils.h"
#include <sstream>
#include <fstream>

#include "../../HadronizationPartitionFunctionStorage/Include/PartitionFunctionArchiveWriter.h"


// TODO: add support for win?
PartitionFunctionArchiveWriter::PartitionFunctionArchiveWriter(const string& i_partitionFunctionArchiveFolder,
                                                               const MicrocanonicalParameterGridStructure& i_microcanonicalParameterGridStructure,
                                                               const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                               const vector<double>& i_massValues)
                                                              :m_microcanonicalParameterGridStructure(i_microcanonicalParameterGridStructure)
                                                              ,m_massValues(i_massValues)// TODO: add check on mass value set
                                                              ,m_partitionFunctionArchive(NULL)
                                                              ,m_partitionFunctionSummaryFile(NULL)
                                                              ,m_partitionFunctionDataFile(NULL)
                                                              ,m_referenceMassIndex(0)
                                                              ,m_isNewMassData(true)
{
    try
    {
        // TODO: to be factorized for simulation step partition function data handling
        // Check existence of PartitionFunctionArchive file in provided folder
        const vector<string> archiveFolderContent(getDirectoryContent(i_partitionFunctionArchiveFolder));
        const vector<string>::const_iterator archiveFileIt(find(archiveFolderContent.begin(),
                                                                archiveFolderContent.end(),
                                                                s_partitionFunctionArchiveFileName));
        
        const string partitionFunctionArchiveFullFileName(i_partitionFunctionArchiveFolder +
                                                          "/" +
                                                          s_partitionFunctionArchiveFileName);
        if(archiveFileIt!=archiveFolderContent.end())
        {
            unique_ptr<PartitionFunctionArchiveFile> loadedPartitionFunctionArchive;
            try
            {
                // Parse existing PartitionFunctionArchive file
                loadedPartitionFunctionArchive = PartitionFunctionArchiveFile_(partitionFunctionArchiveFullFileName,
                                                                               xml_schema::flags::dont_validate);
            }
            catch(const xml_schema::exception& xmlException)
            {
                // Error during setup file parsing, throw exception
                string errorMessage("Unable to parse the existing partition function archive file. ");
                errorMessage += xmlException.what();
                // TODO: exit code everywhere
                throw HadronizationException(errorMessage,__FUNCTION__,905);
            }
            
            // Release auto pointer and pass ownership to data member
            m_partitionFunctionArchive = loadedPartitionFunctionArchive.release();
            
            // TODO: make here the required checks (for example the microcanonical parameter grid structure)
            
        }
        else
        {
            try
            {
                // PartitionFunctionArchive does not exist, create a new one
                m_partitionFunctionArchive = new PartitionFunctionArchiveFile(i_microcanonicalParameterGridStructure,
                                                                              SingleChargeConfigurationArchiveList(0));
            }
            catch(const xml_schema::exception& xmlException)
            {
                // TODO: add codesynt err. and exit code
                // Retrieve error message
                string errorMessage("Error during partition function archive summary creation. ");
                errorMessage += xmlException.what();
                throw HadronizationException(errorMessage,__FUNCTION__,901);
            }
        }
     
        
        // Create current charge configuration
        // Update current charge configuration
        // TODO: this assignement can be avoided if m_currentChargeConfiguration is already i_chargeConfiguration
        // TODO: can't be this cast avoided?
        const ChargeConfiguration currentChargeConfiguration(static_cast<int>(i_chargeConfiguration.electricCharge),
                                                             static_cast<int>(i_chargeConfiguration.baryonicCharge),
                                                             i_chargeConfiguration.strangeCharge,
                                                             i_chargeConfiguration.charmCharge,
                                                             i_chargeConfiguration.bottomCharge);

        // Build current charge partition data folder full name
        const string currentChargeFolderName(createNewChargeConfigurationPartitionFunctionFolderName(currentChargeConfiguration));
        m_currentChargeFullFolderName = i_partitionFunctionArchiveFolder + "/" + currentChargeFolderName;
        
        // TODO: Add check presence of current charge configuration partition function data in archive file
        // TODO: implement data merging procedure
        
        // Current charge configuration not present, create partition function data folder
        createDirectory(m_currentChargeFullFolderName);

        // Create PartitionFunctionDataSummary
        // TODO: move into a function
        m_partitionFunctionSummaryFile = new PartitionFunctionDataSummaryFile(currentChargeConfiguration,
                                                                              MassGridDataElementList(0));
        m_partitionFunctionSummaryFullFileName = m_currentChargeFullFolderName + "/" + s_partitionFunctionDataSummaryFileName;

        // Update PartitionFunctionArchive
        updatePartitionFunctionArchiveFile(partitionFunctionArchiveFullFileName,
                                           currentChargeConfiguration,
                                           currentChargeFolderName);

        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

PartitionFunctionArchiveWriter::~PartitionFunctionArchiveWriter(void)
{
    if(m_partitionFunctionArchive != NULL)
    {
        delete m_partitionFunctionArchive;
    }
    if(m_partitionFunctionSummaryFile != NULL)
    {
        delete m_partitionFunctionSummaryFile;
    }
    if(m_partitionFunctionDataFile != NULL)
    {
        delete m_partitionFunctionDataFile;    
    }
}

// TODO: refactor with existing write data method
// TODO: specificare che questo metodo deve essere usato in modo "sequenziale" per le densità e
// che nessun controllo viene fatto su superamento del numero id valori di densità
void PartitionFunctionArchiveWriter::updatePartitionFunctionData(const SingleMassEnergyDensityPairPartitionFunctionData& i_partitionFunctionData,
                                                                 const unsigned int i_massIndex)
{
    // Check input index validity
    if(i_massIndex>=m_massValues.size())
    {
        throw HadronizationException("Error during partition function data writing: not identified mass value for the provided index",
                                     __FUNCTION__,903);
    }
    
    // Check input size with respect to stred microcanonical grid size // TODO: to be update with SingleMassEnergyDensityPairPartitionFunctionData
                                                                       // switch into a class
    if((i_partitionFunctionData.partitionFunctions.size() !=
        m_microcanonicalParameterGridStructure.numberOfStrangenessSuppressionParameterValues()) ||
       (i_partitionFunctionData.partitionFunctionErrors.size() !=
        m_microcanonicalParameterGridStructure.numberOfStrangenessSuppressionParameterValues()) ||
       (i_partitionFunctionData.partitionFunctionErrorStatus.size() !=
        m_microcanonicalParameterGridStructure.numberOfStrangenessSuppressionParameterValues()) ||
       (i_partitionFunctionData.numberOfChannelSamplings.size() !=
        m_microcanonicalParameterGridStructure.numberOfStrangenessSuppressionParameterValues()))
    {
        throw HadronizationException("Error during partition function data writing: provided data size not equal to strangeness suppression parameter grid size",
                                     __FUNCTION__,909);
    }
    
    // TODO: refine logic and use file presence check instead
    // Check if required data writing is for a new mass value
    if(m_referenceMassIndex!=i_massIndex)
    {
        m_isNewMassData = true;
    }
    if(m_isNewMassData)
    {
        // Create new partition function data file and update summary
        try
        {
            // Create empty partition function data file
            if(m_partitionFunctionDataFile!=NULL)
            {
                delete m_partitionFunctionDataFile;
            }
            const double& massValue(m_massValues[i_massIndex]);
            m_partitionFunctionDataFile = new PartitionFunctionDataFile(massValue,
                                                                        PartitionFunctionDataList(0));
            // Create filenames
            // TODO: add win support?
            stringstream massGridPointIndexString;
            massGridPointIndexString<<i_massIndex;
            const string partitionFunctionDataFileName(s_partitionFunctionDataFileNamePrefix +
                                                       massGridPointIndexString.str() +
                                                       ".xml");

            
            m_partitionFunctionDataFullFileName = m_currentChargeFullFolderName +
                                                  "/" +
                                                  partitionFunctionDataFileName;
            
            // Update data summary file // TODO: if error are found in next
                                                    // writing steps the summary file should be
                                                    // accordingly corrected
            updatePartitionFunctionDataSummaryFile(massValue,
                                                   partitionFunctionDataFileName);
            
            // Update storage status
            m_referenceMassIndex=i_massIndex;
            m_isNewMassData = false;
            
        }
        catch(const xml_schema::exception& xmlException)
        {
            // TODO: add codesynt err. and exit code
            // Retrieve error message
            string errorMessage("Error during partition function data file creation. ");
            errorMessage += xmlException.what();
            throw HadronizationException(errorMessage,__FUNCTION__,908);
        }
    }
    
    // Update partition function data file
    try
    {
        updatePartitionFunctionDataFile(i_partitionFunctionData);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
    return;
}

void PartitionFunctionArchiveWriter::updatePartitionFunctionArchiveFile(const string& i_partitionFunctionArchiveFullFileName,
                                                                        const ChargeConfiguration& i_chargeConfiguration,
                                                                        const string& i_chargeFolderName)
{
    // Update charge configuration list if needed
    // TODO: implement existence check for current charge configuration
    ++m_partitionFunctionArchive->singleChargeConfigurationArchiveList().length();
    m_partitionFunctionArchive->singleChargeConfigurationArchiveList().singleChargeConfigurationArchive().
        push_back(SingleChargeConfigurationArchive(i_chargeConfiguration,
                                                   i_chargeFolderName));
 
    // Update PartitionFunctionArchiveFile on disk
    // TODO: check for best practice
    ofstream xmlStream(i_partitionFunctionArchiveFullFileName.c_str());
    try
    {
        PartitionFunctionArchiveFile_(xmlStream,
                                      *m_partitionFunctionArchive);
    }
    catch(const xml_schema::exception& xmlException)
    {
        // Retrieve error message
        string errorMessage("Error during partition function archive file writing. ");
        errorMessage += xmlException.what();
        throw HadronizationException(errorMessage,__FUNCTION__,904);
    }

    return;
}

void PartitionFunctionArchiveWriter::updatePartitionFunctionDataSummaryFile(const double i_massValue,
                                                                            const string& i_partitionFunctionDataFileName)
{
    // TODO: add input check
    ++m_partitionFunctionSummaryFile->massGridDataElementList().length();
    m_partitionFunctionSummaryFile->massGridDataElementList().massGridDataElement().
        push_back(MassGridDataElement(i_massValue,
                                      i_partitionFunctionDataFileName));
    
    // Update PartitionFunctionArchiveFile on disk
    // TODO: check for best practice
    ofstream xmlStream(m_partitionFunctionSummaryFullFileName.c_str());
    try
    {
        PartitionFunctionDataSummaryFile_(xmlStream,
                                          *m_partitionFunctionSummaryFile);
    }
    catch(const xml_schema::exception& xmlException)
    {
        // Retrieve error message
        string errorMessage("Error during partition function summary file writing. ");
        errorMessage += xmlException.what();
        throw HadronizationException(errorMessage,__FUNCTION__,904);
    }
    
    return;
}

void PartitionFunctionArchiveWriter::updatePartitionFunctionDataFile(const SingleMassEnergyDensityPairPartitionFunctionData& i_partitionFunctionData)
{

    // Update partition function data
    unsigned int dataGlobalIndex(m_partitionFunctionDataFile->partitionFunctionDataList().length());
    // Update data list lenght
    m_partitionFunctionDataFile->partitionFunctionDataList().length() +=
        m_microcanonicalParameterGridStructure.numberOfStrangenessSuppressionParameterValues();
    // Insert new data elements (For first element in data list the parameters of the sampling function used for calculation are also provided)
    PartitionFunctionGridElement firstElementData(dataGlobalIndex,
                                                  i_partitionFunctionData.partitionFunctions.at(0),
                                                  i_partitionFunctionData.partitionFunctionErrors.at(0),
                                                  i_partitionFunctionData.numberOfChannelSamplings.at(0));
    firstElementData.importanceSamplingFunctionParameters(
                                ImportanceSamplingFunctionParameters(i_partitionFunctionData.grandCanonicalParameters.electricChargeFugacity,
                                                                     i_partitionFunctionData.grandCanonicalParameters.baryonicChargeFugacity,
                                                                     i_partitionFunctionData.grandCanonicalParameters.strangeChargeFugacity,
                                                                     i_partitionFunctionData.grandCanonicalParameters.charmChargeFugacity,
                                                                     i_partitionFunctionData.grandCanonicalParameters.bottomChargeFugacity,
                                                                     i_partitionFunctionData.grandCanonicalParameters.temperature));
        
    m_partitionFunctionDataFile->partitionFunctionDataList().partitionFunctionData().push_back(firstElementData);
    ++dataGlobalIndex;
    for(unsigned int dataIndex=1;
        dataIndex<m_microcanonicalParameterGridStructure.numberOfStrangenessSuppressionParameterValues();
        ++dataIndex,++dataGlobalIndex)
    {
        m_partitionFunctionDataFile->partitionFunctionDataList().partitionFunctionData().
            push_back(PartitionFunctionGridElement(dataGlobalIndex,
                                                   i_partitionFunctionData.partitionFunctions.at(dataIndex),
                                                   i_partitionFunctionData.partitionFunctionErrors.at(dataIndex),
                                                   i_partitionFunctionData.numberOfChannelSamplings.at(dataIndex)));
    }
    
    // Open output file stream
    ofstream xmlStream(m_partitionFunctionDataFullFileName.c_str());
    
    // Update data file
    try
    {
        PartitionFunctionDataFile_(xmlStream,
                                   *m_partitionFunctionDataFile);
    }
    catch(const xml_schema::exception& xmlException)
    {
        // Retrieve error message
        string errorMessage("Error during partition function data file updating. ");
        errorMessage += xmlException.what();
        throw HadronizationException(errorMessage,__FUNCTION__,902);
    }
    
    return;
}


string PartitionFunctionArchiveWriter::createNewChargeConfigurationPartitionFunctionFolderName(const ChargeConfiguration& i_chargeConfiguration) const
{
    string o_partitionFunctionFolderName;
    
    // Set strange charge
    o_partitionFunctionFolderName += createChargeValueString(i_chargeConfiguration.strangeCharge());
    o_partitionFunctionFolderName += "_";
    
    // Set electric charge
    o_partitionFunctionFolderName += createChargeValueString(i_chargeConfiguration.electricCharge());
    o_partitionFunctionFolderName += "_";
    
    // Set baryonic charge
    o_partitionFunctionFolderName += createChargeValueString(i_chargeConfiguration.baryonicCharge());
    o_partitionFunctionFolderName += "_";
    
    // set charm charge
    o_partitionFunctionFolderName += createChargeValueString(i_chargeConfiguration.charmCharge());
    o_partitionFunctionFolderName += "_";
    
    // Set bottom charge
    o_partitionFunctionFolderName += createChargeValueString(i_chargeConfiguration.bottomCharge());
    
    return o_partitionFunctionFolderName;
}


string PartitionFunctionArchiveWriter::createChargeValueString(const int i_chargeValue) const
{
    string o_chargeValueString;
    
    // Set charge sign char
    if(i_chargeValue<0)
    {
        o_chargeValueString += "1";
    }
    else
    {
        o_chargeValueString += "0";
    }

    // Set charge absolute value chars
    stringstream chargeAbsoluteValue;
    chargeAbsoluteValue<<abs(i_chargeValue);
    o_chargeValueString += chargeAbsoluteValue.str();
    
    return o_chargeValueString;
}
