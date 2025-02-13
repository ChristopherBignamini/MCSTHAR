#include "../Include/PartitionFunctionHandling.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../HadronizationPartitionFunctionStorage/Include/partitionFunctionDataConstants.h"
#include "../../../Utilities/Include/fileUtils.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
// TODO: build something like part func data IO element and do the detailed factorization!
#include "../../HadronizationPartitionFunctionStorage/Include/PartitionFunctionDataSummaryFile.h"
#include "../../HadronizationPartitionFunctionStorage/Include/PartitionFunctionDataFile.h"
#include "../../HadronizationPartitionFunctionStorage/Include/PartitionFunctionArchiveWriter.h"

PartitionFunctionHandling::PartitionFunctionHandling(const string& i_partitionFunctionDataFolder,
                                                     const double i_energyDensity,
                                                     const double i_gammaS,
                                                     const bool isNewMethod)
                                                    :m_partitionFunctionDataFolder(i_partitionFunctionDataFolder)
                                                    ,m_energyDensity(i_energyDensity)
                                                    ,m_gammaS(i_gammaS)
                                                    ,m_strangenessSuppressionParameterPointsNumber(0)
                                                    ,m_grid2DDimension(0)
                                                    ,m_index00(0)
                                                    ,m_index10(0)
                                                    ,m_index01(0)
                                                    ,m_index11(0)
                                                    ,m_weight00(0.)
                                                    ,m_weight10(0.)
                                                    ,m_weight01(0.)
                                                    ,m_weight11(0.)
                                                    ,m_isInterpolationDataAvailable(false)
                                                    ,m_is1DGridDataAvailable(false)
{
    try
    {
        if(m_energyDensity<=0.)
        {
            throw HadronizationException("Error during partition function handling object creation, non positive energy density parameter provided",__FUNCTION__,131);
        }
        if(m_gammaS<0. || m_gammaS>1.)
        {
            throw HadronizationException("Error during partition function handling object creation, strangeness suppression parameter outside [0,1] range provided",__FUNCTION__,131);
        }
        
        // Load partition function archive data (e.g., available charge configuration and microcanonical parameter grid structure)
        loadPartitionFunctionArchiveStructureData();
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void PartitionFunctionHandling::loadPartitionFunctionArchiveStructureData(void)
{
    // TODO: to be factorized for calculation step partition function data writing
    // Check existence of PartitionFunctionArchive file
    const vector<string> archiveFolderContent(getDirectoryContent(m_partitionFunctionDataFolder));
    const vector<string>::const_iterator archiveFileIt(find(archiveFolderContent.begin(),
                                                            archiveFolderContent.end(),
                                                            s_partitionFunctionArchiveFileName));
    
    if(archiveFileIt!=archiveFolderContent.end())
    {
        const string partitionFunctionArchiveFullFileName(m_partitionFunctionDataFolder +
                                                          "/" +
                                                          s_partitionFunctionArchiveFileName);
        unique_ptr<PartitionFunctionArchiveFile> partitionFunctionArchiveFile;
        try
        {
            // Parse existing PartitionFunctionArchive file
            // TODO: win support?
            partitionFunctionArchiveFile = PartitionFunctionArchiveFile_(partitionFunctionArchiveFullFileName,
                                                                         xml_schema::flags::dont_validate);
        }
        catch(const xml_schema::exception& xmlException)
        {
            // Error during setup file parsing, throw exception
            string errorMessage("Error during parsing of partition function archive file " +
                                partitionFunctionArchiveFullFileName +
                                ". ");
            errorMessage += xmlException.what();
            // TODO: exit code everywhere
            throw HadronizationException(errorMessage,__FUNCTION__,132);
        }
        
        // TODO: make here the required checks (for example the microcanonical parameter grid structure)
        
        // Set single charge configuration data set folder names
        setPartitionFunctionDataSetFileList(partitionFunctionArchiveFile);
        
        // Set available microcanonical parameter grid structure
        setMicrocanonicalParameterGridStructure(partitionFunctionArchiveFile);
        
        // Compute interpolation data on microcanonical parameter grid
        computeMicrocanonicalParameterGridInterpolationData();
                
    }
    else
    {
        throw HadronizationException("Error during partition function archive structure retrieval, missing archive summary file",
                                     __FUNCTION__,133);
    }
    
    return;
}

void PartitionFunctionHandling::setPartitionFunctionDataSetFileList(const unique_ptr<PartitionFunctionArchiveFile>& i_partitionFunctionArchiveFile)
{
    // Check number of single charge configuration partition function data set
    if(i_partitionFunctionArchiveFile->singleChargeConfigurationArchiveList().length()==0)
    {
        throw HadronizationException("Error during partition function data set structure retrieval, no single charge configuration partition function data set found",
                                     __FUNCTION__,134);
    }
    
    // Loop over available single charge configuration partition function data sets
    for(unsigned int folderIndex=0;folderIndex<i_partitionFunctionArchiveFile->singleChargeConfigurationArchiveList().length();++folderIndex)
    {                
        // Store charge configuration and data set path
        const SingleChargeConfigurationArchive& singleChargeConfigurationData(
                        i_partitionFunctionArchiveFile->singleChargeConfigurationArchiveList().singleChargeConfigurationArchive().at(folderIndex));
        
        // TODO: Extend to Windows case
        // TODO: avoid static cast
        // TODO: avoid duplicated charge configuration structure
        m_partitionFunctionDataSetFolders[MCSTHAR::Utilities::ChargeConfiguration(
                                                        singleChargeConfigurationData.chargeConfiguration().strangeCharge(),
                                                        singleChargeConfigurationData.chargeConfiguration().charmCharge(),
                                                        singleChargeConfigurationData.chargeConfiguration().bottomCharge(),
                                                        static_cast<double>(singleChargeConfigurationData.chargeConfiguration().electricCharge()),
                                                        static_cast<double>(singleChargeConfigurationData.chargeConfiguration().baryonicCharge()))]
                                                        = m_partitionFunctionDataFolder +
                                                          "/" +
                                                          singleChargeConfigurationData.folderName();
    }
                                                                                  
    return;
}

// TODO: check for code refactoring with partition function writing class
// TODO: required, no check is actually performed on microcanonical grid structure!
void PartitionFunctionHandling::setMicrocanonicalParameterGridStructure(const unique_ptr<PartitionFunctionArchiveFile>& i_partitionFunctionArchiveFile)
{
    // TODO: warning, the same grid is assumed for all charge configurations
    // Set microcanonical parameter grid dimensions
    const unsigned int energyDensityPointsNumber(i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().numberOfEnergyDensityValues());
    if(energyDensityPointsNumber==0)
    {
        throw HadronizationException("Error in microcanonical parameter grid structure, zero energy density available values",
                                     __FUNCTION__,135);
    }
    m_strangenessSuppressionParameterPointsNumber = (i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().numberOfStrangenessSuppressionParameterValues());
    if(m_strangenessSuppressionParameterPointsNumber==0)
    {
        throw HadronizationException("Error in microcanonical parameter grid structure, zero strangeness suppression parameter available values",
                                     __FUNCTION__,135);
    }
    m_energyDensityValues.resize(energyDensityPointsNumber);
    m_gammaSValues.resize(m_strangenessSuppressionParameterPointsNumber);
    m_grid2DDimension = energyDensityPointsNumber*m_strangenessSuppressionParameterPointsNumber;
    
    // Set energy density grid values
    // Check parameter minimum value
    if(i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minEnergyDensity()<=0.)
    {
        throw HadronizationException("Error in microcanonical parameter grid structure, non positive minimum energy density value provided",
                                     __FUNCTION__,135);
    }
    m_energyDensityValues[0] = i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minEnergyDensity();
    if(energyDensityPointsNumber>1)
    {
        const double parameterStep((i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().maxEnergyDensity() -
                                   i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minEnergyDensity())/(energyDensityPointsNumber-1));
        
        // Check parameter min/max value ordering
        if(parameterStep<=0.)
        {
            throw HadronizationException("Error in microcanonical parameter grid structure, wrong ordering in energy density min/max values",
                                         __FUNCTION__,135);
        }
        for(unsigned int gridPointIndex=1;gridPointIndex<energyDensityPointsNumber;++gridPointIndex)
        {
            m_energyDensityValues[gridPointIndex] = i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minEnergyDensity() +
                                                    gridPointIndex*parameterStep;
        }
    }

    // Set strangeness suppression parameter grid values
    // Check parameter minimum value
    if(i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minStrangenessSuppressionParameter()<0.)
    {
        throw HadronizationException("Error in microcanonical parameter grid structure, negative minimum strangeness suppression parameter value provided",
                                     __FUNCTION__,135);
    }
    m_gammaSValues[0] = i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minStrangenessSuppressionParameter();
    if(m_strangenessSuppressionParameterPointsNumber>1)
    {
    
        // Check max value
        if(i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().maxStrangenessSuppressionParameter()>1.)
        {
            throw HadronizationException("Error in microcanonical parameter grid structure, strangeness suppression parameter max values larger than 1 provided",
                                         __FUNCTION__,135);
        }

        const double parameterStep((i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().maxStrangenessSuppressionParameter() -
                                    i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minStrangenessSuppressionParameter())/(m_strangenessSuppressionParameterPointsNumber-1));
        
        // Check parameter min/max value ordering
        if(parameterStep<=0.)
        {
            throw HadronizationException("Error in microcanonical parameter grid structure, wrong ordering in strangeness suppression parameter min/max values",
                                         __FUNCTION__,135);
        }
        
        for(unsigned int gridPointIndex=1;gridPointIndex<m_strangenessSuppressionParameterPointsNumber;++gridPointIndex)
        {
            m_gammaSValues[gridPointIndex] =
                i_partitionFunctionArchiveFile->microcanonicalParameterGridStructure().minStrangenessSuppressionParameter() +
                gridPointIndex*parameterStep;
        }
    }

    return;
}

void PartitionFunctionHandling::computeMicrocanonicalParameterGridInterpolationData(void)
{

    bool interpolationErrorStatus(false);
    const SingleParameterInterpolationData energyInterpDensityData(computeInterpolationData(m_energyDensityValues,
                                                                                            m_energyDensity,
                                                                                            interpolationErrorStatus));
    // Check interpolation error status
    if(interpolationErrorStatus)
    {
        throw HadronizationException("Error during interpolation data calculation (energy density), interpolation value outside the available range",
                                     __FUNCTION__,136);
    }

    const SingleParameterInterpolationData strangenessSuppressionInterpData(computeInterpolationData(m_gammaSValues,
                                                                                                     m_gammaS,
                                                                                                     interpolationErrorStatus));
    // Check interpolation error status
    if(interpolationErrorStatus)
    {
        throw HadronizationException("Error during interpolation data calculation (gamma_S), interpolation value outside the available range",
                                     __FUNCTION__,136);
    }

    // Compute 2D grid interpolation indexes
    m_index00 = strangenessSuppressionInterpData.gridPoints.first +
    m_strangenessSuppressionParameterPointsNumber*energyInterpDensityData.gridPoints.first;
    m_index10 = m_index00+1;
    m_index01 = m_index00+m_strangenessSuppressionParameterPointsNumber;
    m_index11 = m_index01+1;

    // Compute 2D grid interpolation weights
    m_weight00 = strangenessSuppressionInterpData.weights.first*energyInterpDensityData.weights.first;
    m_weight10 = strangenessSuppressionInterpData.weights.second*energyInterpDensityData.weights.first;
    m_weight01 = strangenessSuppressionInterpData.weights.first*energyInterpDensityData.weights.second;
    m_weight11 = strangenessSuppressionInterpData.weights.second*energyInterpDensityData.weights.second;
    
    m_isInterpolationDataAvailable=true;

    return;
}

// TODO: cluster mass value seems to differ wrt reference one at some precision level, causing a different value in mass interp. coeff.
double PartitionFunctionHandling::getPartitionFunctionValue(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                            const double i_mass,
                                                            unsigned int& io_calculationErrorStatus,
                                                            const bool isNewMethod)
{
    try
    {
        MCSTHAR::Utilities::ChargeConfiguration chargeConfiguration(i_chargeConfiguration);

        
        // Look for provided charge configuration
        map<MCSTHAR::Utilities::ChargeConfiguration,PartitionFunctionMassGridData>::const_iterator massGridIt =
        m_partitionFunctionMassGrids.find(chargeConfiguration);
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Look for anti configuration
            chargeConfiguration = chargeConfiguration.getAntiConfiguration();
            massGridIt = m_partitionFunctionMassGrids.find(chargeConfiguration);
        }
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Data currently not available for the provided charge configuration.
            // Look for required data set existence
            map<MCSTHAR::Utilities::ChargeConfiguration,string>::iterator dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            
            if(dataSetIt==m_partitionFunctionDataSetFolders.end())
            {
                // Look for original configuration
                chargeConfiguration = chargeConfiguration.getAntiConfiguration();
                dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            }
            
            if(dataSetIt!=m_partitionFunctionDataSetFolders.end())
            {
                // Load new data set and build partition function subset
                updatePartitionFunctionSubSet(dataSetIt->first,
                                              dataSetIt->second,
                                              isNewMethod);
            }
            else
            {
                // Required data not found in provided partition function data set
                io_calculationErrorStatus = 1;
                return 0.;
            }
        }
        
        // Perform 1D linear interpolation
        PartitionFunctionMassGridData* massGridData = &m_partitionFunctionMassGrids[chargeConfiguration];
        bool interpolationError(false);
        const SingleParameterInterpolationData massInterpData(computeInterpolationData(massGridData->mass,
                                                                                       i_mass,
                                                                                       interpolationError));
                
        // Check interpolation error status
        if(interpolationError)
        {
            // Interpolation mass value outside available range
            io_calculationErrorStatus = 2;
            return 0.;
        }
        
        // Interpolation parameter calculation correctly performed
        io_calculationErrorStatus = 0;
        
        
//        cout<<"massInterpData.gridPoints.first "<<massInterpData.gridPoints.first
//            <<" massInterpData.gridPoints.second "<<massInterpData.gridPoints.second
//            <<" massGridData->partitionFunction[massInterpData.gridPoints.first] "<<massGridData->partitionFunction[massInterpData.gridPoints.first]
//            <<" massGridData->partitionFunction[massInterpData.gridPoints.second] "<<massGridData->partitionFunction[massInterpData.gridPoints.second]
//            <<" massInterpData.weights.first "<<massInterpData.weights.first
//            <<" massInterpData.weights.second "<<massInterpData.weights.second
//            <<" massGridData->mass[massInterpData.gridPoints.first] "<<massGridData->mass[massInterpData.gridPoints.first]
//            <<" massGridData->mass[massInterpData.gridPoints.second] "<<massGridData->mass[massInterpData.gridPoints.second]<<endl;
        
        return massGridData->partitionFunction[massInterpData.gridPoints.first]*massInterpData.weights.first +
        massGridData->partitionFunction[massInterpData.gridPoints.second]*massInterpData.weights.second;
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

double PartitionFunctionHandling::getPartitionFunctionValueErr(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                            const double i_mass,
                                                            unsigned int& io_calculationErrorStatus,
                                                            const bool isNewMethod)
{
    try
    {
        MCSTHAR::Utilities::ChargeConfiguration chargeConfiguration(i_chargeConfiguration);
        
        
        // Look for provided charge configuration
        map<MCSTHAR::Utilities::ChargeConfiguration,PartitionFunctionMassGridData>::const_iterator massGridIt =
        m_partitionFunctionMassGrids.find(chargeConfiguration);
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Look for anti configuration
            chargeConfiguration = chargeConfiguration.getAntiConfiguration();
            massGridIt = m_partitionFunctionMassGrids.find(chargeConfiguration);
        }
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Data currently not available for the provided charge configuration.
            // Look for required data set existence
            map<MCSTHAR::Utilities::ChargeConfiguration,string>::iterator dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            
            if(dataSetIt==m_partitionFunctionDataSetFolders.end())
            {
                // Look for original configuration
                chargeConfiguration = chargeConfiguration.getAntiConfiguration();
                dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            }
            
            if(dataSetIt!=m_partitionFunctionDataSetFolders.end())
            {
                // Load new data set and build partition function subset
                updatePartitionFunctionSubSet(dataSetIt->first,
                                              dataSetIt->second,
                                              isNewMethod);
            }
            else
            {
                // Required data not found in provided partition function data set
                io_calculationErrorStatus = 1;
                return 0.;
            }
        }
        
        // Perform 1D linear interpolation
        PartitionFunctionMassGridData* massGridData = &m_partitionFunctionMassGrids[chargeConfiguration];
        bool interpolationError(false);
        const SingleParameterInterpolationData massInterpData(computeInterpolationData(massGridData->mass,
                                                                                       i_mass,
                                                                                       interpolationError));
        
        // Check interpolation error status
        if(interpolationError)
        {
            // Interpolation mass value outside available range
            io_calculationErrorStatus = 2;
            return 0.;
        }
        
        // Interpolation parameter calculation correctly performed
        io_calculationErrorStatus = 0;
        
        return massGridData->partitionFunctionErr[massInterpData.gridPoints.first]*massInterpData.weights.first +
        massGridData->partitionFunctionErr[massInterpData.gridPoints.second]*massInterpData.weights.second;
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}


// TODO: check all new function for optimization, useless code, etc..
void PartitionFunctionHandling::updatePartitionFunctionSubSet(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                              const string& i_dataFolder,
                                                              const double isNewMethod)
{
    // TODO: check refactoring with writing class
    
    // Load partition function data summary for current charge configuration
    // TODO: add support for win?
    const string filePath(i_dataFolder + "/");
    const string partitionFunctionDataSummaryFullFileName(filePath +
                                                          s_partitionFunctionDataSummaryFileName);
    unique_ptr<PartitionFunctionDataSummaryFile> partitionFunctionDataSummaryFile;
    try
    {
        partitionFunctionDataSummaryFile = PartitionFunctionDataSummaryFile_(partitionFunctionDataSummaryFullFileName,
                                                                             xml_schema::flags::dont_validate);
    }
    catch(const xml_schema::exception& xmlException)
    {
        // Error during setup file parsing, throw exception
        string errorMessage("Error during parsing the partition function data summary file " +
                            partitionFunctionDataSummaryFullFileName +
                            ". ");
        errorMessage += xmlException.what();
        // TODO: exit code everywhere
        throw HadronizationException(errorMessage,__FUNCTION__,137);
    }
    
    // TODO: add charge configuration check
    
    // Check number of single mass partition function data set
    if(partitionFunctionDataSummaryFile->massGridDataElementList().length()==0)
    {
        throw HadronizationException("Error during single charge configuration partition function data set structure retrieval, no single mass partition function data set found",
                                     __FUNCTION__,137);
    }
    
    // Loop over single mass partition function data files
    try
    {
        vector<double> partitionFunctionGrid;// TODO: serve ancora?
        vector<double> errpartitionFunctionGrid;// TODO: serve ancora?
        PartitionFunctionMassGridData& singleChargeConfigurationGrid(m_partitionFunctionMassGrids[i_chargeConfiguration]);
        for(unsigned int fileIndex=0;fileIndex<partitionFunctionDataSummaryFile->massGridDataElementList().length();++fileIndex)
        {
            
            // TODO: add check on mass value (positive, oredered, etc...)
            // TODO: avoid push back
            // Update mass value grid
            const MassGridDataElement& massGridDataElement(
                                                    partitionFunctionDataSummaryFile->massGridDataElementList().massGridDataElement().at(fileIndex));
            singleChargeConfigurationGrid.mass.push_back(massGridDataElement.massValue());
            
            // Load partition function data file
            loadPartitionFunctionDataFile(string(filePath +
                                                 massGridDataElement.dataFileName()),
                                          i_chargeConfiguration,
                                          partitionFunctionGrid,
                                          errpartitionFunctionGrid,
                                          isNewMethod);
        }
        
        // Update mass grid data for the considered charge configuration
        singleChargeConfigurationGrid.massValueNumber = partitionFunctionDataSummaryFile->massGridDataElementList().length();
        
        // Build partition function subset
        buildPartitionFunctionSubSet(i_chargeConfiguration,partitionFunctionGrid,errpartitionFunctionGrid);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }

    return;
}


void PartitionFunctionHandling::loadPartitionFunctionDataFile(const string& i_partitionFunctionDataFullFileName,
                                                              const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                              vector<double>& io_partitionFunctionGrid,
                                                              vector<double>& io_partitionFunctionErrGrid,
                                                              const bool isNewMethod)
{
    
    // Load partition function data 
    // TODO: add support for win?
    unique_ptr<PartitionFunctionDataFile> partitionFunctionDataFile;
    try
    {
        partitionFunctionDataFile = PartitionFunctionDataFile_(i_partitionFunctionDataFullFileName,
                                                               xml_schema::flags::dont_validate);
    }
    catch(const xml_schema::exception& xmlException)
    {
        // Error during setup file parsing, throw exception
        string errorMessage("Error during parsing the partition function data file " +
                            i_partitionFunctionDataFullFileName +
                            ". ");
        errorMessage += xmlException.what();
        // TODO: exit code everywhere
        throw HadronizationException(errorMessage,__FUNCTION__,138);
    }

    // Check number of partition function values versus microcanonical grid dimension
    if(partitionFunctionDataFile->partitionFunctionDataList().length()!=m_grid2DDimension)
    {
        const string errorMessage("Error during " +
                                  i_partitionFunctionDataFullFileName +
                                  " partition function data file, microcanonical parameter grid dimension and data multiplicity mismatch");
        throw HadronizationException(errorMessage,__FUNCTION__,139);
    }
    // Loop over single partition function data item
    for(unsigned int dataIndex=0;dataIndex<m_grid2DDimension;++dataIndex)
    {
        // TODO: avoid push_back, optimize everywhere
        io_partitionFunctionGrid.push_back(
                            partitionFunctionDataFile->partitionFunctionDataList().partitionFunctionData().at(dataIndex).partitionFunctionValue());
        io_partitionFunctionErrGrid.push_back(partitionFunctionDataFile->partitionFunctionDataList().partitionFunctionData().at(dataIndex).partitionFunctionErrorValue());
    }

    return;
}



































PartitionFunctionHandling::PartitionFunctionHandling(const string& i_partitionFunctionDataFolder,
                                                     const double i_energyDensity,
                                                     const double i_gammaS)
                                                    :m_partitionFunctionDataFolder(i_partitionFunctionDataFolder)
                                                    ,m_energyDensity(i_energyDensity)
                                                    ,m_gammaS(i_gammaS)
                                                    ,m_energyDensityPointNumber(0)
                                                    ,m_grid2DDimension(0)
                                                    ,m_index00(0)
                                                    ,m_index10(0)
                                                    ,m_index01(0)
                                                    ,m_index11(0)
                                                    ,m_weight00(0.)
                                                    ,m_weight10(0.)
                                                    ,m_weight01(0.)
                                                    ,m_weight11(0.)
                                                    ,m_isInterpolationDataAvailable(false)
                                                    ,m_is1DGridDataAvailable(false)
{
    try
    {
        if(m_energyDensity<=0.)
        {
            throw HadronizationException("Error during partition function handling object creation, non positive energy density parameter provided",__FUNCTION__,131);
        }
        if(m_gammaS<0. || m_gammaS>1.)
        {
            throw HadronizationException("Error during partition function handling object creation, strangeness suppression parameter outside [0,1] range provided",__FUNCTION__,131);
        }
        
        // Load partition function full set file list
        setPartitionFunctionDataSetFileList();
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

PartitionFunctionHandling::~PartitionFunctionHandling(void)
{
}

void PartitionFunctionHandling::setPartitionFunctionDataSetFileList(void)
{
    try
    {
        // Check main directory existence and retrieve subfolder list
        const vector<string> partitionFunctionFolderList(getDirectoryContent(m_partitionFunctionDataFolder));
        const unsigned int dataFolderNumber(partitionFunctionFolderList.size());
        if(dataFolderNumber==0)
        {
            throw HadronizationException("Error during partition function data set structure retrieval, no file found",
                                         __FUNCTION__,132);
        }
        
        // Loop over available folders (corresponding to single charge configuration partition function data)
        for(unsigned int folderIndex=0;folderIndex<dataFolderNumber;++folderIndex)
        {
            // Skip non data directories
            // TODO: refine this selection
            if(partitionFunctionFolderList[folderIndex].size()==10)
            {
                // Retrieve charge configuration from folder name
                const MCSTHAR::Utilities::ChargeConfiguration chargeConfiguration(
                                                                        retrieveChargeConfiguration(partitionFunctionFolderList[folderIndex]));
                
                // Store charge configuration and data set path
                // TODO: Extend to Windows case
                m_partitionFunctionDataSetFolders[chargeConfiguration] =
                    m_partitionFunctionDataFolder+"/"+partitionFunctionFolderList[folderIndex];
            }
            
        }
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

// TODO: cluster mass value seems to differ wrt reference one at some precision level, causing a different value in mass interp. coeff.
double PartitionFunctionHandling::getPartitionFunctionValue(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                            const double i_mass,
                                                            unsigned int& io_calculationErrorStatus)
{
    try
    {
        MCSTHAR::Utilities::ChargeConfiguration chargeConfiguration(i_chargeConfiguration);

        // Look for provided charge configuration
        map<MCSTHAR::Utilities::ChargeConfiguration,PartitionFunctionMassGridData>::const_iterator massGridIt =
            m_partitionFunctionMassGrids.find(chargeConfiguration);
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Look for anti configuration
            chargeConfiguration = chargeConfiguration.getAntiConfiguration();
            massGridIt = m_partitionFunctionMassGrids.find(chargeConfiguration);
        }
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Data currently not available for the provided charge configuration.
            // Look for required data set existence
            map<MCSTHAR::Utilities::ChargeConfiguration,string>::iterator dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            
            if(dataSetIt==m_partitionFunctionDataSetFolders.end())
            {
                // Look for original configuration
                chargeConfiguration = chargeConfiguration.getAntiConfiguration();
                dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            }
            
            if(dataSetIt!=m_partitionFunctionDataSetFolders.end())
            {
                // Load new data set and build partition function subset
                updatePartitionFunctionSubSet(dataSetIt->first,dataSetIt->second);
            }
            else
            {
                // Required data not found in provided partition function data set
                io_calculationErrorStatus = 1;
                return 0.;
            }
        }
        
        // Perform 1D linear interpolation
        PartitionFunctionMassGridData* massGridData = &m_partitionFunctionMassGrids[chargeConfiguration];
        bool interpolationError(false);
        const SingleParameterInterpolationData massInterpData(computeInterpolationData(massGridData->mass,i_mass,interpolationError));
        
        // Check interpolation error status
        if(interpolationError)
        {            
            // Interpolation mass value outside available range
            io_calculationErrorStatus = 2;
            return 0.;
        }
        
        // Interpolation parameter calculation correctly performed
        io_calculationErrorStatus = 0;
        
        return massGridData->partitionFunction[massInterpData.gridPoints.first]*massInterpData.weights.first +
               massGridData->partitionFunction[massInterpData.gridPoints.second]*massInterpData.weights.second;

    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

// TODO: cluster mass value seems to differ wrt reference one at some precision level, causing a different value in mass interp. coeff.
double PartitionFunctionHandling::getPartitionFunctionValueErr(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                            const double i_mass,
                                                            unsigned int& io_calculationErrorStatus)
{
    try
    {
        MCSTHAR::Utilities::ChargeConfiguration chargeConfiguration(i_chargeConfiguration);
        
        // Look for provided charge configuration
        map<MCSTHAR::Utilities::ChargeConfiguration,PartitionFunctionMassGridData>::const_iterator massGridIt =
        m_partitionFunctionMassGrids.find(chargeConfiguration);
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Look for anti configuration
            chargeConfiguration = chargeConfiguration.getAntiConfiguration();
            massGridIt = m_partitionFunctionMassGrids.find(chargeConfiguration);
        }
        
        if(massGridIt==m_partitionFunctionMassGrids.end())
        {
            // Data currently not available for the provided charge configuration.
            // Look for required data set existence
            map<MCSTHAR::Utilities::ChargeConfiguration,string>::iterator dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            
            if(dataSetIt==m_partitionFunctionDataSetFolders.end())
            {
                // Look for original configuration
                chargeConfiguration = chargeConfiguration.getAntiConfiguration();
                dataSetIt = m_partitionFunctionDataSetFolders.find(chargeConfiguration);
            }
            
            if(dataSetIt!=m_partitionFunctionDataSetFolders.end())
            {
                // Load new data set and build partition function subset
                updatePartitionFunctionSubSet(dataSetIt->first,dataSetIt->second);
            }
            else
            {
                // Required data not found in provided partition function data set
                io_calculationErrorStatus = 1;
                return 0.;
            }
        }
        
        // Perform 1D linear interpolation
        PartitionFunctionMassGridData* massGridData = &m_partitionFunctionMassGrids[chargeConfiguration];
        bool interpolationError(false);
        const SingleParameterInterpolationData massInterpData(computeInterpolationData(massGridData->mass,i_mass,interpolationError));
        
        // Check interpolation error status
        if(interpolationError)
        {
            // Interpolation mass value outside available range
            io_calculationErrorStatus = 2;
            return 0.;
        }
        
        // Interpolation parameter calculation correctly performed
        io_calculationErrorStatus = 0;
        
        return massGridData->partitionFunctionErr[massInterpData.gridPoints.first]*massInterpData.weights.first +
        massGridData->partitionFunctionErr[massInterpData.gridPoints.second]*massInterpData.weights.second;
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}


void PartitionFunctionHandling::updatePartitionFunctionSubSet(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                              const string& i_dataFolder)
{
    // Retrieve data folder file content
    // TODO: patch for correct file ordering (Linux), to be removed with
    // final partition function loading strategy
    vector<string> dataFileList(getDirectoryContent(i_dataFolder));
    unsigned int dataFileNumber(dataFileList.size());
    // const vector<string> dataFileList(getDirectoryContent(i_dataFolder));
    // const unsigned int dataFileNumber(dataFileList.size());
    // TODO: end patch for correct file ordering (Linux)

    if(dataFileNumber==0)
    {
        throw HadronizationException("Error during single charge configuration partition function data set loading, no file found",
                                     __FUNCTION__,132);
    }
    
    try
    {
        vector<double> partitionFunctionGrid;
        vector<double> partitionFunctionGridErr;
        const string filePath(i_dataFolder + "/");
        
        // TODO: patch for correct file ordering (Linux), to be removed with
        // final partition function loading strategy
        // Build ordered list
        map<unsigned int,string> fileMap;
        for(unsigned int fileIndex = 0;fileIndex<dataFileNumber;++fileIndex)
        {
            if(dataFileList[fileIndex].size()==8)
            {
                stringstream fileNumberStr(dataFileList[fileIndex]);
                unsigned int fileNumber;
                fileNumberStr>>fileNumber;
                fileMap[fileNumber] = dataFileList[fileIndex];
            }
        }
        // Fill data file list vector wit ordered map
        dataFileNumber = fileMap.size();
        dataFileList.resize(dataFileNumber);
        const map<unsigned int,string>::const_iterator fileMapItEnd(fileMap.end());
        unsigned int fileIndex(0);
        for(map<unsigned int,string>::iterator fileMapIt = fileMap.begin();
            fileMapIt != fileMapItEnd;++fileMapIt,++fileIndex)
        {
            dataFileList[fileIndex] = fileMapIt->second;
        }
        // TODO: end patch for correct file ordering (Linux)
        
        
        // Loop over data folder file set
        for(unsigned int fileIndex = 0;fileIndex<dataFileNumber;++fileIndex)
        {
            // Skip non data file
            // TODO: remove after file list file usage
            if(dataFileList[fileIndex].size()==8)
            {
                // TODO: Extend to Windows case
                string fullFileName(filePath+dataFileList[fileIndex]);
                
                // Load partition function data file
                loadPartitionFunctionDataFile(fullFileName,i_chargeConfiguration,partitionFunctionGrid,partitionFunctionGridErr);
            }
        }
        
        // Update mass grid data for the considered charge configuration
        m_partitionFunctionMassGrids[i_chargeConfiguration].massValueNumber = dataFileNumber;
        
        // Compute interpolation data if required
        if(m_isInterpolationDataAvailable==false)
        {
            // TODO: does a faster STL algorithm exist?
            // Remove duplicated energy density grid parameter values
            sort(m_energyDensityValues.begin(),m_energyDensityValues.end());
            double referenceValue = m_energyDensityValues[0];
            for(unsigned int index = 1;index<m_energyDensityValues.size();++index)
            {
                if(m_energyDensityValues[index]==referenceValue)
                {
                    m_energyDensityValues.erase(m_energyDensityValues.begin()+index);
                    --index;
                }
                else
                {
                    referenceValue = m_energyDensityValues[index];
                }
            }
                
            // Remove duplicated gammaS grid parameter values
            sort(m_gammaSValues.begin(),m_gammaSValues.end());
            referenceValue = m_gammaSValues[0];
            for(unsigned int index = 1;index<m_gammaSValues.size();++index)
            {
                if(m_gammaSValues[index]==referenceValue)
                {
                    m_gammaSValues.erase(m_gammaSValues.begin()+index);
                    --index;
                }
                else
                {
                    referenceValue = m_gammaSValues[index];
                }
            }
            
            bool interpolationErrorStatus(false);
            const SingleParameterInterpolationData energyInterpDensityData(computeInterpolationData(m_energyDensityValues,
                                                                                                    m_energyDensity,
                                                                                                    interpolationErrorStatus));
            // Check interpolation error status
            if(interpolationErrorStatus)
            {
                throw HadronizationException("Error during interpolation data calculation (energy density), interpolation value outside the available range",
                                             __FUNCTION__,134);
            }
            
            const SingleParameterInterpolationData strangenessSuppressionInterpData(computeInterpolationData(m_gammaSValues,
                                                                                                             m_gammaS,
                                                                                                             interpolationErrorStatus));
            // Check interpolation error status
            if(interpolationErrorStatus)
            {
                throw HadronizationException("Error during interpolation data calculation (gamma_S), interpolation value outside the available range",
                                             __FUNCTION__,134);
            }

            // Compute 2D grid interpolation indexes
            m_energyDensityPointNumber = m_energyDensityValues.size();
            m_grid2DDimension = m_energyDensityPointNumber*m_gammaSValues.size();
            m_index00 = energyInterpDensityData.gridPoints.first +
                        m_energyDensityPointNumber*strangenessSuppressionInterpData.gridPoints.first;
            m_index10 = m_index00+1;
            m_index01 = m_index00+m_energyDensityPointNumber;
            m_index11 = m_index01+1;
            
            // Compute 2D grid interpolation weights
            m_weight00 = energyInterpDensityData.weights.first*strangenessSuppressionInterpData.weights.first;
            m_weight10 = energyInterpDensityData.weights.second*strangenessSuppressionInterpData.weights.first;
            m_weight01 = energyInterpDensityData.weights.first*strangenessSuppressionInterpData.weights.second;
            m_weight11 = energyInterpDensityData.weights.second*strangenessSuppressionInterpData.weights.second;
            
            m_isInterpolationDataAvailable=true;
        }
        
        // Build partition function subset
        buildPartitionFunctionSubSet(i_chargeConfiguration,partitionFunctionGrid,partitionFunctionGridErr);

    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void PartitionFunctionHandling::loadPartitionFunctionDataFile(const string& i_dataFileName,
                                                              const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                              vector<double>& io_partitionFunctionGrid,
                                                              vector<double>& io_partitionFunctionGridErr)
{
    // TODO: QUESTO METODO E' DA RIPULIRE DOPO AVER DECISO COME GESTIRE IL DATA SET
    // Open data file
    ifstream dataFile;
    dataFile.open(i_dataFileName.c_str());
    if(dataFile.is_open())
    {
        // TODO: add check for data file contents (nan, inf, ect)
        unsigned int numberOfIntegrationIterations;
        double energyDensitySamplingValue;
        
        // TODO: differenza tra sampling mass e mass value???
        double massSamplingValue;
        
        // TODO: esiste funzione per contare le linee prima della lettura???
        unsigned int numberOfFileLines = 0;
        double mass;
        double strangenessSuppression;
        double energyDensity;
        double partitionFunctionValue;
        double partitionFunctionVarianceValue;
        vector<double> partitionFunctionGrid;
            
        // Load data file
        dataFile>>numberOfIntegrationIterations;
        dataFile>>energyDensitySamplingValue>>massSamplingValue;
        for(;!(dataFile.eof());++numberOfFileLines)
        {
            // TODO: questa è una stronzata
            mass = strangenessSuppression = energyDensity = partitionFunctionValue = partitionFunctionVarianceValue = 0.;
            
            // TODO: la varianza sarebbe da eliminare visto il livello di precisione...
            dataFile>>mass>>strangenessSuppression>>energyDensity>>partitionFunctionValue>>partitionFunctionVarianceValue;
            if(dataFile.eof())
            {
                break;
            }
            
            // TODO: find a different strategy for this check (eof-like)
            if(mass==0.)
            {
                string errorMessage = "Error during " + i_dataFileName + " partition function data file loading, corrupted file";
                throw HadronizationException(errorMessage,__FUNCTION__,133);
            }
            
            if(numberOfFileLines==0)
            {
                m_partitionFunctionMassGrids[i_chargeConfiguration].mass.push_back(mass);
            }
            io_partitionFunctionGrid.push_back(partitionFunctionValue);
            io_partitionFunctionGridErr.push_back(partitionFunctionVarianceValue);
            
            if(m_is1DGridDataAvailable==false)
            {
                m_energyDensityValues.push_back(energyDensity);
                m_gammaSValues.push_back(strangenessSuppression);
            }
        }
        
        if(m_is1DGridDataAvailable==false)
        {
            m_is1DGridDataAvailable = true;
        }
        
        // Check file structure
        if(numberOfFileLines==0)
        {
            const string errorMessage("Error during " + i_dataFileName +
                                      " partition function data file loading, missing partition function values");
            throw HadronizationException(errorMessage,__FUNCTION__,133);
        }
        dataFile.close();
    }
    else
    {
        const string errorMessage("Error during " + i_dataFileName + " partition function data file loading");
        throw HadronizationException(errorMessage,__FUNCTION__,133);
    }
}

void PartitionFunctionHandling::buildPartitionFunctionSubSet(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                             const vector<double>& i_partitionFunctionGrid,
                                                             const vector<double>& i_partitionFunctionErrGrid)
{
    
    // Retrieve partition function subset for the provided charge configuration
    map<MCSTHAR::Utilities::ChargeConfiguration,PartitionFunctionMassGridData>::iterator massGridIt =
        m_partitionFunctionMassGrids.find(i_chargeConfiguration);

    // Loop over mass grid points
    const unsigned int massPointNumber(massGridIt->second.massValueNumber);
    massGridIt->second.partitionFunction.resize(massPointNumber);
    massGridIt->second.partitionFunctionErr.resize(massPointNumber);
    
    for(unsigned int massIndex=0;massIndex<massPointNumber;++massIndex)
    {
        const unsigned int offsetIndex(m_grid2DDimension*massIndex);
        
//        cout<<"massIndex "<<massIndex
//            <<" i_partitionFunctionGrid.at(offsetIndex+m_index00) "<<i_partitionFunctionGrid.at(offsetIndex+m_index00)
//            <<" i_partitionFunctionGrid.at(offsetIndex+m_index10) "<<i_partitionFunctionGrid.at(offsetIndex+m_index10)
//            <<" i_partitionFunctionGrid.at(offsetIndex+m_index01) "<<i_partitionFunctionGrid.at(offsetIndex+m_index01)
//        <<" i_partitionFunctionGrid.at(offsetIndex+m_index11) "<<i_partitionFunctionGrid.at(offsetIndex+m_index11)<<endl;
        
        massGridIt->second.partitionFunction[massIndex] = i_partitionFunctionGrid.at(offsetIndex+m_index00)*m_weight00 +
                                                          i_partitionFunctionGrid.at(offsetIndex+m_index10)*m_weight10 +
                                                          i_partitionFunctionGrid.at(offsetIndex+m_index01)*m_weight01 +
                                                          i_partitionFunctionGrid.at(offsetIndex+m_index11)*m_weight11;
        massGridIt->second.partitionFunctionErr[massIndex] = i_partitionFunctionErrGrid.at(offsetIndex+m_index00)*m_weight00 +
                                                          i_partitionFunctionErrGrid.at(offsetIndex+m_index10)*m_weight10 +
                                                          i_partitionFunctionErrGrid.at(offsetIndex+m_index01)*m_weight01 +
                                                          i_partitionFunctionErrGrid.at(offsetIndex+m_index11)*m_weight11;
    }
    return;
}

// TODO: does a faster interpolation coefficients/points exist?
SingleParameterInterpolationData PartitionFunctionHandling::computeInterpolationData(const vector<double>& i_valueSet,
                                                                                     const double i_interpolationValue,
                                                                                     bool& io_interpolationErrorStatus) const
{
    // Check inclusion of provided interpolation point wrt available range
    if(i_interpolationValue<*(i_valueSet.begin()) || i_interpolationValue>*(i_valueSet.end()-1))
    {
        io_interpolationErrorStatus = true;
        return SingleParameterInterpolationData();
    }
    
    // Find interpolation point right side grid point
    double rightBound;
    const unsigned int numberOfPoints(i_valueSet.size());
    unsigned int pointIndex;
    for(pointIndex = 0;pointIndex<numberOfPoints;++pointIndex)
    {
        rightBound = i_valueSet[pointIndex];
        if(i_interpolationValue<=rightBound)
        {
            break;
        }
    }
    
    // Build interpolation data
    SingleParameterInterpolationData o_interpolationData;
    if(pointIndex>0)
    {
        o_interpolationData.gridPoints = pair<unsigned int,unsigned int>(pointIndex-1,pointIndex);
        const double leftBound(i_valueSet[pointIndex-1]);
        o_interpolationData.weights.second = (i_interpolationValue - leftBound)/(rightBound - leftBound);
    }
    else
    {
        o_interpolationData.gridPoints = pair<unsigned int,unsigned int>(0,0);
        o_interpolationData.weights.second = 0.;
    }
    o_interpolationData.weights.first = 1. - o_interpolationData.weights.second;
    
    // Interpolation correctly performed
    io_interpolationErrorStatus = false;
    
    return o_interpolationData;
}

// TODO: move to utils togheter with its opposite method (in partition function writer class)
MCSTHAR::Utilities::ChargeConfiguration PartitionFunctionHandling::retrieveChargeConfiguration(const string& i_dataFolderName) const
{
    // TODO: add check for othe possible errors and find final strategy
    if(i_dataFolderName.size() != 10)
    {
        const string errorMessage = "Error during " +
                                    i_dataFolderName +
                                    " partition function data file loading, unidentified charge configuration";
        throw HadronizationException(errorMessage,__FUNCTION__,135);
    }
    
    try
    {
        MCSTHAR::Utilities::ChargeConfiguration o_chargeConfiguration;
        int chargeValue;
        
        // Retrieve strange charge
        retrieveChargeData(i_dataFolderName.substr(0,2),chargeValue);
        o_chargeConfiguration.strangeCharge = chargeValue;
        
        // Retrieve electric charge
        retrieveChargeData(i_dataFolderName.substr(2,2),chargeValue);
        o_chargeConfiguration.electricCharge = static_cast<double>(chargeValue);

        // Retrieve baryonic charge
        retrieveChargeData(i_dataFolderName.substr(4,2),chargeValue);
        o_chargeConfiguration.baryonicCharge = static_cast<double>(chargeValue);

        // Retrieve charm charge
        retrieveChargeData(i_dataFolderName.substr(6,2),chargeValue);
        o_chargeConfiguration.charmCharge = chargeValue;

        // Retrieve bottom charge
        retrieveChargeData(i_dataFolderName.substr(8,2),chargeValue);
        o_chargeConfiguration.bottomCharge = chargeValue;

        return o_chargeConfiguration;
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void PartitionFunctionHandling::retrieveChargeData(const string& i_chargeDataString,int& io_chargeValue) const
{
    // Single charge string is a 2 char string XY with X storing charge sign and Y charge value
    // TODO: find final strategy
    if(i_chargeDataString.size()!=2)
    {
        throw HadronizationException("Error during data folder charge configuration retrieval, unidentified charge configuration",
                                     __FUNCTION__,135);
    }
    
    stringstream chargeData;
    chargeData<<i_chargeDataString.substr(1,1);
    chargeData>>io_chargeValue;
    
    int chargeDataNumber;
    chargeData.clear();
    chargeData<<i_chargeDataString.substr(0,1);
    chargeData>>chargeDataNumber;
    
    if(chargeDataNumber == 1)
    {
        io_chargeValue *= -1;
    }
    else
    {
        if(chargeDataNumber != 0)
        {
            throw HadronizationException("Error during data folder charge configuration retrieval, unidentified charge configuration",
                                         __FUNCTION__,135);
        }
    }
    
    return;
}

// TODO: check and find a more elegant solution...
bool MCSTHAR::Utilities::ChargeConfiguration::operator<(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration) const
{
    if(strangeCharge<i_chargeConfiguration.strangeCharge)
    {
        return true;
    }
    else if(strangeCharge==i_chargeConfiguration.strangeCharge)
    {
        if(charmCharge<i_chargeConfiguration.charmCharge)
        {
            return true;
        }
        else if(charmCharge==i_chargeConfiguration.charmCharge)
        {
            if(bottomCharge<i_chargeConfiguration.bottomCharge)
            {
                return true;
            }
            else if(bottomCharge==i_chargeConfiguration.bottomCharge)
            {
                if(electricCharge<i_chargeConfiguration.electricCharge)
                {
                    return true;
                }
                else if(electricCharge==i_chargeConfiguration.electricCharge)
                {
                    if(baryonicCharge<i_chargeConfiguration.baryonicCharge)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                    
                }
                else
                {
                    return false;
                }
                
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}
