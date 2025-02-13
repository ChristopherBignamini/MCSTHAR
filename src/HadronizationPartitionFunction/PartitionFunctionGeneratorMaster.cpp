#include "../Include/PartitionFunctionGeneratorMaster.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/logMessage.h"
#include "../../../Utilities/Include/Constants.h"
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;

// TODO: check everywhere exit code
// TODO: check everywhere input parameter constness

PartitionFunctionGeneratorMaster::PartitionFunctionGeneratorMaster(const double i_minGammaS,
                                                                   const double i_maxGammaS,
                                                                   const unsigned int i_numberOfGammaSValues,
                                                                   const HadronSamplingGroups& i_hadronSamplingGroups,
                                                                   const unsigned long long int i_maxChannelSamplingNumber,
                                                                   const double i_integrationErrorThreshold,
                                                                   const unsigned int i_randomNumberGeneratorSeed,
                                                                   const MPI_Comm& i_mpiCommunicator)
                                                                  :PartitionFunctionGenerator(i_minGammaS,
                                                                                              i_maxGammaS,
                                                                                              i_numberOfGammaSValues,
                                                                                              i_hadronSamplingGroups,
                                                                                              i_randomNumberGeneratorSeed,
                                                                                              i_mpiCommunicator)
                                                                  ,m_isPartitionFunctionAvailable(false)
                                                                  ,m_partitionFunctionData(i_numberOfGammaSValues)
                                                                  ,m_cumulativeWeightSumData(NULL)
                                                                  ,m_totalMaxSamplingNumber(i_maxChannelSamplingNumber)// TODO: Store somehwere with other analogous parameters?
                                                                  ,m_localMaxSamplingNumber(m_totalMaxSamplingNumber)
                                                                  ,m_minimumStabilityCheckSamplingNumber(1000)//100000
                                                                  ,m_integrationSquaredErrorThreshold(i_integrationErrorThreshold*
                                                                                                      i_integrationErrorThreshold)
                                                                  ,m_numberOfTasks(1)
{
    if(i_maxChannelSamplingNumber==0)
    {
        throw HadronizationException("Error in partition function calculation setup, maximum number of channel samplings set to zero",
                                     __FUNCTION__,
                                     861);
    }
    
    // Check provided number of minimum and maximum samplings
    if(m_totalMaxSamplingNumber<m_minimumStabilityCheckSamplingNumber)
    {
        // Throw exception for bad sampling number configuration provided
        throw HadronizationException(
                                     "Error in partition function calculation setup, minimum number of channel samplings larger than the provided maximum one",
                                     __FUNCTION__,862);
    }
    
    if(i_integrationErrorThreshold<=0.)
    {
        throw HadronizationException("Error in partition function calculation setup, negative error integration relative error threshold provided",
                                     __FUNCTION__,
                                     863);
    }

    // Retrieve total number of tasks
    int exitStatus(MPI_Comm_size(m_mpiCommunicator,
                                 &m_numberOfTasks));
    if(exitStatus!=0)
    {
        throw HadronizationException("Error during number of tasks retrieval",
                                     __FUNCTION__,
                                     864);
    }

    // Compute local integration parameters
    if(m_numberOfTasks>1)
    {
        // Compute single task minimum sampling number
        if(m_minimumStabilityCheckSamplingNumber<m_numberOfTasks)
        {
            // If number of minimum sampling number is smaller than number of tasks reset it to number of threads itself
            m_minimumStabilityCheckSamplingNumber = m_numberOfTasks;
            m_localStabilityCheckSamplingNumberStart = 1;
            
            // Log warning for bad sampling number configuration
            stringstream numberOfSamplingsStr;
            numberOfSamplingsStr<<m_minimumStabilityCheckSamplingNumber;
            logMessage("Minimum number of hadronization channel (full) configuration samplings reset to " +
                       numberOfSamplingsStr.str() +
                       "(one per task).",
                       WARNING);
        }
        else
        {
            // Check if minimum sampling number can be evenly distributed among tasks, otherwise recompute it
            const unsigned long long int samplingNumberPerTask(m_minimumStabilityCheckSamplingNumber/m_numberOfTasks);
            if((samplingNumberPerTask*m_numberOfTasks)!=m_minimumStabilityCheckSamplingNumber)
            {
                m_localStabilityCheckSamplingNumberStart = samplingNumberPerTask + 1;
                m_minimumStabilityCheckSamplingNumber = m_localStabilityCheckSamplingNumberStart*m_numberOfTasks;
                
                // Log warning for unbalanced sampling number configuration reset
                stringstream numberOfSamplingsStr;
                stringstream numberOfSamplingsPerTaskStr;
                numberOfSamplingsStr<<m_minimumStabilityCheckSamplingNumber;
                numberOfSamplingsPerTaskStr<<m_localStabilityCheckSamplingNumberStart;
                logMessage("Minimum number of hadronization channel (full) configuration samplings reset to " +
                           numberOfSamplingsStr.str() +
                           " (" +
                           numberOfSamplingsPerTaskStr.str() +
                           " per task).",
                           WARNING);
                
            }
            else
            {
                m_localStabilityCheckSamplingNumberStart = samplingNumberPerTask;
            }
        }
                
        // Compute single task maximum sampling number
        if(m_totalMaxSamplingNumber<m_numberOfTasks)
        {
            // If number of maximum sampling number is smaller than number of tasks reset it to number of tasks itself
            m_totalMaxSamplingNumber = m_numberOfTasks;
            
            // Log warning for bad sampling number configuration
            stringstream numberOfSamplingsStr;
            numberOfSamplingsStr<<m_totalMaxSamplingNumber;
            logMessage("Maximum number of hadronization channel (full) configuration samplings reset to " +
                       numberOfSamplingsStr.str() +
                       "(one per task). Bad task configuration provided.",
                       WARNING);// TODO: throw exception in this case
        }
        else
        {
            // Check if maximum sampling number can be evenly distributed among threads, otherwise recompute it
            const unsigned long long int samplingNumberPerTask(m_totalMaxSamplingNumber/m_numberOfTasks);
            if((samplingNumberPerTask*m_numberOfTasks)!=m_totalMaxSamplingNumber)
            {
                m_localMaxSamplingNumber = samplingNumberPerTask+1;
                m_totalMaxSamplingNumber = m_localMaxSamplingNumber*m_numberOfTasks;
                
                // Log warning for unbalanced sampling number configuration reset
                stringstream numberOfSamplingsStr;
                stringstream numberOfSamplingsPerTaskStr;
                numberOfSamplingsStr<<m_totalMaxSamplingNumber;
                numberOfSamplingsPerTaskStr<<m_localMaxSamplingNumber;
                logMessage("Maximum number of hadronization channel (full) configuration samplings reset to " +
                           numberOfSamplingsStr.str() +
                           " (" +
                           numberOfSamplingsPerTaskStr.str() +
                           " per task).",
                           WARNING);
            }
            else
            {
                m_localMaxSamplingNumber = samplingNumberPerTask;
            }
        }
    }
    else
    {
        m_localStabilityCheckSamplingNumberStart = m_minimumStabilityCheckSamplingNumber;
    }
    
    // Send to slave tasks stability check sampling number data
    exitStatus = MPI_Bcast(&m_localStabilityCheckSamplingNumberStart,
                           1,
                           MPI_UNSIGNED_LONG_LONG,
                           0,
                           m_mpiCommunicator);    
    // Set data members
    if(exitStatus==0)
    {
        m_isPartitionFunctionGeneratorReady = true;
    }
    else
    {
        throw HadronizationException("Error during stability check data broadcasting",
                                     __FUNCTION__,
                                     865);
    }
    
    // Allocate memory for partition function calculation
    m_cumulativeWeightSumData = new double[m_weightSumDataSize];// TODO: add allocation success
    
    // TODO: temporary!
    m_hadronSamplingNew.m_verboseLog = true;
    
}

PartitionFunctionGeneratorMaster::~PartitionFunctionGeneratorMaster(void)
{
    if(m_cumulativeWeightSumData != NULL)
    {
        delete [] m_cumulativeWeightSumData;
        m_cumulativeWeightSumData = NULL;
    }
}

const SingleMassEnergyDensityPairPartitionFunctionData& PartitionFunctionGeneratorMaster::getPartitionFunctionData(void) const
{
    if(m_isPartitionFunctionAvailable)
    {
        return m_partitionFunctionData;
    }
    else
    {
        throw HadronizationException("Error during partition function data retrieval, partition function data not available",
                                     __FUNCTION__,869);
    }
}

// TODO: remove the energy density input parameter with the new m_hadronSampling interface
// TODO: refactoring required
void PartitionFunctionGeneratorMaster::runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                                       const double i_mass,
                                                                       const double i_energyDensity)
{
    // Reset global cumulative calculation data
    fill(m_cumulativeWeightSumData,m_cumulativeWeightSumData+m_weightSumDataSize,0.);
    
    // TODO: unsigned long int or unsigned long long int? Check counters everywhere
    // TODO: Is there a better calculation method to avoid huge numbers?
    // TODO: debug! unsigned long int or unsigned long long int? Is there a better calculation method to avoid huge numbers?
    unsigned long long int stabilityStatusPrintSamplingNumber(1000);//100000// TODO: remove
    unsigned long long int cumulativeLocalSamplingNumber(0);
        
    // TODO: check alternative strategy (no communication)
    m_isPartitionFunctionAvailable = false;
    bool newIterationsRequired(true);
    double maximumRelativeError(-1.);
    unsigned int dataIndex(0);
    unsigned long long int totalSamplingNumber(0);
    bool doHadronizationChannelsExist(true);
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
        
        // Check stability status
        if(m_cumulativeWeightSumData[0] > 0.)// TODO: switch to MPI derived datatype
        {
            // Update total sampling numbers
            cumulativeLocalSamplingNumber += m_localStabilityCheckSamplingNumber;
            totalSamplingNumber = cumulativeLocalSamplingNumber*m_numberOfTasks;
            // Loop over strangeness suppression parameter grid points and update weight values
            // TODO: avoid this loop by passing a weight vector to computeStrangenessSuppression method
            maximumRelativeError = -1.;
            dataIndex = 0;
            for(unsigned int gammaSIndex=0;gammaSIndex<m_numberOfGammaSValues;
                ++gammaSIndex,dataIndex+=2)
            {
                // Reset stability status for current strangeness supression parameter value
                // (the state "true" could have been set in a previous check)
                m_partitionFunctionData.partitionFunctionErrorStatus[gammaSIndex] = false;
                
                // Check stability status for current strangeness suppression parameter value
                const double& weightSumElement(m_cumulativeWeightSumData[dataIndex]);
                const double& squaredWeightSumElement(m_cumulativeWeightSumData[dataIndex+1]);
                
                // TODO: can the usage of squared error and squared sum be avoided?
                // TODO: add numerical check for /0, etc...
                // TODO: optimize
                // TODO: add comments on calculations
                // Check if channel exists (//TODO: use flag instead)
                if(weightSumElement>0.)
                {
                    const double oneOnWeightSumElementSquared(1./(weightSumElement*weightSumElement));
                    double channelSumSquaredRelativeError((totalSamplingNumber*squaredWeightSumElement*oneOnWeightSumElementSquared - 1.)/
                                                          (totalSamplingNumber-1));
                    
                    // Check for numerical relative squared error non positiveness
                    if(channelSumSquaredRelativeError<0.)
                    {                        
                        if(-channelSumSquaredRelativeError<=partitionFunctionSquaredErrorZeroRange)
                        {
                            // Relative squared error is approximatively zero, set it to exactly zero
                            channelSumSquaredRelativeError = 0.;
                        }
                        else
                        {
                            // Relative squared error is negative, throw exception
                            throw HadronizationException("Error during channel sampling: negative channel sum relative squared error",
                                                         __FUNCTION__,866);
                        }
                    }
                    
                    // Update maximum relative error
                    if(channelSumSquaredRelativeError>maximumRelativeError)
                    {
                        maximumRelativeError = channelSumSquaredRelativeError;
                    }
                    
                    // TODO: check on relative or absolute error?
                    // TODO: check denominator
                    // Check stability status
                    // Check if unstable channel sum exist
                    if(channelSumSquaredRelativeError<m_integrationSquaredErrorThreshold)
                    {
                        // Stability reached for current strangeness suppression parameter value
                        m_partitionFunctionData.partitionFunctionErrorStatus[gammaSIndex] = true;
                    }
                }
                else
                {
                    throw HadronizationException("Error during partition function calculation, non positive weight sum",
                                                 __FUNCTION__,
                                                 867);
                }
            }

            // Print stability status
//            if(totalSamplingNumber==stabilityStatusPrintSamplingNumber)
//            {
                // Create log message
                // TODO: store log string as const static (?);
                stringstream totalSamplingNumberStr;
                stringstream maximumRelativeErrorValueStr;
                totalSamplingNumberStr<<totalSamplingNumber;
                maximumRelativeErrorValueStr<<sqrt(maximumRelativeError);
                // TODO: use function for message pattern creation
                const string message("Maximum relative error " +
                                     maximumRelativeErrorValueStr.str() +
                                     " after " +
                                     totalSamplingNumberStr.str() +
                                     " hadronization channel samplings.");
                // Log message
                logMessage(message);
                
//                
//                stabilityStatusPrintSamplingNumber += 10000;// TODO: find final solution
//            }
            
            // No errors/anomalies present, check if the maximum number of sampling number or stability has been reached
            if(find(m_partitionFunctionData.partitionFunctionErrorStatus.begin(),
                    m_partitionFunctionData.partitionFunctionErrorStatus.end(),
                    false) == m_partitionFunctionData.partitionFunctionErrorStatus.end())
            {
                // Stability reached for all strangeness suppression parameter values, break while loop
                newIterationsRequired = false;
                m_localStabilityCheckSamplingNumber = 0;
            }
            else
            {
                // TODO: some refactoring required among all the performed checks
                // Check if maximum number of allowed iterations has been reached
                if(cumulativeLocalSamplingNumber==m_localMaxSamplingNumber)
                {
                    newIterationsRequired = false;
                    m_localStabilityCheckSamplingNumber = 0;
                }
                else
                {
                    // Stability not reached for all strangeness suppression parameter values, set next stability check sampling number
                    // TODO: if this is final strategy remove code below
                    const unsigned long int nextLocalSamplingNumberCheck(cumulativeLocalSamplingNumber+
                                                                         m_localStabilityCheckSamplingNumberStart);
                    if(nextLocalSamplingNumberCheck<=m_localMaxSamplingNumber)
                    {
                        m_localStabilityCheckSamplingNumber = m_localStabilityCheckSamplingNumberStart;
                    }
                    else
                    {
                        m_localStabilityCheckSamplingNumber = m_localMaxSamplingNumber - cumulativeLocalSamplingNumber;
                    }
                }
            }
        }
        else
        {
            // Apparently no hadronization channels exist for the provided cluster mass, energy density
            // and charge configuration and for the available hadron set: set partition function
            // to zero and log warning
            resetPartitionFunctionData();
            newIterationsRequired = false;
            m_localStabilityCheckSamplingNumber = 0;
            doHadronizationChannelsExist = false;
            string warningMsg("No hadronization channel found. Partition function values are set to zero for current cluster mass value. Check hadron list file for channel availability");
            logMessage(warningMsg,WARNING);            
        }
        
        // Master task sends calculation status check flag to slaves
        const int exitStatus(MPI_Bcast(&m_localStabilityCheckSamplingNumber,// TODO: use data member for exit status??
                                       1,
                                       MPI_UNSIGNED_LONG_LONG,
                                       0,
                                       m_mpiCommunicator));
                
        if(exitStatus!=0)
        {
            // TODO: add MPI error status
            throw HadronizationException("Error during stability check data broadcasting",
                                         __FUNCTION__,
                                         868);
        }
    }
    
    // Set partition function data
    // Set grand canonical ensemble parameter used by channel sampling function
    m_partitionFunctionData.grandCanonicalParameters = GrandCanonicalParameters(m_hadronSamplingNew.getSamplingElectricChargeFugacity(),
                                                                                m_hadronSamplingNew.getSamplingStrangeChargeFugacity(),
                                                                                m_hadronSamplingNew.getSamplingBaryonicChargeFugacity(),
                                                                                m_hadronSamplingNew.getSamplingCharmChargeFugacity(),
                                                                                m_hadronSamplingNew.getSamplingBottomChargeFugacity(),
                                                                                m_hadronSamplingNew.getSamplingTemperature());
    if(doHadronizationChannelsExist)
    {
        dataIndex = 0;
        for(unsigned int gammaSIndex=0;gammaSIndex<m_numberOfGammaSValues;++gammaSIndex,dataIndex+=2)
        {
            // TODO: optimize calculation below
            // TODO: check unsigned int
            const double& weightSumElement(m_cumulativeWeightSumData[dataIndex]);
            
            // Compute partition function valie
            m_partitionFunctionData.partitionFunctions[gammaSIndex] = weightSumElement/totalSamplingNumber;
            
            // Check for numerical squared error non positiveness
            double partitionFunctionSquaredError((m_cumulativeWeightSumData[dataIndex+1]-
                                                  weightSumElement*weightSumElement/totalSamplingNumber)/
                                                 (totalSamplingNumber*(totalSamplingNumber-1)));
            
            if(partitionFunctionSquaredError<0.)
            {
                if(-partitionFunctionSquaredError<=partitionFunctionSquaredErrorZeroRange)
                {
                    // Squared error is approximatively zero, set it to exactly zero
                    partitionFunctionSquaredError = 0.;
                }
                else
                {                    
                    // Squared error is negative, throw exception
                    throw HadronizationException("Error during channel sampling: negative channel sum squared error",
                                                 __FUNCTION__,867);
                }
            }
            
            // Compute partition function error
            m_partitionFunctionData.partitionFunctionErrors[gammaSIndex] = sqrt(partitionFunctionSquaredError);
            m_partitionFunctionData.numberOfChannelSamplings[gammaSIndex] = totalSamplingNumber;

        }        
    }
    
    m_isPartitionFunctionAvailable = true;
    return;
}

void PartitionFunctionGeneratorMaster::resetPartitionFunctionData(void)
{
    // TODO: factorize with standard partition function data assignement
    // TODO: zero reassignement in case of subsequent empty channel sets
    
    // Set to zero the partition function values and the corresponding error and number of sampled channels. Set to true the
    // their stability status, assuming that no hadronization channel is available for the considered mass, charge configuration
    // and hadron set
    m_partitionFunctionData.partitionFunctions.assign(m_numberOfGammaSValues,0.);
    m_partitionFunctionData.partitionFunctionErrors.assign(m_numberOfGammaSValues,0.);
    m_partitionFunctionData.partitionFunctionErrorStatus.assign(m_numberOfGammaSValues,true);
    m_partitionFunctionData.numberOfChannelSamplings.assign(m_numberOfGammaSValues,0);
    
    return;
}

void PartitionFunctionGeneratorMaster::mergeLocalSums(void)
{
    
    int exitStatus(MPI_Reduce(MPI_IN_PLACE,
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
                                     868);
    }
    
    // Update global cumulative weight and squared weight sum
    for(unsigned int dataIndex=0;dataIndex<m_weightSumDataSize;++dataIndex)
    {
        m_cumulativeWeightSumData[dataIndex] += m_weightSumData[dataIndex];// TODO: find a way to avoid copy
    }
    
    return;
}