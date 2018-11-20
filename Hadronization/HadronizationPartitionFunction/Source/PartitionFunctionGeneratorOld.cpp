#include "../Include/PartitionFunctionGeneratorOld.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/logMessage.h"
#include "../../../Utilities/Include/Constants.h"
#include <cassert>
#include <algorithm>

// TODO: debug
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

// TODO: remove!
bool OrderedHadronIdsComparison::operator()(const vector<int>& i_firstChannelOrderedHadronIds,
                                            const vector<int>& i_secondChannelOrderedHadronIds) const
{
    if(i_firstChannelOrderedHadronIds<i_secondChannelOrderedHadronIds)
    {
        return true;
    }
    return false;
}

PartitionFunctionGeneratorOld::PartitionFunctionGeneratorOld(const MCSTHAR::Setup::MicrocanonicalParameterGrid& i_microcanonicalParameterGrid,
                                                             const double i_hadronSamplingTemperature,
                                                             const HadronSamplingGroups& i_hadronSamplingGroups,
                                                             const unsigned long int i_maxChannelSamplingNumber,
                                                             const double i_integrationErrorThreshold,
                                                             const unsigned long int i_maxPhaseSpaceSamplingNumber,
                                                             const double i_phaseSpaceIntegrationErrorThreshold,
                                                             RandomNumberGenerator& io_randomGenerator)
                                                             :m_isPartitionFunctionGeneratorReady(false)
                                                             ,m_isPartitionFunctionAvailable(false)
                                                             ,m_partitionFunctionData(i_microcanonicalParameterGrid)
                                                             ,m_hadronSampling(i_hadronSamplingGroups,
                                                                               i_hadronSamplingTemperature,
                                                                               m_partitionFunctionData.energyDensityValues[0],
                                                                               io_randomGenerator)
                                                             ,m_phaseSpaceIntegrator(i_maxPhaseSpaceSamplingNumber,
                                                                                     i_phaseSpaceIntegrationErrorThreshold,
                                                                                     io_randomGenerator)
                                                             ,m_maxChannelSamplingNumber(i_maxChannelSamplingNumber)
                                                             ,m_integrationSquaredErrorThreshold(i_integrationErrorThreshold*
                                                                                                 i_integrationErrorThreshold)
                                                             ,m_noChannelAvailability(false)
                                                             ,m_previousRunLastEnergyDensityIndex(0)
                                                             ,m_minimumStabilityCheckSamplingNumber(10000)// TODO: Store somehwere with other analogous parameters?
{
    // TODO: add check on input parameter
    m_isPartitionFunctionGeneratorReady = true;
    m_numberOfEnergyDensityValueStr<<m_partitionFunctionData.microcanonicalParameterGrid.numberOfEnergyDensityValues;
}

PartitionFunctionGeneratorOld::~PartitionFunctionGeneratorOld(void)
{
}

void PartitionFunctionGeneratorOld::run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                        const double i_mass)
{
    if(m_isPartitionFunctionGeneratorReady)
    {
        m_isPartitionFunctionAvailable = false;
        
        if(i_mass<=0.)
        {
            // TODO: error code and other checks
            throw HadronizationException("Error in partition function calculation, non positive mass value provided",
                                         __FUNCTION__,811);
        }
        
        // Reset channel phase space data set
        m_channelPhaseSpaceData.clear();
        
        // Loop over energy density grid points
        for(unsigned int energyDensityIndex=0;
            energyDensityIndex<m_partitionFunctionData.microcanonicalParameterGrid.numberOfEnergyDensityValues;
            ++energyDensityIndex)
        {
            // Create log message
            // TODO: store log string as const static (?);
            stringstream energyDensityIndexStr;
            energyDensityIndexStr<<(energyDensityIndex+1);
            // TODO: use function for message pattern creation
            const string message("Running calculation for energy density value number " +
                                 energyDensityIndexStr.str() +
                                 " of " +
                                 m_numberOfEnergyDensityValueStr.str());
            // Log message
            logMessage(message);
            
            // Update energy density
            const double& energyDensityValue(m_partitionFunctionData.energyDensityValues[energyDensityIndex]);
            if(m_previousRunLastEnergyDensityIndex>0)
            {
                try
                {
                    // Update hadron sampling energy density if required (new energy density value)
                    m_hadronSampling.setHadronSamplingEnergyDensity(energyDensityValue);
                }
                catch(HadronizationException& ex)
                {
                    throw ex;
                }
            }

            // Store current energy density index, the corresponding energy density value has been used
            // to set the hadron sampling object and a new set could be required in a subsequent run in case
            // of different energy density value.
            m_previousRunLastEnergyDensityIndex = energyDensityIndex;

            try
            {
                // Run partition function calculation for current energy density value
                runPartitionFunctionCalculation(i_chargeConfiguration,
                                                i_mass,
                                                energyDensityIndex,
                                                energyDensityValue);

                
                if(m_noChannelAvailability)
                {
                    // No hadronization channel available for the considered mass,
                    // reset partition function values to zero and set channel
                    // availability status accordingly
                    m_noChannelAvailability = false;
                    m_isPartitionFunctionAvailable = true;
                    resetPartitionFunctionData();
                    return;
                }
            }
            catch(HadronizationException& ex)
            {
                throw ex;
            }
        }
        
        // Single mass partition function calculation succesfully executed
        m_isPartitionFunctionAvailable = true;
        return;

    }
    else
    {
        // Class not correctly initialized, throw exception
        throw HadronizationException("Error during partition function generation, partition function generator object not correctly initialized",
                                     __FUNCTION__,812);    
    }
}

// TODO: remove the energy density input parameter with the new m_hadronSampling interface
void PartitionFunctionGeneratorOld::runPartitionFunctionCalculation(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                                    const double i_mass,
                                                                    const unsigned int i_energyDensityIndex,
                                                                    const double i_energyDensity)
{
    
    // TODO: optimization, avoid the usage of these vectors
    // Build partition function and strangeness suppression factor storage vectors
    vector<double> strangenessSuppressionFactors(m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues,1.);
    vector<double> weightSum(m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues,0.);
    vector<double> squaredWeightSum(m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues,0.);
    vector<double> phaseSpaceScaledVarianceSum(m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues,0.);
    vector<bool> errorStatus(m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues,false);
    
    // TODO: debug
//    ofstream channelFile;
//    channelFile.open("events.dat");
    
    // Loop over channel sampling
    // TODO: add check on decay channel availability
    // TODO: unsigned int is enough?
    unsigned long int samplingNumber;
    unsigned long int stabilityCheckSamplingNumber(m_minimumStabilityCheckSamplingNumber);
    unsigned long int stabilityCheckSamplingNumberUpdate(m_minimumStabilityCheckSamplingNumber);
    unsigned long int stabilityCheckSamplingNumberIncrement(m_minimumStabilityCheckSamplingNumber);
    unsigned long int stabilityStatusPrintSamplingNumber(100000);
    bool isStabilityReached(false);
    // TODO: check faster pow calculation method
    // TODO: remove unused variable (check error and mean calculations)
    unsigned int hadronSamplingErrorStatus(0);
    double maximumRelativeError(100.);
    unsigned int numberReusedData(0);
    for(samplingNumber=1;samplingNumber<=m_maxChannelSamplingNumber;++samplingNumber)
    {
        try
        {
            
            // Run channel composition sampling
            hadronSamplingErrorStatus = m_hadronSampling.run(i_chargeConfiguration,
                                                             i_mass,
                                                             i_energyDensity);

            // Check sampling error status
            if(hadronSamplingErrorStatus==0)
            {
                // Retrieve sampling weight
                double weigth(m_hadronSampling.getSamplingWeight());
                                
                // Retrieve channel hadron Id
                const unsigned int numberOfHadrons(m_hadronSampling.getNumberOfHadrons());
//                vector<int> hadronIds(numberOfHadrons);// TODO: optimization
//                for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
//                {
//                    hadronIds[hadronIndex] = m_hadronSampling.getHadronSet()[hadronIndex].getIdCode();
//                }
//                
//                // Reorder Id set
//                sort(hadronIds.begin(),hadronIds.end());
                
                // TODO: extend this procedure to hadron sampling weights and strangeness suppression factor,
                // to avoid their recalculation
                // Check if phase space integration data is already present for this channel
//                const map<vector<int>,PhaseSpaceIntegrationData,OrderedHadronIdsComparison>::const_iterator channelDataIt(
//                                                                                            m_channelPhaseSpaceData.find(hadronIds));
//                const map<vector<int>,PhaseSpaceIntegrationData>::const_iterator channelDataIt(m_channelPhaseSpaceData.find(hadronIds));
                double phaseSpaceVariance;
//                if(channelDataIt!=m_channelPhaseSpaceData.end())
//                {
//                    // Channel data already stored, avoid phase space calculation
//
//                    ++numberReusedData;
//                    // Update weight with phase space integral value
//                    weigth *= channelDataIt->second.phaseSpaceIntegral;
//                    
//                    // Retrieve phase space integration variance
//                    // TODO: optimize avoiding error squaring here
//                    phaseSpaceVariance = channelDataIt->second.phaseSpaceIntegralError*
//                                         channelDataIt->second.phaseSpaceIntegralError;
//
//                }
//                else
//                {
                    // Channel data not present, run phase space calculation
                
                    // Compute phase space integral
                    // TODO: avoid the mass copy below
                    vector<double> hadronMasses(numberOfHadrons);
                    //                channelFile<<setprecision(10);
                    //                channelFile<<"samplingNumber "<<samplingNumber<<endl;
                    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
                    {
                        hadronMasses[hadronIndex] = m_hadronSampling.getHadronSet()[hadronIndex].getMass();
                        //                    channelFile<<"Mass "<<hadronMasses[hadronIndex]<<endl;
                    }
                    try
                    {
                        m_phaseSpaceIntegrator.run(i_mass,hadronMasses);
                    }
                    catch(HadronizationException& ex)
                    {
                        throw ex;
                    }
                
                    // Update weight with phase space integral value
                    weigth *= m_phaseSpaceIntegrator.getIntegrationData().phaseSpaceIntegral;

                    // Retrieve phase space integration variance
                    // TODO: optimize avoiding error squaring here
                    phaseSpaceVariance = m_phaseSpaceIntegrator.getIntegrationData().phaseSpaceIntegralError*
                                         m_phaseSpaceIntegrator.getIntegrationData().phaseSpaceIntegralError;
                
//                    // Update channel data set with new channel
//                    m_channelPhaseSpaceData[hadronIds] = m_phaseSpaceIntegrator.getIntegrationData();
//                }
                
                // Compute strangeness suppression factors
                computeStrangenessSuppression(strangenessSuppressionFactors);
                
                // Loop over strangeness suppression parameter grid points and update weight values
                // TODO: avoid this loop by passing a weight vector to computeStrangenessSuppression method
                for(unsigned int gammaSIndex=0;
                    gammaSIndex<m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues;
                    ++gammaSIndex)
                {
                    // Retrieve strangeness suppression
                    const double& strangenessSuppressionFactor(strangenessSuppressionFactors[gammaSIndex]);
                    
                    // Update integration weigth with strangeness suppression
                    const double suppressedWeight(weigth*strangenessSuppressionFactor);
                    
                    // Update weight, squared weight and phase space scaled error sum
                    double& weightSumElement(weightSum[gammaSIndex]);
                    double& squaredWeightSumElement(squaredWeightSum[gammaSIndex]);
                    double& phaseSpaceScaledVarianceSumElement(phaseSpaceScaledVarianceSum[gammaSIndex]);
                    weightSumElement += suppressedWeight;
                    const double squaredWeight(suppressedWeight*suppressedWeight);
                    squaredWeightSumElement += squaredWeight;
                    phaseSpaceScaledVarianceSumElement += phaseSpaceVariance*squaredWeight;
                    
                    // Check calculation stability
                    if(samplingNumber==stabilityCheckSamplingNumber)
                    {
                        // Check error status for current strangeness suppression parameter value
                        // TODO: can the usage of squared error and squared sum be avoided?
                        // TODO: add numerical check for /0, etc...
                        // TODO: optimize
                        // TODO: add comments on calculations
                        if(weightSumElement>0.)
                        {
                            const double oneOnWeightSumElementSquared(1./(weightSumElement*weightSumElement));
                            double channelSumSquaredRelativeError((samplingNumber*squaredWeightSumElement*oneOnWeightSumElementSquared - 1.)/
                                                                  (samplingNumber-1));
                                                    
                            // Check for numerical relative squared error non positiveness
                            if(std::abs(channelSumSquaredRelativeError)<=partitionFunctionSquaredErrorZeroRange)
                            {
                                // Relative squared error is approximatively zero, set it to exactly zero
                                channelSumSquaredRelativeError = 0.;
                            }
                            else if(channelSumSquaredRelativeError<-partitionFunctionSquaredErrorZeroRange)
                            {
                                // Relative squared error is negative, throw exception
                                throw HadronizationException("Error during channel sampling: negative channel sum relative squared error",
                                                             __FUNCTION__,815);
                            }
                            
                            // Compute integration squared relative error for curren strangeness suppression parameter value
                            const double squaredRelativeError = channelSumSquaredRelativeError +
                            phaseSpaceScaledVarianceSumElement*oneOnWeightSumElementSquared/samplingNumber;
                            
                            // Update maximum relative error
                            if(squaredRelativeError>maximumRelativeError)
                            {
                                maximumRelativeError = squaredRelativeError;
                            }
                            
    //                        channelFile<<samplingNumber<<" "
    //                                   <<m_phaseSpaceIntegrator.getIntegrationData().phaseSpaceIntegral<<" "
    //                                   <<m_phaseSpaceIntegrator.getIntegrationData().phaseSpaceIntegralError<<" "
    //                                    <<m_hadronSampling.getSamplingWeight()<<" "
    //                                    <<weightSumElement/samplingNumber<<" "
    //                                    <<squaredRelativeError*weightSumElement*weightSumElement<<endl;
    //                        for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    //                        {
    //                            channelFile<<m_hadronSampling.getHadronSet()[hadronIndex].getIdCode()<<" ";
    //                        }
    //                        channelFile<<endl;
                            
                            // TODO: check on relative or absolute error?
                            // TODO: check denominator
                            // Check stability status
                            if(squaredRelativeError<m_integrationSquaredErrorThreshold)
                            {
                                // Stability reached for current strangeness suppression parameter value
                                errorStatus[gammaSIndex] = true;
                            }
                        }
                    }
                    
                }
            }
            else
            {
                // Error during hadron channel sampling
                if(hadronSamplingErrorStatus == 3)
                {
                    // Throw exception, unvalid cluster charge configuration
                    string errorMsg("Error during channel sampling: double heavy flavored cluster provided.");
                    throw HadronizationException(errorMsg,__FUNCTION__,813);
                }
                else
                {
                    // Apparently no hadronization channels exist for the provided cluster mass
                    // and charge configuration and for the available hadron set: set partition function
                    // to zero and log warning
                    string warningMsg("No hadronization channel found. Partition function values are set to zero for current cluster mass value. Check hadron list file for channel availability");
                    logMessage(warningMsg,WARNING);
                    m_noChannelAvailability = true;
                    return;
                }
            }
        }
        catch(HadronizationException& ex)
        {            
            // TODO: factorize into a function? Add save status function in exception handler
            stringstream massValueStr;
            stringstream energyDensityIndexStr;
            stringstream samplingNumberStr;
            massValueStr<<i_mass;
            energyDensityIndexStr<<(i_energyDensityIndex+1);
            samplingNumberStr<<samplingNumber;

            string errorMsg(ex.getErrorMessage());
            errorMsg += "\n\nMass value " +
            massValueStr.str() +
            " GeV \nEnergy density value " +
            energyDensityIndexStr.str() +
            "\nIntegration sampling number " +
            samplingNumberStr.str();
            throw HadronizationException(errorMsg,__FUNCTION__,ex.getReturnValue());
        }
        
        // Check error status
        if(samplingNumber==stabilityCheckSamplingNumber)
        {
            if(samplingNumber==stabilityStatusPrintSamplingNumber)
            {
                // Create log message
                // TODO: store log string as const static (?);
                stringstream samplingNumberStr;
                stringstream maximumRelativeErrorValueStr;
                samplingNumberStr<<samplingNumber;
                maximumRelativeErrorValueStr<<sqrt(maximumRelativeError);
                // TODO: use function for message pattern creation
                const string message("Maximum relative error " +
                                     maximumRelativeErrorValueStr.str() +
                                     " after " +
                                     samplingNumberStr.str() +
                                     " hadronization channel samplings.");
                // Log message
                logMessage(message);
                
                
                stabilityStatusPrintSamplingNumber += 100000;

            }
            maximumRelativeError = -1;
            
            if(find(errorStatus.begin(),errorStatus.end(),false)==errorStatus.end())
            {
                // Stability reached for all strangeness suppression parameter values, stop integration loop
                isStabilityReached = true;
                break;
            }
            else
            {
                // TODO: debug
//                stabilityCheckSamplingNumber +=m_minimumStabilityCheckSamplingNumber;// TODO: if this is final strategy remove code below
                stabilityCheckSamplingNumber +=1;// TODO: if this is final strategy remove code below
//                // TODO: check if better conter solution exists
//                // Compute next sampling number value for stability check
//                if(samplingNumber==stabilityCheckSamplingNumberUpdate)
//                {
//                    // TODO: check unsigned (everywhere)
//                    stabilityCheckSamplingNumberUpdate = (samplingNumber*10);
//                    stabilityCheckSamplingNumberIncrement = stabilityCheckSamplingNumberUpdate/2;
//                    stabilityCheckSamplingNumber = stabilityCheckSamplingNumberIncrement;
//                }
//                else
//                {
//                    stabilityCheckSamplingNumber += stabilityCheckSamplingNumberIncrement;
//                }
//                if(stabilityCheckSamplingNumber>m_maxChannelSamplingNumber)
//                {
//                    stabilityCheckSamplingNumber = m_maxChannelSamplingNumber;
//                }
            }
        }
    }
    
    //    channelFile.close();
    
    // Integration procedure finished, set partition function data for output
    if(isStabilityReached==false)
    {
        --samplingNumber;
    }
    
    // Compute start position index of data in global storage structure for current energy density value
    unsigned int dataGlobalIndex(i_energyDensityIndex*m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues);
    
    // Loop over strangeness suppression parameter grid points and set partition function data
    for(unsigned int gammaSIndex=0;
        gammaSIndex<m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues;
        ++gammaSIndex,++dataGlobalIndex)
    {
        // TODO: optimize calculation below
        // TODO: check unsigned int
        const double& weightSumElement(weightSum[gammaSIndex]);
        
        // Compute partition function valie
        m_partitionFunctionData.partitionFunctions[dataGlobalIndex] = weightSumElement/samplingNumber;
        
        // Check for numerical squared error non positiveness
        double channelSumSquaredError((squaredWeightSum[gammaSIndex]-weightSumElement*weightSumElement/samplingNumber)/(samplingNumber-1));
        if(std::abs(channelSumSquaredError)<=partitionFunctionSquaredErrorZeroRange)
        {
            // Squared error is approximatively zero, set it to exactly zero
            channelSumSquaredError = 0.;
        }
        else if(channelSumSquaredError<-partitionFunctionSquaredErrorZeroRange)
        {
            // Squared error is negative, throw exception
            throw HadronizationException("Error during channel sampling: negative channel sum squared error",
                                         __FUNCTION__,816);
        }

        // Compute partition function error
        m_partitionFunctionData.partitionFunctionErrors[dataGlobalIndex] =
        sqrt((channelSumSquaredError + phaseSpaceScaledVarianceSum[gammaSIndex]/(samplingNumber*samplingNumber))/samplingNumber);
        m_partitionFunctionData.partitionFunctionErrorStatus[dataGlobalIndex] = errorStatus[gammaSIndex];
        m_partitionFunctionData.numberOfChannelSamplings[dataGlobalIndex] = samplingNumber;
        
    }
    
//    cout<<" numberReusedData "<<numberReusedData<<endl;
//    cout<<" numberOfChannels "<<m_channelPhaseSpaceData.size()<<endl;
    
    return;
}

const PartitionFunctionData& PartitionFunctionGeneratorOld::getPartitionFunctionData(void) const
{
    if(m_isPartitionFunctionAvailable)
    {
        return m_partitionFunctionData;
    }
    else
    {
        throw HadronizationException("Error during partition function data retrieval, partition function data not available",
                                     __FUNCTION__,814);
    }
}

// TODO: remove after debugging/testing
void PartitionFunctionGeneratorOld::computeStrangenessSuppression(vector<double>& io_strangenessSuppressionFactors)
{
    // TODO: optimize calculation and memory management (e.g., can't we directly update the channel weights?)
    io_strangenessSuppressionFactors = vector<double>(m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues,1.);
    
    // Loop over current channel hadrons
    const vector<HadronData>::const_iterator hadronSetEnd(m_hadronSampling.getHadronSet().end());
    for(vector<HadronData>::const_iterator hadronSetIt = m_hadronSampling.getHadronSet().begin();
        hadronSetIt!=hadronSetEnd;++hadronSetIt)
    {
        // Retrieve single hadron ssbar wave function component
        const double ssbarComponentValue(hadronSetIt->getSSBARComponentValue());
        // Check hadron strange composition
        if(ssbarComponentValue>0.)
        {
            // Open strange hadron // TODO chedk "open strange"
            // Loop over strangeness suppression parameter grid points and update strangeness suppression contribution
            for(unsigned int gammaSIndex=0;
                gammaSIndex<m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues;
                ++gammaSIndex)
            {
                // TODO: optimization
                const double& gammaSValue(m_partitionFunctionData.gammaSValues[gammaSIndex]);
                io_strangenessSuppressionFactors[gammaSIndex] *= (hadronSetIt->getOneMinusSSBARComponentValue() +
                                                                  ssbarComponentValue*gammaSValue*gammaSValue);
            }
        }
        else
        {
            // Possible ssbar hadronic status
            // Retrieve number of strange quarks
            const unsigned int& numberOfStrangeQuarks(hadronSetIt->getNumbeOfStrangeQuarks());
            
            if(numberOfStrangeQuarks>0)
            {
                // Loop over strangeness suppression parameter grid points
                for(unsigned int gammaSIndex=0;
                    gammaSIndex<m_partitionFunctionData.microcanonicalParameterGrid.numberOfGammaSValues;
                    ++gammaSIndex)
                {
                    // TODO: why is the following cast required?
                    const double& gammaSValue(m_partitionFunctionData.gammaSValues[gammaSIndex]);
                    io_strangenessSuppressionFactors[gammaSIndex] *=
                    std::pow(gammaSValue,
                             static_cast<int>(numberOfStrangeQuarks));
                }
            }
        }
    }
    
    return;
}

void PartitionFunctionGeneratorOld::resetPartitionFunctionData(void)
{
    // TODO: factorize with standard partition function data assignement
    // TODO: zero reassignement in case of subsequent empty channel sets
    
    // Set to zero the partition function values and the corresponding error and number of sampled channels. Set to true the
    // their stability status, assuming that no hadronization channel is available for the considered mass, charge configuration
    // and hadron set
    m_partitionFunctionData.partitionFunctions.assign(m_partitionFunctionData.microcanonicalParameterGrid.totalNumberOfGridPoints,0.);
    m_partitionFunctionData.partitionFunctionErrors.assign(m_partitionFunctionData.microcanonicalParameterGrid.totalNumberOfGridPoints,0.);
    m_partitionFunctionData.partitionFunctionErrorStatus.assign(m_partitionFunctionData.microcanonicalParameterGrid.totalNumberOfGridPoints,true);
    m_partitionFunctionData.numberOfChannelSamplings.assign(m_partitionFunctionData.microcanonicalParameterGrid.totalNumberOfGridPoints,0);

    return;
}
