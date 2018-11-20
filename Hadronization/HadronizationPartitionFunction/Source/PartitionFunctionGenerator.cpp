#include "../Include/PartitionFunctionGenerator.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/buildParameterRegularGrid.h"
#include <cmath>
#include <cassert>

#ifdef _USEROOTRANDOM
    #include "../../../Utilities/Include/ROOTRandomNumberGenerator.h"
#else
    #include "../../../Utilities/Include/GSLRandomNumberGenerator.h"// TODO: use preprocessor directives
#endif

// TODO: check everywhere exit code
// TODO: check everywhere input parameter constness

PartitionFunctionGenerator::PartitionFunctionGenerator(const double i_minGammaS,
                                                       const double i_maxGammaS,
                                                       const unsigned int i_numberOfGammaSValues,
                                                       const HadronSamplingGroups& i_hadronSamplingGroups,
                                                       const unsigned int i_randomNumberGeneratorSeed,
                                                       const MPI_Comm& i_mpiCommunicator)
                                                      :m_isPartitionFunctionGeneratorReady(false)
                                                      ,m_gammaSValues(i_numberOfGammaSValues,0.)
                                                      ,m_numberOfGammaSValues(i_numberOfGammaSValues)
                                                      ,m_localStabilityCheckSamplingNumber(0)
                                                      ,m_localStabilityCheckSamplingNumberStart(0)
                                                      #ifdef _USEROOTRANDOM // TODO: move to final code (remove root)
                                                        ,m_randomGenerator(new ROOTRandomNumberGenerator(i_randomNumberGeneratorSeed))
                                                      #else
                                                        ,m_randomGenerator(new GSLRandomNumberGenerator(i_randomNumberGeneratorSeed))
                                                      #endif
#if 0
                                                      ,m_hadronSampling(i_hadronSamplingGroups,
                                                                        0.160,
                                                                        0.2,// TODO: to be removed with new interface
                                                                        *m_randomGenerator)
#endif
                                                      ,m_hadronSamplingNew(i_hadronSamplingGroups,
                                                                           *m_randomGenerator)
                                                      ,m_phaseSpaceSampling(*m_randomGenerator)
                                                      ,m_weightSumData(NULL)
                                                      ,m_weightSumDataSize(2*i_numberOfGammaSValues+1)
                                                      ,m_strangenessSuppressionFactors(i_numberOfGammaSValues,1.)
                                                      ,m_mpiCommunicator(i_mpiCommunicator) 
{    
    // Build strangeness suppression parameter value grid
    try
    {
        buildParameterRegularGrid(i_minGammaS,
                                  i_maxGammaS,
                                  i_numberOfGammaSValues,
                                  m_gammaSValues);
    }
    catch (HadronizationException& ex)
    {
        throw ex;
    }
    
    // Allocate memory for partition function calculation
    m_weightSumData = new double[m_weightSumDataSize];
}

PartitionFunctionGenerator::~PartitionFunctionGenerator(void)
{
    if(m_weightSumData!=NULL)
    {
        delete [] m_weightSumData;
        m_weightSumData = NULL;
    }
    if(m_randomGenerator!=NULL)
    {
        delete m_randomGenerator;
    }
}

void PartitionFunctionGenerator::run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                     const double i_mass,
                                     const double i_energyDensity)
{
    if(m_isPartitionFunctionGeneratorReady)
    {
        if(i_mass<=0.)
        {
            // TODO: error code and other checks
            throw HadronizationException("Error in partition function calculation, non positive mass value provided",
                                         __FUNCTION__,811);
        }
        
        if(i_energyDensity<=0.)
        {
            // TODO: error code and other checks
            throw HadronizationException("Error in partition function calculation, non positive energy density value provided",
                                         __FUNCTION__,811);
        }
        
        // Update energy density (TODO: the update should be performed only in case of new energy density values)
        try
        {
            // TODO: to be removed with new hadron sampling interface (sampling density = physical density)
//            m_hadronSampling.setHadronSamplingEnergyDensity(i_energyDensity);
        }
        catch(HadronizationException& ex)
        {
            // TODO: add exception handling for MPI
            throw ex;
        }
        
        try
        {            
            // Reset stability check sampling number value
            m_localStabilityCheckSamplingNumber = m_localStabilityCheckSamplingNumberStart;// TODO: design error??

//            // TEST
//            m_hadronSampling.updateRapidities();
//            // TEST
            
            // Run partition function calculation for current energy density value
            runPartitionFunctionCalculation(i_chargeConfiguration,
                                            i_mass,
                                            i_energyDensity);
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
        
        return;
        
    }
    else
    {
        // Class not correctly initialized, throw exception
        throw HadronizationException("Error during partition function generation, partition function generator object not correctly initialized",
                                     __FUNCTION__,812);
    }
}

// TODO: this method can be factorized with the code used during event generation
void PartitionFunctionGenerator::computeStrangenessSuppression(const vector<double>& i_ssbarComponentValues,
                                                               const vector<double>& i_oneMinusSSBARComponentValues,
                                                               const vector<unsigned int>& i_numbeOfStrangeQuarks)
{
    // TODO: after final design solution, change into error check
    assert(i_ssbarComponentValues.size() == i_oneMinusSSBARComponentValues.size());
    assert(i_ssbarComponentValues.size() == i_numbeOfStrangeQuarks.size());
    
    // TODO: optimize calculation and memory management (e.g., can't we directly update the channel weights?)
    fill(m_strangenessSuppressionFactors.begin(),m_strangenessSuppressionFactors.end(),1.);
    
    // Loop over current channel hadrons
    for(unsigned int particleIndex=0;particleIndex<i_ssbarComponentValues.size();++particleIndex)
    {
        // Retrieve single hadron ssbar wave function component
        const double& ssbarComponentValue(i_ssbarComponentValues.at(particleIndex));
        // Check hadron strange composition
        if(ssbarComponentValue>0.)
        {
            // Open strange hadron // TODO chedk "open strange"
            // Loop over strangeness suppression parameter grid points and update strangeness suppression contribution
            const double& oneMinusSSBARComponentValue(i_oneMinusSSBARComponentValues.at(particleIndex));
            for(unsigned int gammaSIndex=0;gammaSIndex<m_numberOfGammaSValues;++gammaSIndex)
            {
                // TODO: optimization
                const double& gammaSValue(m_gammaSValues.at(gammaSIndex));
                m_strangenessSuppressionFactors.at(gammaSIndex) *= (oneMinusSSBARComponentValue +
                                                                    ssbarComponentValue*gammaSValue*gammaSValue);
            }
        }
        else
        {
            // Possible ssbar hadronic status
            // Retrieve number of strange quarks
            const unsigned int& numberOfStrangeQuarks(i_numbeOfStrangeQuarks.at(particleIndex));
            
            if(numberOfStrangeQuarks>0)
            {
                // Loop over strangeness suppression parameter grid points
                for(unsigned int gammaSIndex=0;gammaSIndex<m_numberOfGammaSValues;++gammaSIndex)
                {
                    // TODO: why is the following cast required?
                    const double& gammaSValue(m_gammaSValues.at(gammaSIndex));
                    m_strangenessSuppressionFactors.at(gammaSIndex) *=
                    std::pow(gammaSValue,static_cast<int>(numberOfStrangeQuarks));
                }
            }
        }
    }
    
    return;
}

void PartitionFunctionGenerator::runNewCalculationIterations(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                                             const double i_mass,
                                                             const double i_energyDensity)
{
    
    // Reset local calcuation data
    fill(m_weightSumData,m_weightSumData+m_weightSumDataSize,0.);
    
    // Loop over channel sampling
    // TODO: add check on decay channel availability
    // TODO: unsigned int is enough?
    unsigned int hadronSamplingErrorStatus(0);
    for(unsigned long long int localSamplingNumber = 0;localSamplingNumber<m_localStabilityCheckSamplingNumber;
        ++localSamplingNumber)
    {
        try
        {
 
            // Run channel composition sampling
//            hadronSamplingErrorStatus = m_hadronSampling.run(i_chargeConfiguration,
//                                                             i_mass,
//                                                             i_energyDensity);
            hadronSamplingErrorStatus = m_hadronSamplingNew.run(i_chargeConfiguration,
                                                             i_mass,
                                                             i_energyDensity);
            
            // Check sampling error status
            if(hadronSamplingErrorStatus==0)
            {
                // Retrieve sampling weight
//                double weigth(m_hadronSampling.getSamplingWeight());
                double weigth(m_hadronSamplingNew.getSamplingWeight());
                
                // Retrieve number of hadrons
//                const unsigned int numberOfHadrons(m_hadronSampling.getNumberOfHadrons());
                const unsigned int numberOfHadrons(m_hadronSamplingNew.getNumberOfHadrons());
                
                // Retrieve hadron data (mass and strangeness composition)
                // TODO: avoid the mass copy below
                vector<double> hadronMasses(numberOfHadrons);
                vector<double> ssbarComponentValues(numberOfHadrons);
                vector<double> oneMinusSSBARComponentValues(numberOfHadrons);
                vector<unsigned int> numbeOfStrangeQuarks(numberOfHadrons);
                for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
                {
//                    const HadronData& currentHadronData(m_hadronSampling.getHadronSet().at(hadronIndex));
                    const HadronData& currentHadronData(m_hadronSamplingNew.getHadronSet().at(hadronIndex));
                    hadronMasses.at(hadronIndex) = currentHadronData.getMass();// TODO: remove .at everywhere
                    ssbarComponentValues.at(hadronIndex) = currentHadronData.getSSBARComponentValue();
                    oneMinusSSBARComponentValues.at(hadronIndex) = currentHadronData.getOneMinusSSBARComponentValue();
                    numbeOfStrangeQuarks.at(hadronIndex) = currentHadronData.getNumbeOfStrangeQuarks();
                }
                
                // Run phase space configuration sampling and update weight
                try
                {
                    weigth *= m_phaseSpaceSampling.run(i_mass,hadronMasses);
                }
                catch(HadronizationException& ex)
                {
                    throw ex;
                }
                
                // Compute strangeness suppression factors
                // TODO: find final interface (ssbarComponentValues and oneMinusSSBARComponentValues are not indepentent)
                computeStrangenessSuppression(ssbarComponentValues,
                                              oneMinusSSBARComponentValues,
                                              numbeOfStrangeQuarks);
                
                // Loop over strangeness suppression parameter grid points and update weight values
                // TODO: avoid this loop by passing a weight vector to computeStrangenessSuppression method
                unsigned int dataIndex(0);
                for(unsigned int gammaSIndex=0;gammaSIndex<m_numberOfGammaSValues;
                    ++gammaSIndex,dataIndex+=2)
                {
                    // Retrieve strangeness suppression
                    // Update integration weigth with strangeness suppression
                    const double suppressedWeight(weigth*m_strangenessSuppressionFactors.at(gammaSIndex));
                    
                    // TODO: Non positive weights???!
                    // Update weight and squared weight sum
                    m_weightSumData[dataIndex] += suppressedWeight;
                    m_weightSumData[dataIndex+1] += suppressedWeight*suppressedWeight;
                }
                
            }
            else
            {
                // Error during hadron channel sampling
                if(hadronSamplingErrorStatus == 3)
                {
                    // Throw exception, not valid cluster charge configuration
                    throw HadronizationException("Error during channel sampling: double heavy flavored cluster provided.",
                                                 __FUNCTION__,
                                                 813);
                }
                else
                {
                    // Apparently no hadronization channels exist for the provided cluster mass,
                    // reset data and break
                    fill(m_weightSumData,
                         m_weightSumData+m_weightSumDataSize,
                         0.);
                    m_weightSumData[m_weightSumDataSize-1] = 1.;
                    
                    break;
                }
            }
            
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
    }
    
    // TODO: should we check for non zero weight sum in each task?
    // Master task receives local calculation data from slaves and makes reduction on global data
    try
    {
        mergeLocalSums();
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
    return;
}

