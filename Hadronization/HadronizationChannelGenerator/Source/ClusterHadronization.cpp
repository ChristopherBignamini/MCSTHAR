#include "../Include/ClusterHadronization.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include "../../HadronizationHadronSet/Include/HadronData.h"
#include <vector>

ClusterHadronization::ClusterHadronization(const double i_energyDensity,
                                           const double i_gammaS,
                                           const double i_hadronSamplingTemperature,
                                           const double i_hadronSamplingEnergyDensity,
                                           const string& i_partitionFunctionDataPath,
                                           const HadronSamplingGroups& i_hadronSamplingGroups,
                                           RandomNumberGenerator& io_randomGenerator)
                                          :m_isHadronizationChannelAvailable(false)
                                          ,m_isClusterHadronizationReady(false)
                                          ,m_hadronSampling(i_hadronSamplingGroups,
                                                            i_hadronSamplingTemperature,
                                                            i_hadronSamplingEnergyDensity,
                                                            io_randomGenerator)
                                          ,m_phaseSpaceSampling(io_randomGenerator)
                                          ,m_gammaS(i_gammaS)
                                          ,m_partitionFunction(i_partitionFunctionDataPath,
                                                               i_energyDensity,
                                                               i_gammaS,
                                                               true)
{
    
    if(m_gammaS<0. || m_gammaS>1.)
    {
        throw HadronizationException("Error during cluster hadronization object creation, strangeness suppression parameter outside [0,1] range provided",__FUNCTION__,101);
    }
    
    m_isClusterHadronizationReady = true;
}

ClusterHadronization::~ClusterHadronization(void)
{
}

unsigned int ClusterHadronization::runHadronization(const Cluster& i_cluster)
{    
    if(m_isClusterHadronizationReady)
    {
        m_isHadronizationChannelAvailable = false;
        
        // TODO: add input cluster check, maybe only in assert (mass!=0, etc)
        try
        {
            // Partition function calculation error status
            unsigned int partitionFunctionErrorStatus;

            
            // Retrieve partition function value
            const MCSTHAR::Utilities::ChargeConfiguration clusterCharges(i_cluster.getStrangeCharge(),
                                                                         i_cluster.getCharmCharge(),
                                                                         i_cluster.getBottomCharge(),
                                                                         i_cluster.getElectricCharge(),
                                                                         i_cluster.getBaryonicCharge());
            
            const double partitionFunctionValue(m_partitionFunction.
                                                getPartitionFunctionValue(clusterCharges,
                                                                          i_cluster.getMass(),
                                                                          partitionFunctionErrorStatus,
                                                                          true));
            
            // Check hadron sampling execution error status
            if(partitionFunctionErrorStatus!=0)
            {
                return partitionFunctionErrorStatus;
            }
            
            // Check partition function value
            // TODO: add tolerance?
            if(partitionFunctionValue>0.0)
            {                
                // Run channel composition sampling
                const unsigned int hadronSamplingErrorStatus(m_hadronSampling.run(clusterCharges,
                                                                                  i_cluster.getMass(),
                                                                                  i_cluster.getEnergyDensity()));

                if(hadronSamplingErrorStatus==0)
                {
                    // Run phase space configuration sampling
                    m_phaseSpaceSampling.run(i_cluster,m_hadronSampling.getHadronSet());
                    
                    // Create cluster hadronization event
                    buildHadronizationEvent(i_cluster.getIndex(),partitionFunctionValue);

                }
                else
                {
                    // Check hadron sampling execution error status
                    if(hadronSamplingErrorStatus==1)
                    {
                        return 4;// TODO: switch to exit code static table/enumerator (also for exception)
                    }
                    if(hadronSamplingErrorStatus==2)
                    {
                        return 5;// TODO: switch to exit code static table/enumerator (also for exception)
                    }
                    if(hadronSamplingErrorStatus==3)
                    {
                        return 6;// TODO: switch to exit code static table/enumerator (also for exception)
                    }
                }                
            }
            else
            {
                // Null partition function value, return non zero error status
                return 3;// TODO: switch to exit code static table/enumerator (also for exception)
            }
            
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
        
        m_isHadronizationChannelAvailable = true;
        
        // Single cluster hadronization succesfully executed
        return 0;

    }
    else
    {
        // Class not correctly initialized, throw exception
        throw HadronizationException("Error during cluster hadronization, cluster hadronization object not correctly initialized",
                                     __FUNCTION__,102);
    }
}

const HadronizationChannel& ClusterHadronization::getHadronizationChannel(void) const
{
    if(m_isHadronizationChannelAvailable)
    {
        return m_hadronizationChannel;
    }
    else
    {
        throw HadronizationException("Error during hadronization channel retrieval, hadronization channel not available",
                                     __FUNCTION__,103);
    }
}

void ClusterHadronization::buildHadronizationEvent(const unsigned int i_clusterIndex,const double i_partitionFunctionValue)
{
    try
    {
        // Set hadronization channel cluster index
        m_hadronizationChannel.clusterIndex = i_clusterIndex;

        // Build final state particles and compute strangeness suppression weight
        double strangenessSuppressionFactor = 1.;
        m_hadronizationChannel.particles.resize(m_hadronSampling.getNumberOfHadrons());
        unsigned int hadronIndex = 0;
        const vector<HadronData>::const_iterator hadronSetEnd(m_hadronSampling.getHadronSet().end());
        vector<TLorentzVector>::const_iterator hadronMomentaIt(m_phaseSpaceSampling.getPhaseSpace().begin());
        // TODO: optimize the for loop (maybe using copy)
        for(vector<HadronData>::const_iterator hadronSetIt = m_hadronSampling.getHadronSet().begin();
            hadronSetIt!=hadronSetEnd;++hadronSetIt,++hadronMomentaIt,++hadronIndex)
        {
            const double ssbarComponentValue(hadronSetIt->getSSBARComponentValue());
            if(ssbarComponentValue>0.)
            {
                strangenessSuppressionFactor *= (hadronSetIt->getOneMinusSSBARComponentValue() +
                                                 ssbarComponentValue*m_gammaS*m_gammaS);
            }
            else
            {
                strangenessSuppressionFactor *= pow(m_gammaS,hadronSetIt->getNumbeOfStrangeQuarks());
            }
            
            
            m_hadronizationChannel.particles[hadronIndex].setDataAndMomentum(*hadronSetIt,*hadronMomentaIt);
            m_hadronizationChannel.particles[hadronIndex].setIndex(hadronIndex);
            m_hadronizationChannel.particles[hadronIndex].setParentIndex(i_clusterIndex);
        }
        
        // Set hadronization channel weight
        m_hadronizationChannel.weight = m_hadronSampling.getSamplingWeight()*
                                        m_phaseSpaceSampling.getSamplingWeight()*
                                        strangenessSuppressionFactor/i_partitionFunctionValue;        
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

