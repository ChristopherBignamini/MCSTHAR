#include "../Include/HadronizationHandler.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/ROOTRandomNumberGenerator.h"

// TODO: check in all classes the exception chain behavior
// TODO: Contructor must set object in a meaningful status... Avoid the m_clusterHadronization.set calls below
//       with a better constructor!
// TODO: build hadronization setup structure with the required parameters only

HadronizationHandler::HadronizationHandler(const HadronizationSetup& i_hadronizationSetup)
                                          :m_runHadronization(true)
                                          ,m_hadronizationSetup(i_hadronizationSetup)
                                          ,m_randomGenerator(new ROOTRandomNumberGenerator(m_hadronizationSetup.randomNumberGeneratorSeed))
                                          ,m_hadronSamplingGroups(m_hadronizationSetup.hadronDataSetFileName,
                                                                  m_hadronizationSetup.lightHadronMaxMass)
                                          ,m_clusterHadronization(m_hadronizationSetup.energyDensity,
                                                                  m_hadronizationSetup.gammaS,
                                                                  m_hadronizationSetup.samplingTemperature,
                                                                  m_hadronizationSetup.samplingEnergyDensity,
                                                                  m_hadronizationSetup.partitionFunctionDataSetPath,
                                                                  m_hadronSamplingGroups,
                                                                  *m_randomGenerator)
                                          ,m_clusterMerging(NULL)
{
    try
    {
        // If needed build pointer to ClusterMerging
        if(m_hadronizationSetup.clusterMergingFlag)
        {
                m_clusterMerging = new ClusterMerging(m_hadronizationSetup.clusterMergingMinMass,
                                                      m_hadronizationSetup.charmClusterMergingMinMass,
                                                      m_hadronizationSetup.bottomClusterMergingMinMass,
                                                      m_hadronizationSetup.maxClusterMass);
        }
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

// TODO: fix random number generator creation
HadronizationHandler::HadronizationHandler(const HadronizationHandler& i_hadronizationHandler)
                                          :m_runHadronization(i_hadronizationHandler.m_runHadronization)
                                          ,m_hadronizationSetup(i_hadronizationHandler.m_hadronizationSetup)
                                          ,m_randomGenerator(new ROOTRandomNumberGenerator(
                                                    *(static_cast<ROOTRandomNumberGenerator*>(i_hadronizationHandler.m_randomGenerator))))
                                          ,m_hadronSamplingGroups(i_hadronizationHandler.m_hadronSamplingGroups)
                                          ,m_eventRecord(i_hadronizationHandler.m_eventRecord)
                                          ,m_clusterHadronization(i_hadronizationHandler.m_clusterHadronization)
                                          ,m_clusterMerging(NULL)
{    
    // If needed build pointer to ClusterMerging
    if(m_hadronizationSetup.clusterMergingFlag)
    {
        m_clusterMerging = new ClusterMerging(*i_hadronizationHandler.m_clusterMerging);
    }
}

HadronizationHandler::~HadronizationHandler(void)
{
    if(m_clusterMerging)
    {
        delete m_clusterMerging;
        m_clusterMerging = NULL;
    }
    
    delete m_randomGenerator;
}

void HadronizationHandler::loadClusters(const vector<Cluster>& i_clusters)
{
    try
    {
        // Import cluster in event record
        m_eventRecord.setInputClusters(i_clusters);
        // Run hadronization on new set of clusters
        m_runHadronization = true;
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}
    
unsigned int HadronizationHandler::runHadronization(void)
{
    try
    {        
        if(m_runHadronization)
        {
            // Check if cluster merging is required
            if(m_hadronizationSetup.clusterMergingFlag)
            {
                // Run merging
                m_clusterMerging->runMerging(m_eventRecord.getEventClusters());
                
                // Update event record with clusters created by merging procedure
                if(m_clusterMerging->getCreatedClusterNumber()>0)
                {                    
                    m_eventRecord.addClustersFromMerging(m_clusterMerging->getCreatedClusters());
                }
            }
            
            // Check availability of hadronizable clusters
            if(m_eventRecord.getNumberOfHadronizableClusters()>0)
            {
                // Loop over hadronizable clusters stored in event record
                const vector<Cluster> hadronizableClusters(m_eventRecord.getHadronizableClusters());
                for(vector<Cluster>::const_iterator clusterIt = hadronizableClusters.begin();
                    clusterIt!= hadronizableClusters.end();++clusterIt)
                {
                    // Run cluster hadronization
                    const unsigned int hadronizationErrorStatus(m_clusterHadronization.runHadronization(*clusterIt));
                    
                    // Single cluster hadronization succesfully executed
                    if(hadronizationErrorStatus==0)
                    {
                        // Update event record
                        m_eventRecord.addHadronizationChannel(m_clusterHadronization.getHadronizationChannel());
                        
                    }
                    else
                    {
                        // Break hadronization loop and return error code
                        return hadronizationErrorStatus;
                    }
                }
                
                // Close event record
                m_eventRecord.closeEventRecord();
                                
            }
            else
            {
                // No cluster available for hadronization, throw exception
                throw HadronizationException("No cluster available for hadronization",
                                             __FUNCTION__,512);
            }
            
            m_runHadronization = false;
        }
        
        // Cluster set hadronization succesfully executed
        return 0;
    }
    catch(HadronizationException ex)
    {
        throw ex;
    }
}

const HadronizationEventRecord& HadronizationHandler::getEventRecord(void) const
{
    if(m_runHadronization)
    {
        // Hadronization not performed
        throw HadronizationException("Error during hadronization event record retrieval, hadronization not performed yet",
                                     __FUNCTION__,511);
    
    }
    return m_eventRecord;
}

HadronizationHandler& HadronizationHandler::operator=(const HadronizationHandler& i_hadronizationHandler)
{
    m_runHadronization = i_hadronizationHandler.m_runHadronization;
    m_hadronizationSetup = i_hadronizationHandler.m_hadronizationSetup;
    *m_randomGenerator = *(i_hadronizationHandler.m_randomGenerator);
    m_hadronSamplingGroups = i_hadronizationHandler.m_hadronSamplingGroups;
    m_eventRecord = i_hadronizationHandler.m_eventRecord;
    m_clusterHadronization = i_hadronizationHandler.m_clusterHadronization;
    m_clusterMerging = NULL;
    
    
    // If needed build pointer to ClusterMerging
    if(m_hadronizationSetup.clusterMergingFlag)
    {
        m_clusterMerging = new ClusterMerging(*i_hadronizationHandler.m_clusterMerging);
    }
    
    return *this;
}
