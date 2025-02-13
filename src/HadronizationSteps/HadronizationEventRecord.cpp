#include "../Include/HadronizationEventRecord.h"
#include "../../../Utilities/Include/HadronizationException.h"

HadronizationEventRecord::HadronizationEventRecord(void)
                                                  :m_hadronizableClusterNumber(0)
                                                  ,m_weight(1.)
                                                  ,m_isEventRecordClosed(false)
{
}

HadronizationEventRecord::HadronizationEventRecord(const Cluster& i_cluster)
                                                  :m_hadronizableClusterNumber(0)
                                                  ,m_weight(1.)
                                                  ,m_isEventRecordClosed(false)
{
    try
    {
        addInputClusters(i_cluster);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

HadronizationEventRecord::HadronizationEventRecord(const vector<Cluster>& i_clusters)
                                                  :m_hadronizableClusterNumber(0)
                                                  ,m_weight(1.)
                                                  ,m_isEventRecordClosed(false)
{
    try
    {
        addInputClusters(i_clusters);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

HadronizationEventRecord::~HadronizationEventRecord(void)
{
}

void HadronizationEventRecord::addHadronizationChannel(const HadronizationChannel& i_hadronizationChannel)
{
    // Check cluster id existence in intermediate custer set
    const unsigned int clusterIndex(i_hadronizationChannel.clusterIndex);
    const map<unsigned int,bool>::iterator mapIt(m_clusterAvailabilityFlag.find(clusterIndex));
    if(mapIt == m_clusterAvailabilityFlag.end())
    {
        // Cluster not found in available cluster set, throw exception
        throw HadronizationException("Error during hadronization channel adding, hadronized cluster id not found",
                                     __FUNCTION__,502);
    
    }
    
    // Set cluster as not available
    mapIt->second = false;
    --m_hadronizableClusterNumber;
        
    // Add hadronization channel (check for channel duplication not required)
    m_hadronizationChannels[clusterIndex] = i_hadronizationChannel;

    // Update weight
    m_weight *= i_hadronizationChannel.weight;
    
    if(m_isEventRecordClosed)
    {
        m_isEventRecordClosed = false;
    }
}

void HadronizationEventRecord::closeEventRecord(void)
{
    // Update single cluster hadronization channel particle indexes
    unsigned int particleStartIndex(m_eventClusters.rbegin()->second.getIndex()+1);
    map<unsigned int,HadronizationChannel>::const_iterator channelItEnd(m_hadronizationChannels.end());
    for(map<unsigned int,HadronizationChannel>::iterator channelIt = m_hadronizationChannels.begin();
        channelIt != channelItEnd; ++channelIt)
    {
        // Loop over channel particles
        const unsigned int numberOfParticles(channelIt->second.particles.size());
        for(unsigned int hadronIndex=0;hadronIndex<numberOfParticles;++hadronIndex)
        {
            channelIt->second.particles[hadronIndex].setIndex(particleStartIndex);
            ++particleStartIndex;
        }
    }
    m_isEventRecordClosed = true;
}

void HadronizationEventRecord::addClusterFromMerging(const Cluster& i_cluster)
{
    
    const pair<unsigned int,unsigned int> mergedPair(i_cluster.getParentIndexes());
    // Check merging pair consistency
    if(mergedPair.first == mergedPair.second)
    {
        // Duplicated merged cluster ids
        throw HadronizationException("Error during cluster merging event adding in hadronization event, merged cluster have the same ids",
                                     __FUNCTION__,503);
    }
    
    // Check presence of first merged clusters in cluster set
    map<unsigned int,bool>::iterator mapIt = m_clusterAvailabilityFlag.find(mergedPair.first);
    if(mapIt == m_clusterAvailabilityFlag.end())
    {
        // Duplicated merged cluster ids
        throw HadronizationException("Error during cluster merging event adding in hadronization event, first merged cluster id not found",
                                     __FUNCTION__,502);
    }
    
    // Set cluster as not available
    mapIt->second = false;

    // Check presence of second merged clusters in available cluster set
    mapIt = m_clusterAvailabilityFlag.find(mergedPair.second);
    if(mapIt == m_clusterAvailabilityFlag.end())
    {
        // Duplicated merged cluster ids
        throw HadronizationException("Error during cluster merging event adding in hadronization event, second merged cluster id not found",
                                     __FUNCTION__,502);
    }
    
    // Set cluster as not available
    mapIt->second = false;

    // Update cluster event set
    try
    {
        // Add new cluster
        updateEventClusterMap(i_cluster);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    // Change hadronizability satus of the merged clusters in the event cluster set
    if(m_eventClusters[mergedPair.first].isHadronizable())
    {
        m_eventClusters[mergedPair.first].setHadronizabilityFlag(false);
        --m_hadronizableClusterNumber;
    }
    if(m_eventClusters[mergedPair.second].isHadronizable())
    {
        m_eventClusters[mergedPair.second].setHadronizabilityFlag(false);
        --m_hadronizableClusterNumber;
    }
    
    // TODO: use const!
    // Update merging map
    unsigned int clusterIndex = i_cluster.getIndex();
    m_clusterMergingMap[mergedPair] = clusterIndex;
    
    // Update cluster availability map
    if(i_cluster.isHadronizable())
    {
        m_clusterAvailabilityFlag[clusterIndex] = true;
        ++m_hadronizableClusterNumber;
    }
    else
    {
        m_clusterAvailabilityFlag[clusterIndex] = false;
    }

}

void HadronizationEventRecord::addClustersFromMerging(const vector<Cluster>& i_clusters)
{
    try
    {
        for(unsigned int clusterIndex=0;clusterIndex<i_clusters.size();++clusterIndex)
        {
            addClusterFromMerging(i_clusters[clusterIndex]);
        }
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }

}

void HadronizationEventRecord::addInputClusters(const Cluster& i_cluster)
{
    // Check if the provided cluster has the correct hadronization flag
    if(!i_cluster.isHadronizable())
    {
        // Not hadronizable new cluster provided
        throw HadronizationException("Error during new cluster adding in hadronization event, provided cluster is not hadronizable",
                                     __FUNCTION__,507);
    }

    try
    {
        updateEventClusterMap(i_cluster);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
    // Update cluster availability map
    m_clusterAvailabilityFlag[i_cluster.getIndex()] = true;
    ++m_hadronizableClusterNumber;
}

void HadronizationEventRecord::addInputClusters(const vector<Cluster>& i_clusters)
{
    // Add clusters to event record checking insert failure
    for(vector<Cluster>::const_iterator vecIt = i_clusters.begin(); vecIt != i_clusters.end();++vecIt)
    {
        try
        {
            addInputClusters(*vecIt);
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
    }
}

void HadronizationEventRecord::setInputCluster(const Cluster& i_cluster)
{
    cleanHadronizationEventRecord();
    try
    {
        addInputClusters(i_cluster);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void HadronizationEventRecord::setInputClusters(const vector<Cluster>& i_clusters)
{
    cleanHadronizationEventRecord();
    try
    {
        addInputClusters(i_clusters);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }        
}

double HadronizationEventRecord::getEventWeight(void) const
{
    if(m_hadronizationChannels.size() == 0)
    {
        // No hadronization channel stored, throw exception
        throw HadronizationException("No hadronization channel stored in hadronization event, weight not available",
                                     __FUNCTION__,504);
    }
    
    if(!m_isEventRecordClosed)
    {
        // Access to hadronization channel on not closed event record, throw exceptio
        throw HadronizationException("Error during hadronization event weight retrieval, not closed event record",
                                     __FUNCTION__,509);
    }
    
    return m_weight;
}

const map<unsigned int,HadronizationChannel>& HadronizationEventRecord::getHadronizationChannels(void) const
{
    if(m_hadronizationChannels.size() == 0)
    {
        // No hadronization channel stored, throw exception
        throw HadronizationException("No hadronization channel stored in hadronization event",
                                     __FUNCTION__,504);
    }
    
    if(!m_isEventRecordClosed)
    {
        // Access to hadronization channel on not closed event record, throw exceptio
        throw HadronizationException("Error during hadronization channel retrievals, not closed event record",
                                     __FUNCTION__,509);
    }
    
    return m_hadronizationChannels;
}

const HadronizationChannel& HadronizationEventRecord::getHadronizationChannel(const unsigned int i_clusterId) const
{
    if(m_hadronizationChannels.size() == 0)
    {
        // No hadronization channel stored, throw exception
        throw HadronizationException("No hadronization channel stored in hadronization event",
                                     __FUNCTION__,504);
    }
    
    // Check cluster id existence in hadronization channel map
    map<unsigned int,HadronizationChannel>::const_iterator mapIt = m_hadronizationChannels.find(i_clusterId);
    if(mapIt == m_hadronizationChannels.end())
    {
        // Cluster not found in hadronized cluster set , throw exception
        throw HadronizationException("Error during hadronization channel retrieval, hadronized cluster id not found",
                                     __FUNCTION__,502);
    }

    if(!m_isEventRecordClosed)
    {
        // Access to hadronization channel on not closed event record, throw exceptio
        throw HadronizationException("Error during hadronization channel retrievals, not closed event record",
                                     __FUNCTION__,509);
    }

    return mapIt->second;
}

vector<Cluster> HadronizationEventRecord::getEventClusters(void) const
{
    if(m_eventClusters.size() == 0)
    {
        // No hadronization channel stored, throw exception
        throw HadronizationException("No cluster stored in hadronization event",
                                     __FUNCTION__,505);
    }
    
    vector<Cluster> o_eventClusters(m_eventClusters.size());
    unsigned int clusterIndex = 0;
    for(map<unsigned int,Cluster>::const_iterator mapIt = m_eventClusters.begin();
        mapIt != m_eventClusters.end(); ++mapIt)
    {
        o_eventClusters[clusterIndex] = mapIt->second;
        ++clusterIndex;
    }
    return o_eventClusters;
}

const Cluster& HadronizationEventRecord::getCluster(const unsigned int i_clusterId) const
{
    if(m_eventClusters.size() == 0)
    {
        // No hadronization channel stored, throw exception
        throw HadronizationException("No cluster stored in hadronization event",
                                     __FUNCTION__,505);
    }

    // Check cluster id existence in event cluster map
    map<unsigned int,Cluster>::const_iterator mapIt = m_eventClusters.find(i_clusterId);
    if(mapIt == m_eventClusters.end())
    {
        // Cluster not found in hadronized cluster set , throw exception
        throw HadronizationException("Error during event cluster retrieval, cluster id not found",
                                     __FUNCTION__,502);
    }
    
    return mapIt->second;
    
}

const map<pair<unsigned int,unsigned int>,unsigned int>& HadronizationEventRecord::getClusterMergingIdMap(void) const
{

    if(m_clusterMergingMap.size() == 0)
    {
        // No cluster merging case stored, throw exception
        throw HadronizationException("No cluster merging event stored in hadronization event, merging map not available",
                                     __FUNCTION__,506);
    }
    return m_clusterMergingMap;
    
}

vector<Cluster> HadronizationEventRecord::getHadronizableClusters(void) const
{
    if(m_hadronizableClusterNumber>0)
    {
        vector<Cluster> o_clusters(m_hadronizableClusterNumber);
        const map<unsigned int,Cluster>::const_iterator mapItEnd(m_eventClusters.end());
        unsigned int clusterIndex = 0;
        for(map<unsigned int,Cluster>::const_iterator mapIt = m_eventClusters.begin(); mapIt!=mapItEnd; ++mapIt)
        {
            if(mapIt->second.isHadronizable())
            {
                o_clusters[clusterIndex] = mapIt->second;
                ++clusterIndex;
            }
        }
        return o_clusters;
    }
    else
    {
        // No cluster available for hadronization, throw exception
        throw HadronizationException("No cluster available for hadronization",
                                     __FUNCTION__,508);
    }
}

void HadronizationEventRecord::updateEventClusterMap(const Cluster& i_newCluster)
{
    pair<map<unsigned int,Cluster>::iterator,bool> checkInsert;
        
    // Add clusters to event record checking insert failure
    checkInsert = m_eventClusters.insert(pair<unsigned int,Cluster>(i_newCluster.getIndex(),i_newCluster));
    if(!checkInsert.second)
    {
        // Duplicated cluster id, throw exception
        throw HadronizationException("Error during cluster storing in hadronization event, duplicated cluster id",
                                     __FUNCTION__,501);
    }
}

void HadronizationEventRecord::cleanHadronizationEventRecord(void)
{
    if(m_eventClusters.size())
    {
        m_eventClusters.clear();
        m_hadronizationChannels.clear();
        m_clusterAvailabilityFlag.clear();
        m_hadronizableClusterNumber = 0;
        m_clusterMergingMap.clear();
        m_weight = 1.;
        if(m_isEventRecordClosed)
        {
            m_isEventRecordClosed = false;
        }
    }
}