#ifndef HADRONIZATIONEVENTRECORD_H
#define HADRONIZATIONEVENTRECORD_H

#include "../../HadronizationChannelGenerator/Include/HadronizationChannel.h"
#include "../../HadronizationObjects/Include/Cluster.h"
#include <vector>
#include <map>
#include <set>

using namespace std;

/**
* @brief Hadronization event record storage class
*
* This class is used to store a complete hadronization event,
* from input cluster to intermediate clusters produced by the 
* to hadronization channels
*
* @author Christopher Bignamini
*/
class HadronizationEventRecord
{
    public:
    
        /**
        * @brief Constructor
        *
        */
        HadronizationEventRecord(void);

        /**
        * @brief Constructor
        *
        * @param i_cluster Input clusters initially available for hadronization procedure
        * @throw HadronizationException if the provided cluster is not hadronizable 
        */
        HadronizationEventRecord(const Cluster& i_cluster);
    
        /**
        * @brief Constructor
        *
        * @param i_clusters Input clusters initially available for hadronization procedure
        * @throw HadronizationException in case of cluster id duplication or if the provided
        *        cluster is not hadronizable
        */
        HadronizationEventRecord(const vector<Cluster>& i_clusters);
    
        // Default copy constructor and overloaded assignement operator are being used
        
        /**
        * @brief Destructor
        *
        */
        ~HadronizationEventRecord(void);
    
        // TODO: test all these methods
        /**
        * @brief Hadronization channel event add method
        *
        * Method used to add to the hadronization record a new hadronization 
        * channel for an already stored cluster. The hadronized cluster id is 
        * retrieved from the HadronizationChannel object.
        *
        * @param i_hadronizationChannel Hadronization channel
        * @throw HadronizationException in case of hadronized cluster not present in event record
        */
        void addHadronizationChannel(const HadronizationChannel& i_hadronizationChannel);
    
        /**
        * @brief Cluster add method
        *
        * Method used to add a new cluster (not coming from cluster merging)
        * to the hadronization event record. A hadronization channel can then be 
        * provided for the added cluster.
        *
        * @param i_cluster New cluster
        * @throw HadronizationException if new cluster already exists or if the provided
        *        cluster is not hadronizable
        */
        void addInputClusters(const Cluster& i_cluster);

        /**
        * @brief Cluster add method
        *
        * Method used to add a set of clusters (not coming from cluster merging)
        * to the hadronization event record. A hadronization channel can then be
        * provided for each one of the added cluster.
        *
        * @param i_clusters New clusters
        * @throw HadronizationException in case of cluster id duplication or if at least
        *        one of the provided clusters is not hadronizable
        */
        void addInputClusters(const vector<Cluster>& i_clusters);

        /**
        * @brief Cluster add (from merging) method
        *
        * Method used to add a new cluster coming from cluster merging
        * to the hadronization event record. It must be noted that the merged clusters must
        * be already present in the hadronization event record.
        *
        * @param i_cluster New cluster created from cluster merging
        * @throw HadronizationException if new cluster already exists or if cluster 
        *        parents are not found in the available cluster set or duplicated
        */
        void addClusterFromMerging(const Cluster& i_cluster);
    
        /**
        * @brief Cluster add (from merging) method
        *
        * Method used to add a set of new clusters coming from cluster merging
        * to the hadronization event record. It must be noted that the merged clusters must
        * be already present in the hadronization event record.
        *
        * @param i_clusters New clusters created from cluster merging
        * @throw HadronizationException if new clusters already exists or if cluster 
        *        parents are not found in the available cluster set or duplicated
        */
        void addClustersFromMerging(const vector<Cluster>& i_clusters);
    
        /**
        * @brief Cluster set method
        *
        * Method used to load a new cluster into the hadronization event record. 
        * A hadronization channel can then be associated to the cluster itself. 
        * It must be noted that calling this method all the other event record 
        * informations are deleted
        *
        * @param i_clusters New cluster available for hadronization procedure
        * @throw HadronizationException if the provided cluster is not hadronizable
        */
        void setInputCluster(const Cluster& i_cluster);
    
        /**
        * @brief Cluster set method
        *
        * This method can be used to load a new set of clusters into the hadronization 
        * event record. A hadronization channel can then be provided for each of them. 
        * It must be noted that calling this method all the other event record informations 
        * are deleted
        *
        * @param i_clusters New set of clusters available for hadronization procedure
        * @throw HadronizationException in case of cluster id duplication or if at least
        *        one of the provided clusters is not hadronizable
        */
        void setInputClusters(const vector<Cluster>& i_clusters);
    
        /**
        * @brief Hadronization event weight get method
        *
        * @return Total hadronization event weight
        * @throw HadronizationException if no hadronization channel has been stored in the current event record
        *        or if event record has not been closed
        */
        double getEventWeight(void) const;

        /**
        * @brief Event record closing method
        * 
        * This metod must be called at the end of the hadronization channel
        * insertion procedure, previously to any call to hadronization channel get methods
        * 
        */
        void closeEventRecord(void);
    
        /**
        * @brief Hadronization channel get method
        * 
        * @return Hadronization channel map composed of (cluster id,hadronization channel) pairs
        * @throw HadronizationException if no hadronization channel has been stored in the current 
        *        event record or if event record has not been closed
        */
        const map<unsigned int,HadronizationChannel>& getHadronizationChannels(void) const;

        /**
        * @brief Single cluster hadronization channel get method
        *
        * @return Hadronization channel corresponding to the provided cluster id
        * @throw HadronizationException if no hadronization channel has been stored in the current 
        *        event record (at all or for the provided id) or if event record has not been closed
        */
        const HadronizationChannel& getHadronizationChannel(unsigned int i_clusterId) const;

        /**
        * @brief Hadronization event record cluster get method
        *
        * @return Full set of clusters stored for the currend event record (Clusters are ordered according to their insertion sequence)
        * @throw HadronizationException if no cluster has been stored in the current event record
        */
        vector<Cluster> getEventClusters(void) const;
    
        /**
        * @brief Single cluster get method
        *
        * @param i_clusterId Cluster id
        * @return Event cluster corresponding to the provided cluster id 
        * @throw HadronizationException if no cluster has been stored in the current event record (at all or for the provided cluster id)
        */
        const Cluster& getCluster(unsigned int i_clusterId) const;
    
        /**
        * @brief Cluster merging (id) map get method
        *
        * @return Cluster merging id map composed of ((merged cluster id1,merged cluster id2),created cluster id) pairs
        * @throw HadronizationException if no cluster merging event has been stored in the current event record
        */
        const map<pair<unsigned int,unsigned int>,unsigned int>& getClusterMergingIdMap(void) const;
    
        /**
        * @brief Number of cluster merging events get method
        *
        * @return Number of stored cluster merging events
        */
        inline unsigned int getNumberOfMergingEvent(void) const { return m_clusterMergingMap.size(); }
    
        /**
        * @brief Number of input clusters (not coming from merging events) get method
        *
        * @return Number of input clusters stored in the event
        */
        inline unsigned int getNumberOfInputClusters(void) const { return (m_eventClusters.size() - m_clusterMergingMap.size()); }

        /**
        * @brief Number of hadronization channels get method
        *
        * @return Number of stored hadronization channels (or hadronized clusters)
        */
        inline unsigned int getNumberOfHadronizationChannels(void) const { return m_hadronizationChannels.size(); }    

        /**
        * @brief Number of hadronizable clusters get method
        *
        * @return Number of clusters available for hadronization
        */
        inline unsigned int getNumberOfHadronizableClusters(void) const { return m_hadronizableClusterNumber; }

        /**
        * @brief Hadronizable clusters get method
        *
        * @return Clusters available for hadronization
        * @throw HadronizationException if no hadronizable clusters exist in the current event
        */
        vector<Cluster> getHadronizableClusters(void) const;

    private:

        /**
        * @brief Event cluster set update method
        *
        * @param i_newCluster New cluster to be stored in the event record
        * @throw HadronizationException if new cluster already exists
        */
        void updateEventClusterMap(const Cluster& i_newCluster);
    
        /**
        * @brief Event record cleaning method
        */
        void cleanHadronizationEventRecord(void);

        /**
        * Full event cluster set
        */
        map<unsigned int,Cluster> m_eventClusters;
    
        /**
        * Single cluster hadronization channel map (composed of (cluster id, channel) pairs)
        */
        map<unsigned int,HadronizationChannel> m_hadronizationChannels;
    
        /**
        * Cluster availability (for hadronization or merging) map
        */
        map<unsigned int,bool> m_clusterAvailabilityFlag;

        /**
        * Number of clusters available for hadronization
        */
        unsigned int m_hadronizableClusterNumber;
    
        /**
        * Cluster merging event map composed of ((merged cluster id1,merged cluster id2),created cluster id) pairs
        */
        map<pair<unsigned int,unsigned int>,unsigned int> m_clusterMergingMap;
    
        /**
        * Total event weigth
        */
        double m_weight;
    
        /**
        * Event record status
        */
        bool m_isEventRecordClosed;    
};

#endif

