#ifndef CLUSTERMERGING_H
#define CLUSTERMERGING_H

#include "../../HadronizationObjects/Include/Cluster.h"
#include <vector>
#include <map>

using namespace std;

/**
* @brief Cluster pair merging data
*
* This structure is used to store the information
* of a cluster pair merging candidate
*
* @author Christopher Bignamini
*/
struct MergingData
{
    
    /**
    * @brief Constructor
    */
    MergingData(void)
               :isAvailable(false)
               ,mergingMeasure(0.)
               ,parentIndexes(0,0)
    {
    }

    /**
    * @brief Constructor
    *
    * @param i_isAvailable Availability status of the present cluster pair for final merging
    * @param i_mergingMeasure Merging measure of the present cluster pair used for merging candidate pair ordering
    * @param i_firstParentIndex Merging pair first cluster position index
    * @param i_secondParentIndex Merging pair second cluster position index
    */
    MergingData(const bool i_isAvailable,
                const double i_mergingMeasure,
                const unsigned int i_firstParentIndex,
                const unsigned int i_secondParentIndex)
               :isAvailable(i_isAvailable)
               ,mergingMeasure(i_mergingMeasure)
               ,parentIndexes(pair<unsigned int,unsigned int>(i_firstParentIndex,i_secondParentIndex))
    {
    }
    
    // Default copy constructor and assignement operator are being used
    
    /**
    * Flag storing the availability status of the present cluster pair for final merging
    */
    bool isAvailable;

    /**
    * Merging measure of the present cluster pair (e.g., 4-momenta product) used for 
    * merging candidate pair ordering
    */
    double mergingMeasure;

    /**
    * Position indexes of the clusters involved in the present merging
    */
    pair<unsigned int,unsigned int> parentIndexes;

};


/**
* @brief Cluster merging class
*
* This class is used as merging operator working on a set of externally 
* provided clusters: this gives some control on the mass spectrum of the
* hadronized cluster, in particular through the usage of the i_clusterThresholdMass,
* i_charmClusterThresholdMass and i_bottomClusterThresholdMass threshold parameters.
* The full set of clusters generated during the merging procedure itself can be 
* retrieved using the getCreatedClusters method and provided to a HadronizationEventRecord
* instance to update the current event record.
*
* @author Christopher Bignamini
*/
class ClusterMerging
{
 
    public:
    
        // TODO: refine description of the thresholds
        /**
        * @brief Constructor
        *
        * @param i_clusterThresholdMass Cluster threshold mass for merging activation (for both light and heavy flavored clusters)
        * @param i_charmClusterThresholdMass Charmed cluster threshold mass for merging activation
        * @param i_bottomClusterThresholdMass Bottom cluster threshold mass for merging activation
        * @param i_maximumClusterMass Maximum allowed cluster mass for clusters produced by merging (for both light and heavy 
        *        flavored clusters)
        * @throw HadronizationException in case of non positive maximum cluster mass provided or incompatibility between cluster 
        *        mass merging threshold and maximum value
        */
        ClusterMerging(double i_clusterThresholdMass,
                       double i_charmClusterThresholdMass,
                       double i_bottomClusterThresholdMass,
                       double i_maximumClusterMass);
        
        /**
        * @brief Destructor
        */
        ~ClusterMerging(void);
        
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Cluster threshold mass set method
        *
        * @param i_clusterThresholdMass Cluster merging threshold mass
        * @throw HadronizationException in case of incompatibility between cluster mass merging threshold and maximum value
        */
        void setTresholdMass(double i_clusterThresholdMass);
    
        /**
        * @brief Charm cluster threshold mass set method
        *
        * @param i_charmClusterThresholdMass Charm cluster merging threshold mass
        * @throw HadronizationException in case of incompatibility between charm cluster mass merging threshold and maximum value
        */
        void setCharmTresholdMass(double i_charmClusterThresholdMass);
    
        /**
        * @brief Bottom cluster threshold mass set method
        *
        * @param i_bottomClusterThresholdMass Bottom cluster merging threshold mass
        * @throw HadronizationException in case of incompatibility between bottom cluster mass merging threshold and maximum value
        */
        void setBottomTresholdMass(double i_bottomClusterThresholdMass);

        /**
        * @brief Maximum mass allowed for clusters produced by merging procedure set method
        *
        * @param i_maximumClusterMass Cluster (produced by the merging procedure) maximum mass
        * @throw HadronizationException in case of non positive maximum cluster mass provided or incompatibility between cluster
        *        mass merging threshold and maximum value
        */
        void setMaximumClusterMass(double i_maximumClusterMass);
    
        /**
        * @brief Merging run method
        *
        * Given the provided input clusters, this method performs the merging
        * procedure on the clusters themselves and making use of the previously 
        * set threshold and maximum mass values 
        *
        * @param i_inputClusters Clusters to be processed and merged
        * @throw HadronizationException if the provided cluster set is empty
        */
        void runMerging(const vector<Cluster>& i_inputClusters);
    
        /**
        * @brief Created cluster get method
        * 
        * This method provides access to the set of new clusters created
        * during the merging procedure.
        *
        * @return Vector of created clusters
        */
        vector<Cluster> getCreatedClusters(void) const;

        /**
        * @brief Created cluster number get method
        *
        * This method returns the number of new clusters created
        * during the merging procedure.
        *
        * @return Number of created clusters
        */
        inline unsigned int getCreatedClusterNumber(void) const { return m_numberOfNewClusters; }
    
    private:
    
        /**
        * @brief Cluster pair check for merging
        *
        * This method is called by the runMerging method to check if the
        * cluster pair corresponding to the provided position indexes is 
        * a candidate for merging, evaluating cluster masses and merging 
        * resulting mass.
        *
        * @param i_firstIndex First cluster position index wrt m_clusters vector
        * @param i_secondIndex Second cluster position index wrt m_clusters vector
        * @param io_mergingData Current merging pair candidate set  
        */
        void checkClusterPair(unsigned int i_firstIndex,unsigned int i_secondIndex,vector<MergingData>& io_mergingData) const;

        /**
        * Heavy cluster merging procedure activation flag
        */
        bool m_isHeavyMergingNeeded;
    
        /**
        * Cluster threshold mass for merging activation
        */
        double m_clusterThresholdMass;

        /**
        * Charm cluster threshold mass for merging activation
        */
        double m_charmClusterThresholdMass;

        /**
        * Bottom cluster threshold mass for merging activation
        */
        double m_bottomClusterThresholdMass;
    
        /**
        * Maximum allowed mass for clusters produced by merging procedure
        */
        double m_maximumClusterMass;
    
        /**
        * Merging procedure cluster set (input + clusters produced by merging iterations)
        */
        vector<Cluster> m_clusters;
        
        /**
        * Number of created clusters
        */
        unsigned int m_numberOfNewClusters;

        /**
        * Start position of created clusters in m_clusters vector
        */
        unsigned int m_newClustersStartIndex;
    
};

#endif