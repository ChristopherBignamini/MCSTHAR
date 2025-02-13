#include "../Include/ClusterMerging.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include <cassert>

bool mergingDataComparison(const MergingData& i_firstMerging,const MergingData& i_secondMerging)
{
    return (i_firstMerging.mergingMeasure<i_secondMerging.mergingMeasure);
}

ClusterMerging::ClusterMerging(const double i_clusterThresholdMass,
                               const double i_charmClusterThresholdMass,
                               const double i_bottomClusterThresholdMass,
                               const double i_maximumClusterMass)
                              :m_isHeavyMergingNeeded(false)
                              ,m_clusterThresholdMass(i_clusterThresholdMass)
                              ,m_charmClusterThresholdMass(i_charmClusterThresholdMass)
                              ,m_bottomClusterThresholdMass(i_bottomClusterThresholdMass)
                              ,m_maximumClusterMass(i_maximumClusterMass)
                              ,m_numberOfNewClusters(0)
                              ,m_newClustersStartIndex(0)
{
    // TODO: what to do with "=" comparisons with fp variables?
    if(m_clusterThresholdMass>=m_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }
    if(m_charmClusterThresholdMass>=m_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, charm cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }
    if(m_bottomClusterThresholdMass>=m_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, bottom cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }
    if(m_maximumClusterMass<=0.)
    {
        throw HadronizationException("Error during cluster merging setting, non positive cluster maximum mass provided",
                                     __FUNCTION__,521);
    }
}

ClusterMerging::~ClusterMerging(void)
{
}

void ClusterMerging::setTresholdMass(const double i_clusterThresholdMass)
{
    if(i_clusterThresholdMass >= m_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }

    m_clusterThresholdMass = i_clusterThresholdMass;
}

void ClusterMerging::setCharmTresholdMass(const double i_charmClusterThresholdMass)
{
    if(i_charmClusterThresholdMass>=m_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, charm cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }

    m_charmClusterThresholdMass = i_charmClusterThresholdMass;
}

void ClusterMerging::setBottomTresholdMass(const double i_bottomClusterThresholdMass)
{
    if(i_bottomClusterThresholdMass >= m_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, bottom cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }

    m_bottomClusterThresholdMass = i_bottomClusterThresholdMass;
}

void ClusterMerging::setMaximumClusterMass(const double i_maximumClusterMass)
{
    if(m_clusterThresholdMass >= i_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }
    if(m_charmClusterThresholdMass >= i_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, charm cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }
    if(m_bottomClusterThresholdMass >= i_maximumClusterMass)
    {
        throw HadronizationException("Error during cluster merging setting, bottom cluster merging threshold value is larger than cluster maximum mass",__FUNCTION__,521);
    }
    if(i_maximumClusterMass<=0.)
    {
        throw HadronizationException("Error during cluster merging setting, non positive cluster maximum mass provided",
                                     __FUNCTION__,521);
    }

    m_maximumClusterMass = i_maximumClusterMass;
}

void ClusterMerging::runMerging(const vector<Cluster>& i_inputClusters)
{
    // Empty cluster vector provided
    if(i_inputClusters.size()==0)
    {
        throw HadronizationException("Error during cluster merging, provided input cluster set is empty",
                                     __FUNCTION__,522);
    }
    
    // Store input clusters
    m_clusters = i_inputClusters;
        
    // Store new clusters start global index
    unsigned int newClusterGlobalIndex = (m_clusters.end()-1)->getIndex() + 1;
    unsigned int numberOfClusters = m_clusters.size();
    vector<unsigned int> availableClusterIndexes(numberOfClusters);
        
    // Store available cluster indexes and check if heavy flavored cluster merging is required
    m_isHeavyMergingNeeded = false;
    for(unsigned int clusterIndex = 0;clusterIndex<numberOfClusters;++clusterIndex)
    {
        // Check if heavy flavored cluster merging procedure must be activated
        if(m_isHeavyMergingNeeded==false)
        {
            if(m_clusters[clusterIndex].isHeavyFlavored())
            {
                if((m_clusters[clusterIndex].getCharmCharge() !=0) &&
                   (m_clusters[clusterIndex].getMass()<m_charmClusterThresholdMass))
                {                    
                    m_isHeavyMergingNeeded = true;
                }
                else if((m_clusters[clusterIndex].getBottomCharge() !=0) && (m_clusters[clusterIndex].getMass()<m_bottomClusterThresholdMass))
                {
                    m_isHeavyMergingNeeded = true;
                }
            }
        }
        
        // Update available cluster index vector
        availableClusterIndexes[clusterIndex] = clusterIndex;
    }
    
    // Start merging procedure
    // Reset number of new clusters
    m_numberOfNewClusters = 0;
    // Within the do-loop a first run is performed for heavy flavored cluster merging only,
    // if needed, and a second run over the full cluster set. The handling of the two
    // cluster merging runs is realized through the m_isHeavyMergingNeeded flag. In particular,
    // the checkClusterPair is able to identify the two merging runs using the above flag.
    vector<MergingData> mergingData;
    unsigned int numberOfMerging;
    bool runNewIteration = true;
    do
    {
        mergingData.clear();
        numberOfMerging = 0;
        unsigned int availableClusterNumber = availableClusterIndexes.size();

        // Loop over available cluster pairs (merging candidate vector mergingData is updated within checkClusterPair method)
        for(unsigned int firstIndex=0;firstIndex<availableClusterNumber-1;++firstIndex)
        {
            for(unsigned int secondIndex=firstIndex+1;secondIndex<availableClusterNumber;++secondIndex)
            {
                checkClusterPair(availableClusterIndexes[firstIndex],availableClusterIndexes[secondIndex],mergingData);
            }
        }
        
        // Retrieve number of merging pair candidates
        unsigned int mergingCandidateNumber = mergingData.size();
        if(mergingCandidateNumber>0)
        {
            // Sort merging candidate data using merging measures
            sort(mergingData.begin(),mergingData.end(),mergingDataComparison);

            numberOfMerging = 0;
            if(mergingCandidateNumber == 1 and mergingData[0].isAvailable)
            {
                // Update number of merging to be performed
                ++numberOfMerging;
            }
            else
            {
                // Delete multiple cluster merging and store number of merging to be performed
                pair<unsigned int,unsigned int> firstParentIndexes;
                pair<unsigned int,unsigned int> secondParentIndexes;
                // Loop over merging candidate pair
                for(unsigned int firstMergingIndex=0;firstMergingIndex<mergingCandidateNumber;++firstMergingIndex)
                {
                    // Check if merging case is available
                    if(mergingData[firstMergingIndex].isAvailable)
                    {
                        // Store merging case cluster indexes
                        firstParentIndexes = mergingData[firstMergingIndex].parentIndexes;

                        // Loop over subsequent merging cases
                        for(unsigned int secondMergingIndex=firstMergingIndex+1;secondMergingIndex<mergingCandidateNumber;
                            ++secondMergingIndex)
                        {
                            
                            // Check if merging case is available
                            if(mergingData[secondMergingIndex].isAvailable)
                            {
                                // Store merging case cluster indexes
                                secondParentIndexes = mergingData[secondMergingIndex].parentIndexes;

                                // Check the presence of at least one of first merging case
                                // clusters into the second merging case
                                if(firstParentIndexes.first == secondParentIndexes.first ||
                                   firstParentIndexes.first == secondParentIndexes.second ||
                                   firstParentIndexes.second == secondParentIndexes.first ||
                                   firstParentIndexes.second == secondParentIndexes.second)
                                {
                                    // At least one of the clusters involved in the second merging case
                                    // is present also in the first one. Set the second case as unavailable
                                    // for merging, so that the cluster pair with smaller merging measure
                                    // will be merged
                                    mergingData[secondMergingIndex].isAvailable = false;

                                }
                            }
                        }
                        // Update number of merging to be performed
                        ++numberOfMerging;
                    }
                }
            }
            
            // Perform cluster merging
            if(numberOfMerging>0)
            {
                
                // Resize cluster vector
                m_clusters.resize(numberOfClusters + numberOfMerging);
                unsigned int newClusterIndex = numberOfClusters;
                for(unsigned int mergingIndex = 0;mergingIndex<mergingCandidateNumber;++mergingIndex)
                {
                    
                    // Check if merging event has been accepted
                    if(mergingData[mergingIndex].isAvailable)
                    {                        
                        // Retrieve parent indexes
                        pair<unsigned int,unsigned int> parentIndexes(mergingData[mergingIndex].parentIndexes);
                        
                        // Build new cluster and update cluster set
                        Cluster* firstParentCluster(&m_clusters[parentIndexes.first]);
                        Cluster* secondParentCluster(&m_clusters[parentIndexes.second]);
                        m_clusters[newClusterIndex] = *firstParentCluster + *secondParentCluster;
                        
                        // Set new cluster global index
                        m_clusters[newClusterIndex].setIndex(newClusterGlobalIndex);
                        ++newClusterGlobalIndex;
                                                
                        // Remove indexes of the created cluster parents from the list of cluster available for merging
                        // Remove first parent
                        vector<unsigned int>::iterator vecIt =
                            find(availableClusterIndexes.begin(),availableClusterIndexes.end(),parentIndexes.first);
                        assert(vecIt != availableClusterIndexes.end());
                        availableClusterIndexes.erase(vecIt);// TODO: avoid this erase, hadronizability flag is enough to build candidate pairs!!!
                        // Set first parent hadronizability flag to false
                        firstParentCluster->setHadronizabilityFlag(false);
                                   
                        // Remove second parent
                        vecIt = find(availableClusterIndexes.begin(),availableClusterIndexes.end(),
                                     parentIndexes.second);
                        assert(vecIt != availableClusterIndexes.end());
                        availableClusterIndexes.erase(vecIt);
                        // Set second parent hadronizability flag to false
                        secondParentCluster->setHadronizabilityFlag(false);
                        
                        // Add new cluster index to availability vector
                        availableClusterIndexes.push_back(newClusterIndex);
                        
                        // Update cluster index
                        ++newClusterIndex;
                                                
                    }
                }
                // Update number of stored clusters
                numberOfClusters += numberOfMerging;
            }
        }
                
        // Check which kind of merging is needed for next iteration, if any
        if(numberOfMerging == 0)
        {
            // No merging performed, check if the iteration was on
            // heavy flavored cluster only
            if(m_isHeavyMergingNeeded)
            {
                // Last iteration was on heavy flavored clusters only, continue with full cluster set
                m_isHeavyMergingNeeded = false;
            }
            else
            {
                // Last iteration was on full cluster set, no further
                // iterations needed
                runNewIteration = false;
            }
        }
        
        // Update number of new clusters
        m_numberOfNewClusters += numberOfMerging;
        
    }while(runNewIteration);
        
    // Store position for new clusters created by merging procedure
    if(m_numberOfNewClusters>0)
    {
        m_newClustersStartIndex = i_inputClusters.size();
    }
}

vector<Cluster> ClusterMerging::getCreatedClusters(void) const
{
    if(m_numberOfNewClusters>0)
    {
        vector<Cluster> o_newClusters(m_clusters.begin()+m_newClustersStartIndex,m_clusters.end());
        return o_newClusters;
    }
    else
    {
        throw HadronizationException("Error during cluster merging data retrieval, no new cluster created during merging procedure",
                                     __FUNCTION__,523);
    }
}

void ClusterMerging::checkClusterPair(const unsigned int i_firstIndex,
                                      const unsigned int i_secondIndex,
                                      vector<MergingData>& io_mergingData) const
{
    bool addMergingPair = false;
    double mergingMeasure;
    
    const Cluster* firstCluster(&m_clusters[i_firstIndex]);
    const Cluster* secondCluster(&m_clusters[i_secondIndex]);
    
    // Check if merging is currently working on heavy flavored cluster only
    if(m_isHeavyMergingNeeded)
    {        
        // Find cluster pair heavy flavor composition: only pairs
        // with only one heavy flavored cluster will be included in merging
        // candidate list.
        const bool isFirstClusterHeavyFlavored(firstCluster->isHeavyFlavored());
        const bool isSecondClusterHeavyFlavored(secondCluster->isHeavyFlavored());
        
        if(isFirstClusterHeavyFlavored && !isSecondClusterHeavyFlavored)
        {
            // First cluster is the heavy flavored one, check if its mass
            // is smaller than the specified flavor dependent thresholds
            if(firstCluster->getCharmCharge())
            {
                if(firstCluster->getMass()<m_charmClusterThresholdMass)
                {
                    addMergingPair = true;
                }
            }
            else
            {
                if(firstCluster->getMass()<m_bottomClusterThresholdMass)
                {
                    addMergingPair = true;
                }
            }
        }
        else if(!isFirstClusterHeavyFlavored && isSecondClusterHeavyFlavored)
        {
            // Second cluster is the heavy flavored one, check if its mass
            // is smaller than the specified flavor dependent thresholds
            if(secondCluster->getCharmCharge())
            {
                if(secondCluster->getMass()<m_charmClusterThresholdMass)
                {
                    addMergingPair = true;
                }
            }
            else
            {
                if(secondCluster->getMass()<m_bottomClusterThresholdMass)
                {
                    addMergingPair = true;
                }
            }
        }
    }
    else
    {
        // Merging is currently working on full cluster set, check if at least
        // one of the cluster composing the provided pair has mass smaller than the
        // specified threshold
        if((firstCluster->getMass()<m_clusterThresholdMass) || (secondCluster->getMass()<m_clusterThresholdMass))
        {
            addMergingPair = true;
        }
    }
        
    if(addMergingPair)
    {
        // New merging candidate pair must be included in the candidate list
        // Check whether the resulting cluster would have mass larger than the
        // specified maximum cluster mass
        if(TLorentzVector(firstCluster->getP()+secondCluster->getP()).Mag()<m_maximumClusterMass)
        {
            // Compute cluster pair merging measure (4-momenta scalar product)
            mergingMeasure = firstCluster->getP()*secondCluster->getP();
            
            // Update merging candidate list with the new pair
            io_mergingData.push_back(MergingData(true,mergingMeasure,i_firstIndex,i_secondIndex));
        }
    }
}
