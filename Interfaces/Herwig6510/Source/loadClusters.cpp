#include "../Include/loadClusters.h"
#include "../Include/ChargeData.h"
#include "../Include/Herwig6510Wrapper.h"
#include "../Include/Herwig6510Constants.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/Constants.h"
#include "TLorentzVector.h"

vector<Cluster> loadClusters(const double i_energyDensity)
{
    // Store initial space for a cluster only event, to be resized later on (to avoid push back)
    vector<Cluster> clusters(MCSTHAR::hepevt.NHEP);
    unsigned int numberOfClusters = 0;
    int strangeCharge = 0;
    int charmCharge = 0;
    int bottomCharge = 0;
    double electricCharge = 0.;
    double baryonicCharge = 0.;
    
    // Loop over event record entries
    // TODO: optimize code using pointer instead of matrices
    for(int objectIndex=0;objectIndex<MCSTHAR::hepevt.NHEP;++objectIndex)
    {
        // TODO: check herwig code for clusters (164/165)
        // Check if current object is a cluster
        if(MCSTHAR::hepevt.IDHEP[objectIndex] == herwigClusterIdCode)
        {
            const int clusterStatus((MCSTHAR::hepevt.ISTHEP)[objectIndex]);
            if((clusterStatus != herwigBeamClusterBeforeSoftProcessStatus) &&
               (clusterStatus != herwigTargetClusterBeforeSoftProcessStatus))
            {
                const int firstParentIndex((MCSTHAR::hepevt.JMOHEP)[objectIndex][0]-1);
                const int secondParentIndex((MCSTHAR::hepevt.JMOHEP)[objectIndex][1]-1);
                
                try
                {
                    // TODO: store mother indexes for use below
                    retrieveClusterCharges(firstParentIndex,
                                           secondParentIndex,
                                           strangeCharge,
                                           charmCharge,
                                           bottomCharge,
                                           electricCharge,
                                           baryonicCharge);
                }
                catch(HadronizationException& ex)
                {
                    throw ex;
                }
                
                // Add new cluster to cluster set
                clusters[numberOfClusters] =
                    Cluster(true,
                            make_pair(static_cast<unsigned int>(firstParentIndex),static_cast<unsigned int>(secondParentIndex)),
                            strangeCharge,charmCharge,bottomCharge,electricCharge,baryonicCharge,i_energyDensity,
                            TLorentzVector((MCSTHAR::hepevt.PHEP)[objectIndex][0],
                                           (MCSTHAR::hepevt.PHEP)[objectIndex][1],
                                           (MCSTHAR::hepevt.PHEP)[objectIndex][2],
                                           (MCSTHAR::hepevt.PHEP)[objectIndex][3]),
                            objectIndex);
                
                ++numberOfClusters;
            }
        }
    }
    
    // Return resized cluster vector
    return vector<Cluster>(clusters.begin(),clusters.begin() + numberOfClusters);
}

void retrieveClusterCharges(const int i_firstParentIndex,
                            const int i_secondParentIndex,
                            int& io_strangeCharge,
                            int& io_charmCharge,
                            int& io_bottomCharge,
                            double& io_electricCharge,
                            double& io_baryonicCharge)
{
    // Retrieve cluster parents ID codes
    const int firstParentIDCode(MCSTHAR::hepevt.IDHEP[i_firstParentIndex]);
    const int secondParentIDCode(MCSTHAR::hepevt.IDHEP[i_secondParentIndex]);
    int firstParentListIndex = 0;
    int firstParentPDGCode = 0;
    int secondParentListIndex = 0;
    int secondParentPDGCode = 0;
    
    // Retrieve cluster parents IDPDG codes
    unsigned int foundParentNumber=0;
    for(int particleIndex=0;particleIndex<nmxres;++particleIndex)
    {
        // TODO: can't we avoid this step??? and use herwig id???
        // TODO: herwig structure starts from zero, which kind of alignement do we expect???
        const int particlePDGCode(hwprop.IDPDG[particleIndex]);
        
        if(firstParentIDCode == particlePDGCode)
        {
            firstParentListIndex = particleIndex;
            firstParentPDGCode = particlePDGCode;
            ++foundParentNumber;
            if(foundParentNumber==2)
            {
                break;
            }
        }
        
        if(secondParentIDCode == particlePDGCode)
        {
            secondParentListIndex = particleIndex;
            secondParentPDGCode = particlePDGCode;
            ++foundParentNumber;
            if(foundParentNumber==2)
            {
                break;
            }
        }
    }
    
    // Check if parent pair has been identified
    if(foundParentNumber<2)
    {
        throw HadronizationException("Error during Herwig cluster charge retrieval, constituents not identified",
                                     __FUNCTION__,601);
    }
    
    // Compute parent pair total charges
    io_electricCharge = 0;
    io_baryonicCharge = 0.;
    io_strangeCharge = 0;
    io_charmCharge = 0;
    io_bottomCharge = 0;
    
    // Find first parent charge constributions
    map<int,ChargeData>::const_iterator mapIt(quarkChargesMap.find(firstParentPDGCode));
    if(mapIt != quarkChargesMap.end())
    {
        io_electricCharge = mapIt->second.electricCharge;
        io_baryonicCharge = mapIt->second.baryonicCharge;
        io_strangeCharge = mapIt->second.strangeCharge;
        io_charmCharge = mapIt->second.charmCharge;
        io_bottomCharge = mapIt->second.bottomCharge;    
    }
    else
    {
        throw HadronizationException("Error during Herwig cluster charge retrieval, not identified first constituent",
                                     __FUNCTION__,602);
    
    }

    // Find second parent charge constributions
    mapIt = quarkChargesMap.find(secondParentPDGCode);
    if(mapIt != quarkChargesMap.end())
    {
        io_electricCharge += mapIt->second.electricCharge;
        io_baryonicCharge += mapIt->second.baryonicCharge;
        io_strangeCharge += mapIt->second.strangeCharge;
        io_charmCharge += mapIt->second.charmCharge;
        io_bottomCharge += mapIt->second.bottomCharge;
    }
    else
    {
        throw HadronizationException("Error during Herwig cluster charge retrieval, not identified second constituent",
                                     __FUNCTION__,603);
    }
}
