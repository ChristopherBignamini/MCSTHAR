#include "../Include/updateEventRecord.h"
#include "../Include/Herwig6510Wrapper.h"
#include "../Include/Herwig6510Constants.h"
#include "../Include/convertParticleId.h"

void updateEventRecord(const HadronizationEventRecord& i_eventRecord)
{
    // Retrieve current number of event objects
    unsigned int hepevtNHEP(static_cast<unsigned int>(MCSTHAR::hepevt.NHEP));
    // Additional particle number storage variable (see below)
    unsigned int numberOfExtraParticles = 0;
    
    // Check if additional cluster must be added to event record due to merging procedure
    const unsigned int numberOfAdditionalClusters(i_eventRecord.getNumberOfMergingEvent());
    if(numberOfAdditionalClusters)
    {
        
        // Shift already present particles in event record after cluster block (see hwcfor, loop 160) to make
        // place for new clusters produced during mergin procedure
        const unsigned int extraParticleStartIndex(i_eventRecord.getNumberOfInputClusters() +
                                                   i_eventRecord.getEventClusters()[0].getIndex());
        
        // Retrieve number of additional particles
        numberOfExtraParticles = hepevtNHEP - extraParticleStartIndex;
        if(numberOfExtraParticles>0)
        {
            // Shift additional particles in event record 
            unsigned int extraParticleIndex(extraParticleStartIndex);
            unsigned int extraParticleNewIndex(extraParticleStartIndex+numberOfAdditionalClusters);
            for(;extraParticleIndex<hepevtNHEP;++extraParticleNewIndex,++extraParticleIndex)
            {
                // TODO: check if other fields must be updated
                (MCSTHAR::hepevt.ISTHEP)[extraParticleNewIndex] = (MCSTHAR::hepevt.ISTHEP)[extraParticleIndex];
                (MCSTHAR::hepevt.IDHEP)[extraParticleNewIndex] = (MCSTHAR::hepevt.IDHEP)[extraParticleIndex];
                (hwevnt.IDHW)[extraParticleNewIndex] = (hwevnt.IDHW)[extraParticleIndex];
                (MCSTHAR::hepevt.JMOHEP)[extraParticleNewIndex][0] = (MCSTHAR::hepevt.JMOHEP)[extraParticleIndex][0];
                (MCSTHAR::hepevt.JMOHEP)[extraParticleNewIndex][1] = (MCSTHAR::hepevt.JMOHEP)[extraParticleIndex][1];
                (MCSTHAR::hepevt.JDAHEP)[extraParticleNewIndex][0] = (MCSTHAR::hepevt.JDAHEP)[extraParticleIndex][0];
                (MCSTHAR::hepevt.JDAHEP)[extraParticleNewIndex][1] = (MCSTHAR::hepevt.JDAHEP)[extraParticleIndex][1];
                (MCSTHAR::hepevt.PHEP)[extraParticleNewIndex][0] = (MCSTHAR::hepevt.PHEP)[extraParticleIndex][0];
                (MCSTHAR::hepevt.PHEP)[extraParticleNewIndex][1] = (MCSTHAR::hepevt.PHEP)[extraParticleIndex][1];
                (MCSTHAR::hepevt.PHEP)[extraParticleNewIndex][2] = (MCSTHAR::hepevt.PHEP)[extraParticleIndex][2];
                (MCSTHAR::hepevt.PHEP)[extraParticleNewIndex][3] = (MCSTHAR::hepevt.PHEP)[extraParticleIndex][3];
                (MCSTHAR::hepevt.PHEP)[extraParticleNewIndex][4] = (MCSTHAR::hepevt.PHEP)[extraParticleIndex][4];
                (MCSTHAR::hepevt.VHEP)[extraParticleNewIndex][0] = (MCSTHAR::hepevt.VHEP)[extraParticleIndex][0];
                (MCSTHAR::hepevt.VHEP)[extraParticleNewIndex][1] = (MCSTHAR::hepevt.VHEP)[extraParticleIndex][1];
                (MCSTHAR::hepevt.VHEP)[extraParticleNewIndex][2] = (MCSTHAR::hepevt.VHEP)[extraParticleIndex][2];
                (MCSTHAR::hepevt.VHEP)[extraParticleNewIndex][3] = (MCSTHAR::hepevt.VHEP)[extraParticleIndex][3];
            }        
        }
        
        // Add cluster coming from merging to event record
        const map< pair<unsigned int, unsigned int>, unsigned int> mergingMap(i_eventRecord.getClusterMergingIdMap());
        for(map< pair<unsigned int, unsigned int>, unsigned int>::const_iterator mapIt = mergingMap.begin();
            mapIt != mergingMap.end(); ++mapIt)
        {
            // TODO: avoid all minus/plus one
            const unsigned int firstClusterIndex(mapIt->first.first);
            const unsigned int secondClusterIndex(mapIt->first.second);

            // Set merged cluster status
            // TODO: it should be 163, but it conflict with heavy hadron decay... If this is the final
            // strategy the subsequent switch to 183 for hadronized cluster can be avoided.
            (MCSTHAR::hepevt.ISTHEP)[firstClusterIndex] = herwigHadronizedClusterStatus;
            // TODO: optimize using pointers everywhere for structure access
            (MCSTHAR::hepevt.JDAHEP)[firstClusterIndex][0] = mapIt->second+1;
            (MCSTHAR::hepevt.JDAHEP)[firstClusterIndex][1] = mapIt->second+1;
            (MCSTHAR::hepevt.ISTHEP)[secondClusterIndex] = herwigHadronizedClusterStatus;
            (MCSTHAR::hepevt.JDAHEP)[secondClusterIndex][0] = mapIt->second+1;
            (MCSTHAR::hepevt.JDAHEP)[secondClusterIndex][1] = mapIt->second+1;
                        
            // Insert new cluster in Herwig event record
            // TODO: implement "addCluster" function to be reused with particles too
            const unsigned int newClusterIndex(mapIt->second);
            const Cluster newCluster(i_eventRecord.getCluster(mapIt->second));
            (MCSTHAR::hepevt.ISTHEP)[newClusterIndex] = herwigHadronizedClusterStatus;
            (MCSTHAR::hepevt.IDHEP)[newClusterIndex] = newCluster.getIdCode();
            (hwevnt.IDHW)[newClusterIndex] = herwigInternalClusterIdCode;
			(MCSTHAR::hepevt.JMOHEP)[newClusterIndex][0] = newCluster.getParentIndexes().first+1;
			(MCSTHAR::hepevt.JMOHEP)[newClusterIndex][1] = newCluster.getParentIndexes().second+1;
			(MCSTHAR::hepevt.PHEP)[newClusterIndex][0] = newCluster.getP().Px();
			(MCSTHAR::hepevt.PHEP)[newClusterIndex][1] = newCluster.getP().Py();
			(MCSTHAR::hepevt.PHEP)[newClusterIndex][2] = newCluster.getP().Pz();
			(MCSTHAR::hepevt.PHEP)[newClusterIndex][3] = newCluster.getP().E();
			(MCSTHAR::hepevt.PHEP)[newClusterIndex][4] = newCluster.getMass();
            
            // At present the clusters produced during the merging procedure are created in the space-time position of the
            // first mother cluster
            // TODO: after check with old code evaluate the effect of different definitions and discuss with FB
			(MCSTHAR::hepevt.VHEP)[newClusterIndex][0] = (MCSTHAR::hepevt.VHEP)[firstClusterIndex][0];
			(MCSTHAR::hepevt.VHEP)[newClusterIndex][1] = (MCSTHAR::hepevt.VHEP)[firstClusterIndex][1];
			(MCSTHAR::hepevt.VHEP)[newClusterIndex][2] = (MCSTHAR::hepevt.VHEP)[firstClusterIndex][2];
			(MCSTHAR::hepevt.VHEP)[newClusterIndex][3] = (MCSTHAR::hepevt.VHEP)[firstClusterIndex][3];
        }

        // Update number or event record entries
        MCSTHAR::hepevt.NHEP += numberOfAdditionalClusters;
        hepevtNHEP += numberOfAdditionalClusters;
    }
    
    // Insert hadronization channels in event record
    const map<unsigned int,HadronizationChannel> hadronizationChannels(i_eventRecord.getHadronizationChannels());
    unsigned int newHadronIndex = hepevtNHEP;
    unsigned int numberOfHadrons = 0;
    for(map<unsigned int,HadronizationChannel>::const_iterator mapIt = hadronizationChannels.begin();
        mapIt != hadronizationChannels.end();++mapIt)
    {
        const unsigned int numberOfParticles(mapIt->second.particles.size());
        for(unsigned int hadronIndex=0;hadronIndex<numberOfParticles;++hadronIndex)
        {
            const Particle newParticle(mapIt->second.particles[hadronIndex]);

            (MCSTHAR::hepevt.ISTHEP)[newHadronIndex] = herwigUnstableHadronFromClusterStatus;
            (MCSTHAR::hepevt.IDHEP)[newHadronIndex] = newParticle.getIdCode();
            (hwevnt.IDHW)[newHadronIndex] = 0;
            (MCSTHAR::hepevt.JMOHEP)[newHadronIndex][0] = newParticle.getParentIndex()+1;
            (MCSTHAR::hepevt.JMOHEP)[newHadronIndex][1] = (MCSTHAR::hepevt.JMOHEP)[newHadronIndex][0];
            (MCSTHAR::hepevt.JDAHEP)[newHadronIndex][0] = 0;
            (MCSTHAR::hepevt.JDAHEP)[newHadronIndex][1] = 0;
            // Update phase space and production vertex of the new particles
            (MCSTHAR::hepevt.PHEP)[newHadronIndex][0] = newParticle.getP().Px();
            (MCSTHAR::hepevt.PHEP)[newHadronIndex][1] = newParticle.getP().Py();
            (MCSTHAR::hepevt.PHEP)[newHadronIndex][2] = newParticle.getP().Pz();
            (MCSTHAR::hepevt.PHEP)[newHadronIndex][3] = newParticle.getP().E();
            (MCSTHAR::hepevt.PHEP)[newHadronIndex][4] = newParticle.getMass();
            // At present the particle are supposed to be created in the same space-time point
            // of the hadronizing cluster creaction
            // TODO: same as above for cluster creation
            const unsigned int parentIndex(newParticle.getParentIndex());
            (MCSTHAR::hepevt.VHEP)[newHadronIndex][0] = (MCSTHAR::hepevt.VHEP)[parentIndex][0];
            (MCSTHAR::hepevt.VHEP)[newHadronIndex][1] = (MCSTHAR::hepevt.VHEP)[parentIndex][1];
            (MCSTHAR::hepevt.VHEP)[newHadronIndex][2] = (MCSTHAR::hepevt.VHEP)[parentIndex][2];
            (MCSTHAR::hepevt.VHEP)[newHadronIndex][3] = (MCSTHAR::hepevt.VHEP)[parentIndex][3];
            
            ++newHadronIndex;
        }
        
        // Update number of inserted hadrons
        numberOfHadrons += numberOfParticles;
        
        // Update hadronized cluster daughters and status
        (MCSTHAR::hepevt.JDAHEP)[mapIt->first][0] = mapIt->second.particles[0].getIndex()+numberOfExtraParticles+1;
        (MCSTHAR::hepevt.JDAHEP)[mapIt->first][1] = mapIt->second.particles[numberOfParticles-1].getIndex()+numberOfExtraParticles+1;
        (MCSTHAR::hepevt.ISTHEP)[mapIt->first] = herwigHadronizedClusterStatus;
        
    }
    
    // Update the number of record entries
	MCSTHAR::hepevt.NHEP += numberOfHadrons;
    
    // Update the global event weight
	hwevnt.EVWGT *= i_eventRecord.getEventWeight();
    
    // Set Herwig hadron id codes
    convertParticleId();
    
    return;
}