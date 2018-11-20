#ifndef HADRONIZATIONHANDLER_H
#define HADRONIZATIONHANDLER_H

#include "HadronizationEventRecord.h"
#include "ClusterMerging.h"
#include "../../HadronizationChannelGenerator/Include/ClusterHadronization.h"
#include "../../HadronizationChannelGenerator/Include/PartitionFunctionHandling.h"
#include "../../HadronizationSetup/Include/HadronizationSetup.h"
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include <vector>

// TODO: reorder include (check also other classes)
// TODO: switch to abstract "RandomGenerator" class to be used with user defined random generator

using namespace std;

/**
* @brief Hadronization handling class
*
* This class provides the statistical hadronization handler, to be used by
* the external event generator for the simulation of the hadronization step.
* After the creation of an HadronizationHadler instance, provided the required 
* setup parameter, the hadronization step during a event generator can be 
* performed with a call to the public methods "runHadronization". The generated 
* hadronization event can then be retrieved using the "getHadronizationEventRecord" 
* method. The HadronizationHandler instance must be created only once
* in the main executable.
*
* @author Christopher Bignamini
*/
class HadronizationHandler
{
    
    public:
    
        // TODO: update parameter structure and pass to the present class only those really required
        /**
        * @brief Constructor
        *
        * @param i_hadronizationSetup Hadronization handler setup (model free parameters, algorithm switches and resource adresses)
        * @throw HadronizationException in case of errors in hadronization setup 
        */
        HadronizationHandler(const HadronizationSetup& i_hadronizationSetup);
    
        // TODO: check copy constructor
        /**
        * @brief Copy constructor
        *
        * @param i_hadronizationHandler HadronizationHandler instance to be copied
        */
        HadronizationHandler(const HadronizationHandler& i_hadronizationHandler);
    
        /**
        * @brief Destructor
        *
        */
        ~HadronizationHandler(void);
    
        /**
        * @brief Cluster loading method
        *
        * This method is used for the loading of the clusters to be hadronized. 
        * Each time the loadCluster method is called the set of input clusters is 
        * rebuilt deleting the previously loaded ones.
        *
        * @param i_cluster Cluster to be hadronized
        * @throw HadronizationException in case of error during cluster loading
        */
        void loadClusters(const vector<Cluster>& i_cluster);
        
        /**
        * @brief Run hadronization
        *
        * This methods is used to perform the hadronization of a prevously loaded 
        * cluster set. The present method returns 0 in case of hadronization 
        * correctly performed for the provided input cluster set. In case of (non fatal) 
        * error during the hadronization step the exit code listed below are returned.
        * In case of exit code != 0 the corresponding event must be rejected. Fatal
        * errors are instead handled by throwing a corresponding exception.
        *  
        * @return 0 in case of hadronizazion correctly performed,
        *         1 in case of missing partition function charge configuration,
        *         2 in case of partition function interpolation mass value outside available range,
        *         3 in case of null partition function value,
        *         4 in case of maximum number of hadron sampling attempts performed
        *         5 in case of hadron unavailability for a given cluster mass
        *         6 in case of double heavy flavored cluster provided
        * @throw HadronizationException in case of error during cluster merging or hadronization
        */
        unsigned int runHadronization(void);
    
        /**
        * @brief Hadronization event record get method
        *
        * This method provides the access to the full hadronization event record generated 
        * during cluster hadronization.
        *
        * @return HadronizationEventRecord instance for the last performed hadronization procedure
        * @throw HadronizationException in case of empty HadronizationEventRecord (hadronization not performed yet)
        */
        const HadronizationEventRecord& getEventRecord(void) const;
    
        // TODO: check assignement
        /**
        * @brief Assignement operator
        *
        * @param i_hadronizationHandler HadronizationHandler to be used for assignement
        * @return Reference to the created HadronizationHandler instance
        */
        HadronizationHandler& operator=(const HadronizationHandler& i_hadronizationHandler);
        
        /**
        * @brief Cluster energy density parameter get method
        *
        * This method returns the cluster energy density currently in use for event
        * generation. 
        *
        * @return Cluster energy density parameter
        */
        inline const double& getEnergyDensity(void) const { return m_hadronizationSetup.energyDensity; }
    
        // TODO: ass missing set/get method (e.g., for parameters), if needed

    private:
    
        /**
        * Run hadronization flag
        */
        bool m_runHadronization;

        /**
        * Hadronization setup (model free parameters, algorithm switches and resource adresses)
        */
        HadronizationSetup m_hadronizationSetup;

        /**
        * Random number generator
        */
        RandomNumberGenerator* m_randomGenerator;
    
        /**
        * Hadronization channel sampling group handler
        */
        HadronSamplingGroups m_hadronSamplingGroups;
    
        /**
        * Hadronization event record
        */
        HadronizationEventRecord m_eventRecord;
    
        /**
        * Single cluster hadronization handler
        */
        ClusterHadronization m_clusterHadronization;

        /**
        * Cluster merging handler
        */
        ClusterMerging* m_clusterMerging;

};

#endif

