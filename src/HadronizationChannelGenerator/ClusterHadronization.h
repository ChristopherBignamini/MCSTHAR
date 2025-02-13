#ifndef CLUSTERHADRONIZATION_H
#define CLUSTERHADRONIZATION_H

// TODO: the not homogeneous include list which follows is a sign of bad design...
#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "PhaseSpaceSampling.h"
#include "HadronSampling.h"
#include "PartitionFunctionHandling.h"
#include "HadronizationChannel.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include <string>

using namespace std;

/**
* @brief Cluster hadronization class
*
* This class provides the methods for a single cluster hadronization.
* It is used by the HadronizationHandler class for the hadronization of
* the available clusters during event generation. 
*
* @author Christopher Bignamini
*/
class ClusterHadronization
{
    public:
    
        // TODO: the sampling parameters are strongly linked with partition function ones, add compatibility check and
        //       describe difference between phenomenological and sampling parameters
        // TODO: do we need other set/get methods? (considering the subsequent modification in HadronSampling, PartFunctionHandler, etc...)
        /**
        * @brief Constructor
        * 
        * @param i_energyDensity Cluster energy density parameter
        * @param i_gammaS Strangeness suppression parameter
        * @param i_hadronSamplingTemperature Hadron sampling temperature parameter
        * @param i_hadronSamplingEnergyDensity Hadron sampling cluster energy density sampling parameter
        * @param i_partitionFunctionDataPath Partition function data set path
        * @param i_hadronSamplingGroups Hadron sampling group handler 
        * @param io_randomGenerator Random number generator (I/O parameter)
        * @throw HadronizationException in case of sampling/physical parameter outside the allowed ranges 
        *        or error during hadron sampling/partition funciton data handling setup
        */
        ClusterHadronization(double i_energyDensity,
                             double i_gammaS,
                             double i_hadronSamplingTemperature,
                             double i_hadronSamplingEnergyDensity,
                             const string& i_partitionFunctionDataPath,
                             const HadronSamplingGroups& i_hadronSamplingGroups,
                             RandomNumberGenerator& io_randomGenerator);

        /**
        * @brief Destructor
        */
        ~ClusterHadronization(void);

        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Cluster hadronization method
        * 
        * This method performs the microcanonical hadronization of the provided cluster.
        * The generated hadronization channel can be retrieved making use of the getHadronizationChannel
        * method. The present method returns 0 in case of hadronization correctly performed 
        * for the provided cluster. In case of (non fatal) error during the hadronization 
        * step the exit code listed below are returned. In case of exit code != 0 the 
        * corresponding event must be rejected. Fatal errors are instead handled by throwing a 
        * corresponding exception.
        *
        * @param i_cluster Cluster to be hadronized
        * @return 0 in case of hadronizazion correctly performed,
        *         1 in case of missing partition function charge configuration,
        *         2 in case of partition function interpolation mass value outside available range,
        *         3 in case of null partition function value,
        *         4 in case of maximum number of hadron sampling attempts performed
        *         5 in case of hadron unavailability for a given cluster mass
        *         6 in case of double heavy flavored cluster provided
        * @throw HadronizationException in case of ClusterHadronization object not correctly set,
        *        errors during hadronization channel sampling, null cluster partition function value 
        *        or errors during the retrieval of the partition function value itself
        */
        unsigned int runHadronization(const Cluster& i_cluster);
    
        /**
        * @brief Cluster hadronization channel get method
        *
        * @return Hadronization channel
        * @throw HadronizationException in case of hadronization channel unavailability
        */
        const HadronizationChannel& getHadronizationChannel(void) const;
    
    private:
    
        /**
        * @brief Cluster hadronization event building method
        *
        * @param i_clusterIndex Cluster index (object index in event record)
        * @param i_partitionFunctionValue Partition funciton value for the hadronized cluster 
        * @throw HadronizationException in case of error during hadron set or phase space configuration retrieval
        */
        void buildHadronizationEvent(unsigned int i_clusterIndex,double i_partitionFunctionValue);
    
        /**
        * Hadronization channel availability flag
        */
        bool m_isHadronizationChannelAvailable;

        /**
        * Cluster hadronization setup status flag
        */
        bool m_isClusterHadronizationReady;

        /**
        * Cluster hadronization channel
        */
        HadronizationChannel m_hadronizationChannel;
    
        /**
        * Hadron sampling handler
        */
        HadronSampling m_hadronSampling;
    
        /**
        * Phase space sampling handler
        */
        PhaseSpaceSampling m_phaseSpaceSampling;
    
        /**
        * Strangeness suppression parameter
        */
        double m_gammaS;
    
        /**
        * Partition function data handling
        */
        PartitionFunctionHandling m_partitionFunction;
    
};

#endif