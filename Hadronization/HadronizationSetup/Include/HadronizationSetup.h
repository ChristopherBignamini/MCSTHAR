#ifndef HADRONIZATIONSETUP_H
#define HADRONIZATIONSETUP_H

#include <string>

using namespace std;

/**
* @brief Hadronization setup parameter storage structure
*
* @author Christopher Bignamini
*/
struct HadronizationSetup
{
    /**
    * @brief Constructor
    *
    */
    // TODO: use constants for default parameters
    HadronizationSetup(void)
                      :energyDensity(0.35)
                      ,gammaS(0.65)
                      ,clusterMergingFlag(true)
                      ,clusterMergingMinMass(1.8)
                      ,charmClusterMergingMinMass(3.5)
                      ,bottomClusterMergingMinMass(5.7)
                      ,minClusterMass(0.26996)
                      ,maxClusterMass(10.)
                      ,minEnergyDensity(0.2)
                      ,maxEnergyDensity(0.5)
                      ,minGammaS(0.5)
                      ,maxGammaS(1.)
                      ,samplingEnergyDensity(minEnergyDensity)
                      ,samplingTemperature(0.160)
                      ,hadronDataSetFileName("")
                      ,partitionFunctionDataSetPath("")
                      ,lightHadronMaxMass(1.8)
                      ,randomNumberGeneratorSeed(9019984)
    {
    }
    
    
    // Microcanonical model free parameters

    /**
    * Hadronization energy density
    */
	double energyDensity;
    
    /**
    * Strangeness suppression parameter
    */
    double gammaS;
    
    
    
    // Cluster merging parameters
    
    /**
    * Cluster merging activation flag 
    */ 
	bool clusterMergingFlag;
    
    /**
    * Cluster merging mass threshold (for light and heavy flavored clusters)
    */ 
    double clusterMergingMinMass;

    /**
    * Charmed cluster merging mass threshold
    */
    double charmClusterMergingMinMass;

    /**
    * Bottomed cluster merging mass threshold
    */
    double bottomClusterMergingMinMass;
    


    // Partition function interpolation parameters
    // TODO: to be used together with the final handling of partition
    // function data structure and likely moved to another place
    /**
    * Partition function interpolation minimum cluster mass
    */
    double minClusterMass;

    /**
    * Partition function interpolation maximum cluster mass
    */
    double maxClusterMass;
    
    /**
    * Partition function interpolation minimum cluster energy density
    */
    double minEnergyDensity;

    /**
    * Partition function interpolation maximum cluster energy density
    */
    double maxEnergyDensity;

    /**
    * Partition function interpolation minimum strangeness suppression parameter
    */
    double minGammaS;

    /**
    * Partition function interpolation maximum strangeness suppression parameter
    */
    double maxGammaS;
    

    
    // Hadronization channel sampling parameters
    
    /**
     * Hadronization channel sampling energy density parameter
     */
    double samplingEnergyDensity;
    
    /**
     * Hadronization channel sampling temperature parameter
     */
    double samplingTemperature;
    

    
    // Resource data info
    
    /**
    * Hadron list file (with absolute path)
    */
    string hadronDataSetFileName;
    
    /**
    * Partition function grids path
    */ 
    string partitionFunctionDataSetPath;
    

    
    // Hadron selection parameters
    
    /**
    * Maximum light flavored hadron mass
    */
    double lightHadronMaxMass;
    
    
    
    // Random number generator parameters
    
    /**
    * Random number generator initialization seed
    */
	int randomNumberGeneratorSeed;

};

#endif