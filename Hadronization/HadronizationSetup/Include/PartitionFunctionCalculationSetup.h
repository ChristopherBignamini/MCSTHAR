#ifndef PARTITIONFUNCTIONCALCULATIONSETUP_H
#define PARTITIONFUNCTIONCALCULATIONSETUP_H

#include "MicrocanonicalParameterGrid.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include <vector>
#include <string>

using namespace std;

// TODO: switch to structure and add the missing checks
// TODO: factorization with xml parsing structure required!!!
// TODO: factorization with hadronization setup file
// TODO: add sensible parameter to xml part. funct data (e.g., hadron data file, min mass, etc..)
/**
* @brief Partition function calculation setup parameter storage structure
*
* @author Christopher Bignamini
*/
struct PartitionFunctionCalculationSetup
{
    /**
    * @brief Constructor
    *
    */
    // TODO: move default parameter somewhere else
    PartitionFunctionCalculationSetup(void)
                                     :hadronDataSetFileName("")
                                     ,lightHadronMaxMass(1.8)
                                     ,chargeConfiguration(0,
                                                          0,
                                                          0,
                                                          0.,
                                                          0.)
                                     ,microcanonicalParameterGrid(0.35,
                                                                  0.35,
                                                                  1,
                                                                  0.65,
                                                                  0.65,
                                                                  1)
                                     ,massValueList(vector<double>(1,3.0))
                                     ,maxChannelSamplingNumber(100)
                                     ,integrationErrorThreshold(1.e-3)
                                     ,randomNumberGeneratorSeed(9019984)
                                     ,outputFolder("")
    {
    }
    
    // Default copy constructor and overloaded assignement operator are being used
    
    /**
    * Hadron list file full name
    */
    string hadronDataSetFileName;

    /**
    * Maximum light flavored hadron mass
    */
    double lightHadronMaxMass;

    /**
    * Charge configuration
    */
    MCSTHAR::Utilities::ChargeConfiguration chargeConfiguration;
        
    /**
    * Microcanonical parameter grid structure
    */
	MCSTHAR::Setup::MicrocanonicalParameterGrid microcanonicalParameterGrid;
    
    /**
    * Cluster mass values
    */
    vector<double> massValueList;
    
    /**
    * Maximum number of channel samplings performed for partition function calculation
    */
    unsigned long long int maxChannelSamplingNumber;
    
    /**
    * Stability threshold for partition function calculation
    */
    double integrationErrorThreshold;
        
    /**
     * Random number generator initialization seed
     */
	int randomNumberGeneratorSeed;
    
    /**
    * Partition function calculation output file storage folder
    */
	string outputFolder;
    
};

#endif
