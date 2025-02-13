#ifndef PARTITIONFUNCTIONGENERATIONDATA_H
#define PARTITIONFUNCTIONGENERATIONDATA_H

// TODO: header file renaming with serialization and parsing classes
#include <vector>

using namespace std;

/**
* @brief Grand canonical ensemble distribution parameters
*
* @author Christopher Bignamini
*/
struct GrandCanonicalParameters
{
    /**
    * @brief Constructor
    */    
    GrandCanonicalParameters(void)
                            :electricChargeFugacity(1.)
                            ,strangeChargeFugacity(1.)
                            ,baryonicChargeFugacity(1.)
                            ,charmChargeFugacity(1.)
                            ,bottomChargeFugacity(1.)
                            ,temperature(0.160)
    {
    }

    
    /**
    * @brief Constructor
    *
    * @param i_electricChargeFugacity Electric charge fugacity valu
    * @param i_strangeChargeFugacity Strange charge fugacity value
    * @param i_baryonicChargeFugacity Baryonic charge fugacity value
    * @param i_charmChargeFugacity Charm charge fugacity value
    * @param i_bottomChargeFugacity Bottom charge fugacity value
    * @param i_temperature Temperature value (GeV)
    */
    GrandCanonicalParameters(const double i_electricChargeFugacity,// TODO: add parameter check and switch to class
                             const double i_strangeChargeFugacity,
                             const double i_baryonicChargeFugacity,
                             const double i_charmChargeFugacity,
                             const double i_bottomChargeFugacity,
                             const double i_temperature)
                            :electricChargeFugacity(i_electricChargeFugacity)
                            ,strangeChargeFugacity(i_strangeChargeFugacity)
                            ,baryonicChargeFugacity(i_baryonicChargeFugacity)
                            ,charmChargeFugacity(i_charmChargeFugacity)
                            ,bottomChargeFugacity(i_bottomChargeFugacity)
                            ,temperature(i_temperature)
    {
    }

    // Default copy constructor and assignement operator are being used
    
    /**
    * Electric charge fugacity value
    */
    double electricChargeFugacity;
    
    /**
    * Strange charge fugacity value
    */
    double strangeChargeFugacity;
    
    /**
    * Baryonic charge fugacity value
    */
    double baryonicChargeFugacity;
    
    /**
    * Charm charge fugacity value
    */
    double charmChargeFugacity;
    
    /**
    * Bottom charge fugacity value
    */
    double bottomChargeFugacity;
    
    /**
    * Temperature value (GeV)
    */
    double temperature;
};


/**
* @brief Single mass/energy density microcanonical partition function data storage structure
*
* This structure is used to store the partition function data
* corresponding to a single mass value, energy density value and 
* charge configuration for a set of strangeness suppression parameter grid points.
*
* @author Christopher Bignamini
*/
// TODO: find final name
// TODO: switch to class
struct SingleMassEnergyDensityPairPartitionFunctionData
{
    /**
    * @brief Constructor
    *
    * @param i_numberOfGammaSValues Number of strangeness suppression parameter values
    */
    // TODO: add check for mass value set
    SingleMassEnergyDensityPairPartitionFunctionData(const unsigned int i_numberOfGammaSValues)
                                                    :partitionFunctions(i_numberOfGammaSValues,0.)
                                                    ,partitionFunctionErrors(i_numberOfGammaSValues,0.)
                                                    ,partitionFunctionErrorStatus(i_numberOfGammaSValues,false)
                                                    ,numberOfChannelSamplings(i_numberOfGammaSValues,0)
    {
    }

    // Default copy constructor and assignement operator are being used

    /**
    * Partition function values
    */
    vector<double> partitionFunctions;
    
    /**
    * Partition function error values
    */
    vector<double> partitionFunctionErrors;
    
    /**
    * Partition function error status (true if below the defined threshols, false otherwise)
    */
    vector<bool> partitionFunctionErrorStatus;
    
    /**
    * Number of channel samplings performed for the single grid points
    */
    vector<unsigned long long int> numberOfChannelSamplings;
    
    /**
    * Grand canonical ensemble parameters used for partition function set calculation (particularly for sampling function definition)
    */
    GrandCanonicalParameters grandCanonicalParameters;
};

#endif