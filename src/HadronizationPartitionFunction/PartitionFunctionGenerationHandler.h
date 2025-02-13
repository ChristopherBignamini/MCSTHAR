#ifndef PARTITIONFUNCTIONGENERATIONHANDLER_H// TODO: check this practice
#define PARTITIONFUNCTIONGENERATIONHANDLER_H

#include "../../HadronizationHadronSet/Include/HadronSamplingGroups.h"
#include "../../../Utilities/Include/ChargeConfiguration.h"
#include <vector>
#include <string>

using namespace std;

/**
* @brief Microcanonical partition function calculation handling base class
*
* This class is used to handle the partition function calculation, 
* in particular through the PartitionFunctionGenerator (base, slave and master)
* classes. This class provides the base interface and data members for the derived 
* (MPI related) classes PartitionFunctionGenerationHandlerMaster and 
* PartitionFunctionGenerationHandlerSlave
*
* @author Christopher Bignamini
*/
class PartitionFunctionGenerationHandler
{
    public:
        
        /**
        * @brief Constructor
        *
        * @param i_hadronDataSetFileName Hadron list file
        * @param i_lightHadronMaxMass Maximum light flavored hadron mass
        * @param i_minEnergyDensity Minimum energy density parameter value
        * @param i_maxEnergyDensity Maximum energy density parameter value
        * @param i_numberOfEnergyDensityValues Number of values in energy density grid 
        * @param i_massValueList Mass values list
        * @param i_chargeConfiguration Abelian charge configuration
        * @throw HadronizationException in case of error in provided setup (see also data member constructor exceptions)
        */
        PartitionFunctionGenerationHandler(const string& i_hadronDataSetFileName,
                                           double i_lightHadronMaxMass,
                                           double i_minEnergyDensity,
                                           double i_maxEnergyDensity,
                                           double i_numberOfEnergyDensityValues,
                                           const vector<double>& i_massValueList,
                                           const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration);
    
        /**
        * @brief Destructor
        */
        virtual ~PartitionFunctionGenerationHandler(void);// TODO: remove empty destructor
        
        /**
        * @brief Generation run method
        *
        * A call to this method starts the generation procedure
        *
        * @throw HadronizationException in case of error during generation procedure or handler initialization
        */
        virtual void run(void) = 0;
        
    protected:
        
        /**
        * @brief Private copy constructor (uncopiable class) // TODO: check this!
        */
        PartitionFunctionGenerationHandler(const PartitionFunctionGenerationHandler&);
        
        /**
        * @brief Private assignement (uncopiable class) // TODO: check this!
        */
        PartitionFunctionGenerationHandler& operator=(const PartitionFunctionGenerationHandler&);
        
        /**
        * Handler initialization status flag
        */
        bool m_isHandlerAvailable;
        
        /**
        * Hadronization channel sampling group handler
        */
        HadronSamplingGroups m_hadronSamplingGroups;
        
        /**
        * Partition function mass value set
        */
        const vector<double> m_massValues;
        
        /**
        * Partition function energy density value set
        */
        vector<double> m_energyDensityValues;
        
        /**
        * Partition function data set charge configuration
        */
        const MCSTHAR::Utilities::ChargeConfiguration m_chargeConfiguration;    
};

#endif
