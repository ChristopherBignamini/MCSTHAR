#ifndef PARTITIONFUNCTIONCALCULATIONSETUPLOADER_H
#define PARTITIONFUNCTIONCALCULATIONSETUPLOADER_H

#include "PartitionFunctionCalculationSetup.h"
#include <string>

using namespace std;

/**
* @brief Partition function calculation setup load class
*
* @author Christopher Bignamini
*/
class PartitionFunctionCalculationSetupLoader
{
    public:
    
        /**
        * @brief Constructor
        *
        * @param i_partitionFunctionCalculationSetupFileName Partition function calculation setup file
        * @throw HadronizationException in case of error during setup file parsing
        */
        PartitionFunctionCalculationSetupLoader(const string& i_partitionFunctionCalculationSetupFileName);
    
        /**
        * @brief Destructor
        */
        ~PartitionFunctionCalculationSetupLoader(void);
    
        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Partition function calculation setup data get method
        *
        * @return Partition function calculation setup data
        * @throw HadronizationException in case of partition function calculation setup data unavailability
        */
        const PartitionFunctionCalculationSetup& getPartitionFunctionCalculationSetup(void) const;
    
    private:
    
        /**
        * Partition function calculation setup data availability status
        */
        bool m_partitionFunctionCalculationSetupAvailability;
    
        /**
        * Partition function calculation setup data
        */
        PartitionFunctionCalculationSetup m_partitionFunctionCalculationSetup;
};

#endif
