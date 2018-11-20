#ifndef HADRONIZATIONSETUPLOADER_H
#define HADRONIZATIONSETUPLOADER_H

#include "HadronizationSetup.h"
#include <string>

using namespace std;

/**
* @brief Hadronization setup load class
*
* @author Christopher Bignamini
*/
class HadronizationSetupLoader
{
    public:
    
        /**
        * @brief Constructor
        *
        * @param i_hadronizationSetupFileName Hadronization setup file
        * @throw HadronizationException in case of error during setup file parsing
        */
        HadronizationSetupLoader(const string& i_hadronizationSetupFileName);
    
        /**
        * @brief Destructor
        *
        */
        ~HadronizationSetupLoader(void);
    
        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Hadronization setup data get method
        *
        * @return Hadronization setup data
        * @throw HadronizationException in case of hadronization setup data unavailability
        */
        const HadronizationSetup& getHadronizationSetup(void);
    
    private:
    
        /**
        * Hadronization setup data availability status
        */
        bool m_hadronizationSetupAvailability;
    
        /**
        * Hadronization setup data
        */
        HadronizationSetup m_hadronizationSetup;
};

#endif
