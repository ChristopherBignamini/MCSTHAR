#ifndef PHASESPACEINTEGRATIONDATA_H
#define PHASESPACEINTEGRATIONDATA_H

/**
* @brief Phase space integration data storage structure
*
* This structure is used to store the data corresponding
* to a single hadronization channel pahse space integration.
*
* @author Christopher Bignamini
*/
struct PhaseSpaceIntegrationData
{
    /**
    * Phase space integral value
    */
    double phaseSpaceIntegral;
    
    /**
    * Phase space integral error value
    */
    double phaseSpaceIntegralError;
    
    /**
    * Phase space integration error status (true if below the defined threshold, false otherwise)
    */
    bool isErrorUnderThreshold;
    
    /**
    * Number of phase space configuration sampling performed for the integration procedure
    */
    unsigned numberOfSamplings;
};

#endif
