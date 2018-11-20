#ifndef UPDATEEVENTRECORD_H
#define UPDATEEVENTRECORD_H

#include "../../../Hadronization/HadronizationSteps/Include/HadronizationEventRecord.h"

/**
* @brief Herwig6510 event record update function
*
* This function performs the update of Herwig6510 event record 
* (HEPEVT and HWEVNT common data blocks) with the insertion of
* the objects created during the hadronization procedure, namely
* clusters obtained from the merging procedure and hadrons.
*
* // TODO: Any missing check/exception?
* @param Hadronization event record
*
* @author Christopher Bignamini
*/	

void updateEventRecord(const HadronizationEventRecord& i_eventRecord);

#endif
