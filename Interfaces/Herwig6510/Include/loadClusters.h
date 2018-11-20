#ifndef LOADCLUSTERS_H
#define LOADCLUSTERS_H

#include "../../../Hadronization/HadronizationObjects/Include/Cluster.h"
#include <vector>

using namespace std;

/**
* @brief Cluster retrieval from Herwig6510 event record function
*
* This function is used to retrieve Herwig6510 cluster data from
* the available event record. The cluster energy density is provided
* as a function parameter.
*
* @param i_energyDensity Cluster energy density
* @return Herwig6510 clusters
* @throw HadronizationException in case of error during cluster data retrieval (e.g., unidentified cluster constituents)
*
* @author Christopher Bignamini 
*/
vector<Cluster> loadClusters(double i_energyDensity);

/**
* @brief Herwig6510 cluster charge retrieval function
*
*
* This function is used to retrieve Herwig6510 cluster charge
* configuration from the corresponding constituent quarks/diquarks,
* whose positions in the considered event record are provided
* by the first/second parent index parameters.
*
* @param i_firstParentIndex Cluster first constituent position index
* @param i_secondParentIndex Cluster second constituent position index
* @param io_strangeCharge Cluster strange charge (I/O parameter)
* @param io_charmCharge Cluster charm charge (I/O parameter)
* @param io_bottomCharge Cluster bottom charge (I/O parameter)
* @param io_electricCharge Cluster electric charge (I/O parameter)
* @param io_baryonicCharge Cluster baryonic charge (I/O parameter)
* @throw HadronizationException in case of error during cluster data retrieval (e.g., unidentified cluster constituents)
*
* @author Christopher Bignamini 
*/
void retrieveClusterCharges(int i_firstParentIndex,
                            int i_secondParentIndex,
                            int& io_strangeCharge,
                            int& io_charmCharge,
                            int& io_bottomCharge,
                            double& io_electricCharge,
                            double& io_baryonicCharge);

#endif
