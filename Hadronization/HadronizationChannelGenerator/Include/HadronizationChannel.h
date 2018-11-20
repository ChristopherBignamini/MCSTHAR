#ifndef HADRONIZATIONCHANNEL_H
#define HADRONIZATIONCHANNEL_H

#include "../../HadronizationObjects/Include/Particle.h"
#include <vector>

using namespace std;

/**
* @brief Single cluster hadronization channel data
*
* @author Christopher Bignamini
*/
struct HadronizationChannel
{
    /**
    * @brief Constructor
    *
    */
    HadronizationChannel(void)
                        :clusterIndex(0)
                        ,weight(0.)
                        ,weightVariance(0.)
    {
    }

    /**
    * @brief Constructor
    *
    * @param i_clusterIndex Hadronized cluster index in current event record
    * @param i_weight Channel weight
    * @param i_particles Channel particles
    */
    HadronizationChannel(const unsigned int i_clusterIndex,
                         const double i_weight,
                         const vector<Particle> i_particles)
                        :clusterIndex(i_clusterIndex)
                        ,weight(i_weight)
                        ,weightVariance(0.)
                        ,particles(i_particles)
    {
    }
        
    // Default copy constructor and overloaded assignement operator are being used
        
    /**
    * Hadronized cluster index in current event record
    */
    unsigned int clusterIndex;
    
    /**
    * Channel weight
    */
    double weight;
    
    /**
    * Channel weight variance (TODO: at present always set to zero)
    */
    double weightVariance;
    
    /**
    * Channel particles
    */
    vector<Particle> particles;
    
};

#endif