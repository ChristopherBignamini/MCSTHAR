#ifndef PARTITIONFUNCTIONPHASESPACESAMPLING_H
#define PARTITIONFUNCTIONPHASESPACESAMPLING_H

#include "PhaseSpaceIntegrandFunction.h"
#include "PhaseSpaceExactIntegrandFunction.h"
#include "PhaseSpaceApproximatedIntegrandFunction.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include <vector>

using namespace std;

/**
* @brief Phase space sampling class
*
* This class is used for phase space configuration sampling
* during partition function calculation. For N-body configuration sampling,
* with N>2, two algorithm are implemented and used for sampled configuration
* weight calculation, whose choice is determined by the number of the considered phase space
* dimensions. The implemented strategy and algorithm is described in detail
* in // TODO: add reference. It must be noted that only the weight corresponding
* to the sampled configuration is returned by the present class, the momenta set
* not being needed for the partition function calculation. 
*
* @author Christopher Bignamini
*/
class PartitionFunctionPhaseSpaceSampling
{
    public:
    
        /**
        * @brief Constructor
        *
        * @param io_randomGenerator Random number generator
        */
        PartitionFunctionPhaseSpaceSampling(RandomNumberGenerator& io_randomGenerator);

        /**
        * @brief Destructor
        */
        ~PartitionFunctionPhaseSpaceSampling(void);
        
        /**
        * @brief Phase space configuration sampling method
        *
        * This method performs the sampling of a phase space configuration 
        * for the provided available energy and particle masses.
        *
        * @param i_availableEnergy Available phase space energy
        * @param i_particleMasses Particle masses
        * @return Sampled configuration weight 
        * @throw HadronizationException if less than 2 particle masses are provided 
        * or error during sampling procedures (e.g. no phase space availability)
        */
        double run(double i_availableEnergy,
                   const vector<double>& i_particleMasses);
                
    private:
    
        /**
        * @brief Private copy constructor (non copiable class)
        *
        * @param i_phaseSpaceSampling Phase space sampling class to be copied
        */
        PartitionFunctionPhaseSpaceSampling(const PartitionFunctionPhaseSpaceSampling& i_phaseSpaceSampling);
    
        /**
        * @brief Private assignement opertor (non copiable class)
        *
        * @param i_phaseSpaceSampling Phase space sampling class to be copied
        */
        PartitionFunctionPhaseSpaceSampling& operator=(const PartitionFunctionPhaseSpaceSampling& i_phaseSpaceSampling);
    
        /**
        * @brief Two body phase space configuration sampling method
        *
        * @param i_availableEnergy Available phase space energy
        * @param i_particleMasses Particle masses
        * @return Sampled configuration weight
        * @throw HadronizationException in case of no phase space availability
        */
        double runTwoBodySampling(double i_availableEnergy,
                                  const vector<double>& i_particleMasses) const;
    
        /**
        * @brief N (>2) body phase space configuration sampling method
        *
        * This method implements the N body phase space configuration sampling
        * procedure described in //TODO: add reference
        *
        * @param i_availableEnergy Available phase space energy
        * @param i_particleMasses Particle masses
        * @return Sampled configuration weight
        * @throw HadronizationException in case of no phase space availability or
        *   error during the sampling procedure
        */
        double runNBodySampling(double i_availableEnergy,
                                const vector<double>& i_particleMasses);
    
        /**
        * Number of particles involved in phase space sampling
        */
        unsigned int m_numberOfParticles;
                
        /**
        * Maximum number of particles allowed for phase space integrand function exact calculation (see class description)
        */
        // TODO: static const or not? Add input parameter?
        static const unsigned int m_exactIntegrandCalculationMaxMultiplicity;

        /**
        * N-body phase space integrand function in use in a single phase space configuration sampling
        * event for a given momenta set.
        * This object points during sampling to one of the two functors data members m_phaseSpaceExactIntegrandFunction/
        * m_phaseSpaceApproximatedIntegrandFunction implementing the two different integrand function calculation
        * quoted in class description.
        */
        PhaseSpaceIntegrandFunction* m_phaseSpaceIntegrandFunction;
        
        /**
        * Exact integrand function calculation functor
        */
        PhaseSpaceExactIntegrandFunction* m_phaseSpaceExactIntegrandFunction;
        
        /**
        * Approximated integrand function calculation functor
        */
        PhaseSpaceApproximatedIntegrandFunction* m_phaseSpaceApproximatedIntegrandFunction;
    
        /**
        * Random number generator
        */
        RandomNumberGenerator* m_randomGenerator;
};
#endif