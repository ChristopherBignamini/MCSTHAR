#ifndef PHASESPACEINTEGRATOR_H
#define PHASESPACEINTEGRATOR_H

#include "PhaseSpaceIntegrationData.h"
#include "PhaseSpaceIntegrandFunction.h"
#include "PhaseSpaceExactIntegrandFunction.h"
#include "PhaseSpaceApproximatedIntegrandFunction.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include <vector>

using namespace std;

/**
* @brief Phase space integration class
*
* This class is used for phase space integral calculation during partition
* function calculation. Two calculation algorithm are implemented and used
* whose choice is determined by the number of the considered phase space 
* dimensions. The implemented strategy and algorithm is described in detail
* in // TODO: add reference. 
*
* @author Christopher Bignamini
*/
class PhaseSpaceIntegrator
{
    public:
    
        /**
        * @brief Constructor
        *
        * @param i_maxPhaseSpaceSamplingNumber Maximum number of allowed integration phase space samplings
        * @param i_phaseSpaceIntegrationErrorThreshold Phase space integration stability threshold (relative error)
        * @param io_randomGenerator Random number generator
        * @throw HadronizationException in case of zero maximum number of samplings or negative error threshold provided
        */
        PhaseSpaceIntegrator(unsigned int i_maxPhaseSpaceSamplingNumber,
                             double i_phaseSpaceIntegrationErrorThreshold,
                             RandomNumberGenerator& io_randomGenerator);

        /**
        * @brief Destructor
        */
        ~PhaseSpaceIntegrator(void);
        
        /**
        * @brief Phase space integral calculation method
        *
        * After a succesfull call to this method the phase space integral
        * corresponding to the provided input is available. The full data
        * can be retrieved by means of the getIntegrationData method 
        *
        * @param i_availableEnergy Available phase space energy
        * @param i_particleMasses Particle masses
        * @throw HadronizationException in case of error during class initialization, 
        *   less than 2 particle masses provided or error during integration procedures (e.g. no phase space availability)
        */
        void run(double i_availableEnergy,
                 const vector<double>& i_particleMasses);
    
        /**
        * @brief Phase space integration data retrieval method
        *
        * @return Phase space integration data corresponding to the provided integration input
        * @throw HadronizationException in case not available integration data
        */
        const PhaseSpaceIntegrationData& getIntegrationData(void) const;
    
    private:
    
        /**
        * @brief Private copy constructor
        */
        PhaseSpaceIntegrator(const PhaseSpaceIntegrator&);

        /**
        * @brief Private assignement operator
        */
        PhaseSpaceIntegrator& operator=(const PhaseSpaceIntegrator&);
    
        /**
        * @brief Two body phase space integral calculation method
        *
        * @param i_availableEnergy Available phase space energy
        * @param i_particleMasses Particle masses
        * @throw HadronizationException in case of no phase space availability
        */
        void runTwoBodyCalculation(double i_availableEnergy,
                                   const vector<double>& i_particleMasses);
    
        /**
        * @brief N (>2) body phase space integral calculation method
        *
        * This method implements the N body phase space Monte Carlo integration
        * procedure described in //TODO: add reference
        *
        *
        * @param i_availableEnergy Available phase space energy
        * @param i_particleMasses Particle masses
        * @throw HadronizationException in case of no phase space availability or 
        *   error during the integration procedure
        */
        void runNBodyCalculation(double i_availableEnergy,
                                 const vector<double>& i_particleMasses);
    
        /**
        * Class initialization status flag
        */
        bool m_isIntegratorReady;
        
        /**
        * Phase space integration data availability flag
        */
        bool m_isIntegralAvailable;

        /**
        * Number of particles involved in phase space integration
        */
        unsigned int m_numberOfParticles;
        
        /**
        * Phase space integration data
        */ 
        PhaseSpaceIntegrationData m_phaseSpaceIntegrationData;
    
        /**
        * Maximum number of allowed phase space integration samplings
        */
        unsigned int m_maxSamplingNumber;

        /**
        * Integration stability squared (relative) error threshold
        */
        const double m_integrationSquaredErrorThreshold;
    
        /**
        * Minimum number of phase space integration samplings performed before first stability check
        */
        // TODO: static const or not?  Add input parameter?
        unsigned int m_minimumStabilityCheckSamplingNumber;    

        /**
        * Maximum number of particles allowed for phase space integrand function exact calculation (see class description)
        */
        // TODO: static const or not? Add input parameter?
        static const unsigned int m_exactIntegrandCalculationMaxMultiplicity;
    
#ifdef _OPENMP

        // OpenMP based parallel calculation data members
    
        /**
        * Single thread maximum number of samplings
        */
        unsigned int m_localMaxSamplingNumber;

        /**
        * Single thread minimum number of phase space integration samplings performed before first stability check
        */
        unsigned int m_localMinStabilityCheckSamplingNumber;

        // TODO: check best practice!
        /**
        * N-body phase space integrand function in use in a single phase space integration 
        * for a given momenta set. The vector contains pointers to specialized phase space integrand functor,
        * with each vector element associated to a OpenMP thread.    
        * Each element points during integration to a corresponding element in one of the two functors vector data members 
        * m_phaseSpaceExactIntegrandFunction/m_phaseSpaceApproximatedIntegrandFunction, again with the above thread/element
        * correspondence, implementing the two different integrand function calculation quoted in class description. 
        */
        vector<PhaseSpaceIntegrandFunction*> m_phaseSpaceIntegrandFunctions;
    
        /**
        * Exact integrand function calculation functors (each vector element correspond to an OpenMP thread) 
        */
        vector<PhaseSpaceExactIntegrandFunction*> m_phaseSpaceExactIntegrandFunctions;
        
        /**
        * Approximated integrand function calculation functors (each vector element correspond to an OpenMP thread)
        */
        vector<PhaseSpaceApproximatedIntegrandFunction*> m_phaseSpaceApproximatedIntegrandFunctions;
    
        /**
        * Random number generators (each vector element correspond to an OpenMP thread)
        */
        vector<RandomNumberGenerator*> m_randomGenerators;// TODO: check for best practices

        /**
        * Number of OpenMP threads
        */
        unsigned int m_numberOfThreads;
    
#else
        // Serial calculation data members
    
        /**
        * N-body phase space integrand function in use in a single phase space integration
        * for a given momenta set.
        * This object points during integration to one of the two functors data members m_phaseSpaceExactIntegrandFunction/
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
    
#endif

};

#endif