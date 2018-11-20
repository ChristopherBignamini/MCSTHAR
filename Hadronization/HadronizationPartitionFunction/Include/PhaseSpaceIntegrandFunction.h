#ifndef PHASESPACEINTEGRANDFUNCTION_H
#define PHASESPACEINTEGRANDFUNCTION_H

#include <vector>

using namespace std;

/**
* @brief Phase space integrand function base functor 
*
* This class is the base functor for phase space integrand 
* function calculation, whose role in phase space integration is
* described in //TODO: add ref 
*
* @author Christopher Bignamini
*/
class PhaseSpaceIntegrandFunction
{
    public:
        
        /**
        * @brief Constructor
        */
        PhaseSpaceIntegrandFunction(void);
        
        /**
        * @brief Destructor
        */
        virtual ~PhaseSpaceIntegrandFunction(void);
        
        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Function evaluation method
        *
        * @param i_particleMomenta Particle momenta
        * @param io_integrandFunctionValue Function value corresponding to the provided particle momenta (I/O parameter)
        * @return Calculation error status: true if error occured during calculation, false otherwise
        */
        virtual bool operator()(const vector<double>& i_particleMomenta,
                                double& io_integrandFunctionValue) = 0;
    
    protected:
    
        /**
        * Number of particles momenta
        */
        unsigned int m_numberOfParticles;

        /**
        * Number of particles - 3 (stored for faster calculation)
        */
        int m_numberOfParticlesLessThree;

        /**
        * Number of particle momenta involved in previous function evaluation
        */
        unsigned int m_referenceNumberOfParticles;
};

#endif