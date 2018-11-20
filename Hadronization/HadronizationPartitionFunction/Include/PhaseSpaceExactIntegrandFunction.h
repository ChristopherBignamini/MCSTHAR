#ifndef PHASESPACEEXACTINTEGRANDFUNCTION_H
#define PHASESPACEEXACTINTEGRANDFUNCTION_H

#include "PhaseSpaceIntegrandFunction.h"

/**
* @brief Phase space exact integrand function functor
*
* This class is the derived exact phase space integrand function
* calculaiton functor. This functor performs the exact integrand
* function calculation as described in // TODO: add reference Eq..
* The provided output corresponds to the integrand function value
* except for a factorized factor common to the approximated and exact
* calculation procedure.
*
* @author Christopher Bignamini
*/
class PhaseSpaceExactIntegrandFunction : public PhaseSpaceIntegrandFunction
{
    public:
    
        /**
        * @brief Constructor
        */
        PhaseSpaceExactIntegrandFunction(void);
    
        /**
        * @brief Destructor
        */
        ~PhaseSpaceExactIntegrandFunction(void);

        // Default copy constructor and assignement operator are being used
        
        /**
        * @brief Function evaluation method
        *
        * @param i_particleMomenta Particle momenta
        * @param io_integrandFunctionValue Function value corresponding to the provided particle momenta (I/O parameter)
        * @return Calculation error status: true if error occured during calculation (not forseen at present), false otherwise// TODO: no error can occur? Is there a better way to unify the two function interfaces?
        */
        bool operator()(const vector<double>& i_particleMomenta,
                        double& io_integrandFunctionValue);
    
    private:

        /**
        * @brief Exact integrand function calculation constants/data structures update method
        *
        * This method performs the update of all the data and structures only depending on the number 
        * of particle and that can be therefore reused in in a new funtion evaluation call. 
        * The update function is called only in case of a change in number of particles.
        */
        void updateCalculationData(void);
                
        /**
        * Maximum number of sign configurations for current set of particle momenta (2^m_numberOfParticles)
        */
        unsigned int m_maxNumberOfConfiguration;
        
        /**
        * Factor (-1/2^(N+1)*pi*(N-3)!) with N=m_numberOfParticles (see reference quoted in class description)
        */
        double m_exactIntegrandFunctionFactor;
        
        /**
        * Indexes corresponding to the position of the next sign change in each sign configuration
        */
        vector<unsigned int> m_configurationsNewMinusPositionIndex;
        
        /**
        * Particle momenta sums corresponding to different sign configurations
        */
        vector<double> m_particleMomentaSums;
        
        /**
        * Single sign configurations total sign
        */
        vector<int> m_configurationTotalSign;
        
        /**
        * Provided particle momenta sorted in increasing order
        */
        vector<double> m_sortedParticleMomenta;
        
        /**
        * Provided particle momenta sorted in increasing order times 2
        */
        vector<double> m_sortedParticleMomentaTimesTwo;    
};


#endif