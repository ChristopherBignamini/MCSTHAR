#ifndef PHASESPACEAPPROXIMATEDINTEGRANDFUNCTION_H
#define PHASESPACEAPPROXIMATEDINTEGRANDFUNCTION_H

#include "PhaseSpaceIntegrandFunction.h"

/**
* @brief Phase space integrand function calculation auxiliary function
*
* This class implements the integrand (sin product) function of the one dimensional integral
* of Eq. (//TODO: add reference) used for phase space integrand function
* approximated calculation, used in particular in case of large phase
* space dimensions (number of particles).
*
* WARNING: this class can be used only for set of more then four momenta 
*
* @author Christopher Bignamini
*/
class SinusoidalPhaseSpaceFunction
{
    public:
        
        /**
        * @brief Constructor
        *
        * @param i_particleMomenta Particle momenta
        * @param i_mu Particle momenta scaling factor (see reference quoted in class description)
        * @throw HadronizationException if less than four particle momenta are being provided or in case of non positive i_mu factor
        */
        SinusoidalPhaseSpaceFunction(const vector<double> i_particleMomenta,
                                     double i_mu);
        
        /**
        * @brief Destructor
        */
        ~SinusoidalPhaseSpaceFunction(void);
        
        // Default copy constructor and assignement operator are being used
        
        /**
        * @brief Function evaluation method
        *
        * This method returns the function value corresponding to the provided
        * evaluation point if i_x in [0,1]. Zero is returned otherwise.
        *
        * @param i_x Function variable value
        * @param Function value corresponding to the provided evaluation point
        */
        double operator()(double i_x) const;
        
    private:
        
        /**
        * Number of momenta
        */
        const unsigned int m_numberOfMomenta;
        
        /**
        * Number of momenta - 2 (useful for optimization)
        */
        const int m_numberOfMomentaLessTwo;
        
        /**
        * Number of momenta - 4 (useful for optimization)
        */
        const int m_numberOfMomentaLessFour;
        
        /**
        * Momenta scaling factor i_mu^(m_numberOfMomenta-3)/(2pi^2)
        */
        const double m_muFactor;
        
        /**
        * Scaled momenta (see class description)
        */
        vector<double> m_particleScaledMomenta;
        
        // TODO: set final value
        /**
        * Function variable range limits tolerance
        */
        static const double m_xBoundaryTolerance;
};

/**
* @brief Phase space approximated integrand function functor
*
* This class is the derived approximated phase space integrand function
* calculation functor. This functor performs the approximated integrand
* function (trapezoid rule based)calculation as described in // TODO: add reference Eq..
* The provided output corresponds to the integrand function value
* except for a factorized factor common to the approximated and exact
* calculation procedure.
*
* @author Christopher Bignamini
*/
class PhaseSpaceApproximatedIntegrandFunction : public PhaseSpaceIntegrandFunction
{
    public:
    
        /**
        * Constructor
        */
        PhaseSpaceApproximatedIntegrandFunction(void);
    
        /**
        * Destructor
        */
        ~PhaseSpaceApproximatedIntegrandFunction(void);
        
        // Default copy constructor and assignement operator are being used
        
        /**
        * @brief Function evaluation method
        *
        * @param i_particleMomenta Particle momenta
        * @param io_integrandFunctionValue Function value corresponding to the provided particle momenta (I/O parameter)
        * @return Calculation error status: true if error occured during calculation (unstable integral), false otherwise 
        */
        bool operator()(const vector<double>& i_particleMomenta,
                        double& io_integrandFunctionValue);
    
    private:
    
        /**
        * Trapezoid rule based integration stability threshold
        */
        // TODO: static or not? Add input parameter?
        static const double m_trapezoidalIntegralErrorThreshold;
        
        /**
        * Trapezoid rule based integration maximum number of iterations
        */
        // TODO: static or not? Add input parameter?
        static const unsigned int m_trapezoidalIntegralMaxNumberOfIterations;
        
        /**
        * Minimum number of trapezoid rule based integration iterations performed before first stability check
        */
        // TODO: static or not? Add input parameter?
        static const unsigned int m_trapezoidalIntegralMinimumStabilityCheckIterations;

};

#endif