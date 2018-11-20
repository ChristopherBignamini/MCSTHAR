#include "../Include/PhaseSpaceApproximatedIntegrandFunction.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/Constants.h"
#include "../../../Utilities/Include/logMessage.h"
#include <limits>
#include <sstream>

// TODO: debug
#include <iomanip>
using namespace std;

const double SinusoidalPhaseSpaceFunction::m_xBoundaryTolerance(1.e-8);// TODO: set final value

// TODO: check if better design is possible: in relation to the integration procedure, it would
// be useful to define a function base class which can be used with integration module
// making use of polymorphism
// TODO: remove static cast
SinusoidalPhaseSpaceFunction::SinusoidalPhaseSpaceFunction(const vector<double> i_particleMomenta,
                                                           const double i_mu)
                                                          :m_numberOfMomenta(i_particleMomenta.size())
                                                          ,m_numberOfMomentaLessTwo(static_cast<int>(m_numberOfMomenta-2))
                                                          ,m_numberOfMomentaLessFour(static_cast<int>(m_numberOfMomenta-4))
                                                          ,m_muFactor(pow(i_mu,static_cast<int>(m_numberOfMomenta-3))/twoPiSquared)
                                                          ,m_particleScaledMomenta(m_numberOfMomenta)
{
    // TODO: add missing check on class definition parameters (e.g., what happens with 2 particles?)
    if(m_numberOfMomenta<4)
    {
        throw HadronizationException("Error during phase space integrand function calculation, less than 4 particle momenta provided",
                                     __FUNCTION__,841);
    }
    // TODO: mu should be automatically calculated by the present functor class
    if(i_mu<=0.0)
    {
        throw HadronizationException("Error during phase space integrand function calculation, non positive particle scaling factor provided",
                                     __FUNCTION__,842);
    }
    
    for(unsigned int momentaIndex=0;momentaIndex<m_numberOfMomenta;++momentaIndex)
    {
        m_particleScaledMomenta[momentaIndex] = i_particleMomenta[momentaIndex]/i_mu;
    }
    
}

SinusoidalPhaseSpaceFunction::~SinusoidalPhaseSpaceFunction(void)
{
}

double SinusoidalPhaseSpaceFunction::operator()(const double i_x) const
{
    // TODO: add check on x value, in particular to avoid x = 0 and x = 1
    const double oneLessX(1.-i_x);
    if((oneLessX>m_xBoundaryTolerance) && (i_x>m_xBoundaryTolerance))
    {
        const double xOnOneLessX(i_x/oneLessX);
        // TODO: remove static cast
        double o_functionValue(pow(oneLessX,m_numberOfMomentaLessFour)/pow(i_x,m_numberOfMomentaLessTwo));
        for(unsigned int momentaIndex=0;momentaIndex<m_numberOfMomenta;++momentaIndex)
        {
            // TODO: hot spot here on sin calculation, evaluate possibility to switch to look up table/trigonometric formulas, etc..
            o_functionValue *= sin(m_particleScaledMomenta[momentaIndex]*xOnOneLessX);
        }
        return o_functionValue;
    }
    else
    {
        return 0.;
    }
}

const double PhaseSpaceApproximatedIntegrandFunction::m_trapezoidalIntegralErrorThreshold(1.e-2);
const unsigned int PhaseSpaceApproximatedIntegrandFunction::m_trapezoidalIntegralMaxNumberOfIterations(20);
const unsigned int PhaseSpaceApproximatedIntegrandFunction::m_trapezoidalIntegralMinimumStabilityCheckIterations(5);

PhaseSpaceApproximatedIntegrandFunction::PhaseSpaceApproximatedIntegrandFunction(void)
{
}

PhaseSpaceApproximatedIntegrandFunction::~PhaseSpaceApproximatedIntegrandFunction(void)
{
}

bool PhaseSpaceApproximatedIntegrandFunction::operator()(const vector<double>& i_particleMomenta,
                                                         double& io_integrandFunctionValue)
{
    m_numberOfParticles = i_particleMomenta.size();
    
    // Check if update is required to constant data due to change in number of momenta
    if(m_numberOfParticles!=m_referenceNumberOfParticles)
    {
        m_referenceNumberOfParticles = m_numberOfParticles;
        m_numberOfParticlesLessThree = m_numberOfParticles-3;
    }

    // Build sinusoidal auxiliary functor
    // Set mu value (see reference quoted in //TODO:where???)// TODO: move to method for mu calculation?
    double mu(1.);
    unsigned int numberOfMomentaAboveThreshold(0);
    const double momentaMuThreshold(0.1);// TODO: use class parameter
    // TODO: move loop within sinusoidal function calculation for optimization
    for(unsigned int momentaIndex=0;momentaIndex<m_numberOfParticles;++momentaIndex)
    {
        if(i_particleMomenta[momentaIndex]>momentaMuThreshold)
        {
            ++numberOfMomentaAboveThreshold;
        }
    }
    const unsigned int numberOfParticleMuThreshold(10);// TODO: use class parameter
    if(numberOfMomentaAboveThreshold<=numberOfParticleMuThreshold)
    {
        mu = 0.1;
    }
    
    SinusoidalPhaseSpaceFunction sinusoidalFunction(i_particleMomenta,mu);// TODO: optimize by saving construction at each new calculation
    
    // Run integration
    io_integrandFunctionValue = 0.;
    double integrandFunctionError(0.);
    bool isIntegrationStable(false);
    // Loop over trapezoid rule based integration iterations
    // TODO: this call could be avoided by assuming a known value for the integrated function (which is always zero for N>=4)
    try
    {
        io_integrandFunctionValue = (sinusoidalFunction(1.) + sinusoidalFunction(0.))*0.5;//TODO: set to zero without calculation!!!!
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    integrandFunctionError = io_integrandFunctionValue;
    
    for(unsigned int iterationNumber=1;iterationNumber<=m_trapezoidalIntegralMaxNumberOfIterations;++iterationNumber)
    {
        // TODO: pows could be partially avoided by using a class for trapezoidal integration
        // TODO: division by 2 can be transformed in a bitwise operation
        // TODO: avoid static cast
        // TODO: use 0.5 powers where possible to avoid divisions
        const double deltaX(1./pow(2.,static_cast<int>(iterationNumber-1)));
        // TODO: avoid static cast
        // TODO: avoid 2.?
        const unsigned int numberOfIntegrationPoints(static_cast<unsigned int>(pow(2.,static_cast<int>(iterationNumber-1))));
        double integrationPoint(0.5*deltaX);
        
        try
        {
            double newIntegralContribution(sinusoidalFunction(integrationPoint));
            for(unsigned int integrationPointIndex=1;integrationPointIndex<numberOfIntegrationPoints;++integrationPointIndex)
            {
                integrationPoint += deltaX;
                newIntegralContribution += sinusoidalFunction(integrationPoint);
            }
            newIntegralContribution /= numberOfIntegrationPoints;
            integrandFunctionError = 0.5*(newIntegralContribution - io_integrandFunctionValue);
            io_integrandFunctionValue += newIntegralContribution;
            io_integrandFunctionValue *= 0.5;
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
                
        // Check integration stability after m_trapezoidalIntegralMinimumStabilityCheckIterations iterations
        if(iterationNumber>m_trapezoidalIntegralMinimumStabilityCheckIterations)
        {            
            const double absIntegrandFunctionError(std::abs(integrandFunctionError));
            const double absIntegrandFunction(std::abs(io_integrandFunctionValue));
            
            // Check for very small integrand function and error values
//            if(absIntegrandFunctionError<=numeric_limits<double>::epsilon())// TODO: larger values for threshold should be considered
            if(absIntegrandFunctionError<1.e-15 || absIntegrandFunction<1.e-13)// TODO: larger values for threshold should be considered
            {
                integrandFunctionError = 0.;
                isIntegrationStable = true;
                break;
            }
            
            // TODO: check if better stability evaluation condition exist
            if(absIntegrandFunctionError<m_trapezoidalIntegralErrorThreshold*absIntegrandFunction)
            {
                isIntegrationStable = true;
                break;
            }
        }
    }
    
    // Check stability status and return result
    if(isIntegrationStable)
    {
        // Compute missing factor
        // TODO: factorize, this is constant for al fixed channel integration
        io_integrandFunctionValue *= pow(mu,static_cast<int>(m_numberOfParticles-3))/twoPiSquared;
    }
    else
    {
        #ifdef _FULLLOG // TODO: add the possibility to check print level from setup file
            // Unstable integral, log warning
            string warnMessage("Error during approximate integrand function calculation, unstable integral. Particle momenta:\n");
            for(unsigned int momentumIndex=0;momentumIndex<m_numberOfParticles;++momentumIndex)
            {            
                stringstream momentumStr;
                momentumStr<<setprecision(17);// TODO: switch to appropriate double precision
                momentumStr<<i_particleMomenta[momentumIndex];
                warnMessage += "\n" + momentumStr.str();
            }
            warnMessage += "\n";
            logMessage(warnMessage,WARNING);// TODO: handling for multithreaded logging required
        #endif
        
        return true;
    }
    return false;
}

