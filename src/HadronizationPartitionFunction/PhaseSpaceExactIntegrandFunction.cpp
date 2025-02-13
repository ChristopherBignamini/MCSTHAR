#include "../Include/PhaseSpaceExactIntegrandFunction.h"
#include "../../../Utilities/Include/Constants.h"
#include "../../../Utilities/Include/mathFunctions.h"
#include <algorithm>
#include <cassert>
#include <cmath>

PhaseSpaceExactIntegrandFunction::PhaseSpaceExactIntegrandFunction(void)
{
}

PhaseSpaceExactIntegrandFunction::~PhaseSpaceExactIntegrandFunction(void)
{
}

bool PhaseSpaceExactIntegrandFunction::operator()(const vector<double>& i_particleMomenta,
                                                  double& io_integrandFunctionValue)
{
    m_numberOfParticles = i_particleMomenta.size();
    // Check number of involved momenta versus the previous call one
    if(m_numberOfParticles!=m_referenceNumberOfParticles)
    {
        // Update integrand function constant contributions only dependent on particle number
        updateCalculationData();
        m_referenceNumberOfParticles = m_numberOfParticles;
    }
    else
    {
        // The number of momenta is the same of the previous function call, reset data
        m_configurationsNewMinusPositionIndex.assign(m_maxNumberOfConfiguration,0);
        m_configurationTotalSign.assign(m_maxNumberOfConfiguration,1.);
    }
    
    // Sort particle momenta
    m_sortedParticleMomenta = i_particleMomenta;
    sort(m_sortedParticleMomenta.begin(),m_sortedParticleMomenta.end());

    
    // Compute first combination (all positive signs)
    unsigned int storedConfigurationNumber(0);
    int newConfigurationTotalSign(1);
    double particleMomentaSumsElement(m_sortedParticleMomenta[0]);
    for(unsigned int momentaIndex=1;momentaIndex<m_numberOfParticles;++momentaIndex)
    {
        particleMomentaSumsElement += m_sortedParticleMomenta[momentaIndex];
    }
    m_particleMomentaSums[0] = particleMomentaSumsElement;
    
    assert(particleMomentaSumsElement>0.);
    io_integrandFunctionValue = 0.;
    if(particleMomentaSumsElement>0.)
    {
        // Positive momenta combination, update integrand function value.
        io_integrandFunctionValue += pow(particleMomentaSumsElement,m_numberOfParticlesLessThree);
        // Last minus position and stored configuration start index already correctly set
        ++storedConfigurationNumber;
    }
    
    // Compute sorted momenta x 2
    for(unsigned int momentaIndex=0;momentaIndex<m_numberOfParticles;++momentaIndex)
    {
        // TODO: STL function?
        m_sortedParticleMomentaTimesTwo[momentaIndex] = 2.*m_sortedParticleMomenta[momentaIndex];
    }
    
    // Loop over stored accepted configuration
    bool restorePreviousSign(false);
    for(unsigned int configurationIndex=0;configurationIndex<storedConfigurationNumber;++configurationIndex)
    {
        // Get momenta sum corresponding to the current configuration
        particleMomentaSumsElement = m_particleMomentaSums[configurationIndex];
        
        // Set new configuration total sign (a minus sign is added somewhere in the starting configuration)
        newConfigurationTotalSign = -1*m_configurationTotalSign[configurationIndex];// TODO: is not this only flipping at each iteration?
        
        // Loop over all configuration that can be obtained by adding a minus sign to the starting one
        restorePreviousSign = false;
        for(unsigned int signChangePosition=m_configurationsNewMinusPositionIndex[configurationIndex];signChangePosition<m_numberOfParticles;
            ++signChangePosition)
        {
            // Build new configuration
            if(restorePreviousSign)
            {
                // Reset to + the - sign introduced in previous iteration (this operation is performed for all but first iteration)
                particleMomentaSumsElement += m_sortedParticleMomentaTimesTwo[signChangePosition-1];
            }
            else
            {
                // This is the first sign changing iteration on current starting configuration, no reset is required.
                // It will be necessary in the following ones, as signaled by the following flag
                restorePreviousSign = true;
            }
            // Compute particle momenta sum for the considered sign configuration
            particleMomentaSumsElement -= m_sortedParticleMomentaTimesTwo[signChangePosition];
            
            // Check momenta sum value
            if(particleMomentaSumsElement>=0.)
            {
                // Non negative momenta combination, update integrand function value and stored configuration data
                io_integrandFunctionValue += newConfigurationTotalSign*pow(particleMomentaSumsElement,m_numberOfParticlesLessThree);
                
                // Set position of sign change for the configuration derived from the current one
                m_configurationsNewMinusPositionIndex[storedConfigurationNumber] = signChangePosition+1;
                
                // Store total sign of current configuration
                m_configurationTotalSign[storedConfigurationNumber] = newConfigurationTotalSign;
                
                // Store momenta sum corresponding to current configuration
                m_particleMomentaSums[storedConfigurationNumber] = particleMomentaSumsElement;
                
                // Update number of accepted and stored configuration
                ++storedConfigurationNumber;
            }
            else
            {
                // Negative momenta sum, skip all the subsequent configuration obtained from the
                // reference one by adding a negative sign (this can be done thanks to momenta ordering)
                break;
            }
        }
    }
    
    // Update integrand function value with momenta independent factor
    io_integrandFunctionValue *= m_exactIntegrandFunctionFactor;
    
    return false;
}

void PhaseSpaceExactIntegrandFunction::updateCalculationData(void)
{
    m_numberOfParticlesLessThree = m_numberOfParticles-3;
    // TODO: avoid cast
    // Compute maximum number of sign configuration for current number of particles (2^m_numberOfParticles)
    m_maxNumberOfConfiguration = static_cast<unsigned int>(pow(2.,static_cast<int>(m_numberOfParticles)));
    
    // Update calculation storage structures with above maximum number of configuration
    m_configurationsNewMinusPositionIndex.resize(m_maxNumberOfConfiguration,0);
    m_configurationTotalSign.resize(m_maxNumberOfConfiguration,1.);
    m_sortedParticleMomenta.resize(m_numberOfParticles,0.);
    m_sortedParticleMomentaTimesTwo.resize(m_numberOfParticles,0.);
    m_particleMomentaSums.resize(m_maxNumberOfConfiguration);
    
    // Compute factor (-1/2^(N+1)*pi*(N-3)!)
    m_exactIntegrandFunctionFactor = factorial(m_numberOfParticlesLessThree);
    // TODO: avoid cast or factorize (check external factors or avoid recalculation)
    // TODO: Optimization: common factor exist for fixed particle number!
    m_exactIntegrandFunctionFactor *= pow(2.,static_cast<int>(m_numberOfParticles+1));
    m_exactIntegrandFunctionFactor *= pi;
    m_exactIntegrandFunctionFactor = -1./m_exactIntegrandFunctionFactor;
    
    return;
}
