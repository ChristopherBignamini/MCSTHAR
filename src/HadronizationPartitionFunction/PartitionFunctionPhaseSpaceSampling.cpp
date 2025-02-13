#include "../Include/PartitionFunctionPhaseSpaceSampling.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/mathFunctions.h"
#include "../../../Utilities/Include/Constants.h"
#ifdef _USEROOTRANDOM
    #include "../../../Utilities/Include/ROOTRandomNumberGenerator.h"
#else
    #include "../../../Utilities/Include/GSLRandomNumberGenerator.h"// TODO: use preprocessor directives
#endif

// TODO: debug
#include <sstream>
#include <iomanip>
#include <iostream>

const unsigned int PartitionFunctionPhaseSpaceSampling::m_exactIntegrandCalculationMaxMultiplicity(10);

PartitionFunctionPhaseSpaceSampling::PartitionFunctionPhaseSpaceSampling(RandomNumberGenerator& io_randomGenerator)
                                                                        :m_numberOfParticles(0)
                                                                        ,m_phaseSpaceIntegrandFunction(NULL)
                                                                        ,m_phaseSpaceExactIntegrandFunction(
                                                                            new PhaseSpaceExactIntegrandFunction)
                                                                        ,m_phaseSpaceApproximatedIntegrandFunction(
                                                                            new PhaseSpaceApproximatedIntegrandFunction)
                                                                        ,m_randomGenerator(&io_randomGenerator)
{
}

PartitionFunctionPhaseSpaceSampling::~PartitionFunctionPhaseSpaceSampling(void)
{
    delete m_phaseSpaceExactIntegrandFunction;
    delete m_phaseSpaceApproximatedIntegrandFunction;
}

double PartitionFunctionPhaseSpaceSampling::run(const double i_availableEnergy,const vector<double>& i_particleMasses)
{
    try
    {
        // Check number of input particle masses
        m_numberOfParticles = i_particleMasses.size();
        
        if(m_numberOfParticles<2)
        {
            // Phase space configuration not defined, throw exception
            throw HadronizationException("Error during phase space configuration sampling, less than 2 particle masses provided",
                                         __FUNCTION__,821);
        }
        
        double o_weight;
        if(m_numberOfParticles==2)
        {
            // Two-body phase space
            o_weight = runTwoBodySampling(i_availableEnergy,i_particleMasses);
        }
        else
        {
            // N-body phase space
            // Choose integrand function calculation method
            if(m_numberOfParticles<=m_exactIntegrandCalculationMaxMultiplicity)
            {
                m_phaseSpaceIntegrandFunction = m_phaseSpaceExactIntegrandFunction;
            }
            else
            {
                m_phaseSpaceIntegrandFunction = m_phaseSpaceApproximatedIntegrandFunction;
            }

            // Run sampling
            o_weight = runNBodySampling(i_availableEnergy,i_particleMasses);
        }
        
        return o_weight;

        // TODO: Include missing pi factor
        // TODO: can we avoid the calculation of this factor (check presence of matching terms in particle sampling)
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

double PartitionFunctionPhaseSpaceSampling::runTwoBodySampling(const double i_availableEnergy,
                                                               const vector<double>& i_particleMasses) const
{
    // TODO: factorize with other methods
    const double phaseSpaceKineticEnergy(i_availableEnergy-i_particleMasses[0]-i_particleMasses[1]);
    
    // TODO: factorize with phase space generation class
    // Check phase space availability
    if(phaseSpaceKineticEnergy<=0.)
    {
        throw HadronizationException("Error during 2-body phase space configuration sampling, no phase space availability for the provided masses",
                                     __FUNCTION__,822);
    }
    
    const double firstParticleSquaredMass(i_particleMasses[0]*i_particleMasses[0]);
    const double secondParticleSquaredMass(i_particleMasses[1]*i_particleMasses[1]);
    const double particleSquaredMassDifference(firstParticleSquaredMass-secondParticleSquaredMass);
    const double squaredAvailableEnergy(i_availableEnergy*i_availableEnergy);
    const double particleSquaredMomenta(0.25*(squaredAvailableEnergy -
                                              2.*(firstParticleSquaredMass+secondParticleSquaredMass) +
                                              particleSquaredMassDifference*particleSquaredMassDifference/squaredAvailableEnergy));
    
    // TODO: check calculation
    return (fourPi*
            sqrt(particleSquaredMomenta)*
            sqrt(particleSquaredMomenta+firstParticleSquaredMass)*
            sqrt(particleSquaredMomenta+secondParticleSquaredMass)/i_availableEnergy);
}

double PartitionFunctionPhaseSpaceSampling::runNBodySampling(const double i_availableEnergy,
                                                             const vector<double>& i_particleMasses)
{
    // Compute available kinetic energy and store 2*particle mass
    vector<double> particleMassesTimesTwo(m_numberOfParticles);
    vector<double> particleSquaredMasses(m_numberOfParticles);
    double phaseSpaceKineticEnergy(0.);
    for(unsigned int particleIndex=0;particleIndex<m_numberOfParticles;++particleIndex)
    {
        const double& currentParticleMass(i_particleMasses[particleIndex]);
        phaseSpaceKineticEnergy += currentParticleMass;
        particleMassesTimesTwo[particleIndex] = 2.*currentParticleMass;
        particleSquaredMasses[particleIndex] = currentParticleMass*currentParticleMass;
    }
    phaseSpaceKineticEnergy = i_availableEnergy - phaseSpaceKineticEnergy;
    
    // TODO: factorize with 2 body method
    // Check phase space availability
    if(phaseSpaceKineticEnergy<=0.)
    {        
        // No phase space available, throw exception
        string errorMessage("Error during N-body phase space configuration sampling, no phase space availability for the provided masses. Mass values (GeV):\n");
        for(unsigned int massIndex=0;massIndex<m_numberOfParticles;++massIndex)
        {
            stringstream massStr;
            massStr<<setprecision(17);// TODO: switch to appropriate double precision
            massStr<<i_particleMasses[massIndex];// TODO: add same message for 2 body case
            errorMessage += "\n" + massStr.str();
        }
        stringstream massStr;
        massStr<<i_availableEnergy;
        errorMessage += "\n\nAvailable energy (GeV): " + massStr.str();
        
        throw HadronizationException(errorMessage,__FUNCTION__,823);
    }
                
    const unsigned int numberOfParticlesLessOne(m_numberOfParticles-1);
    const unsigned int numberOfParticlesLessTwo(numberOfParticlesLessOne-1);
    vector<double> powExponential(numberOfParticlesLessOne);// TODO: store somewhere precalculated pow?
    unsigned int particleIndexPlusOne(numberOfParticlesLessOne);
    for(unsigned int particleIndex=numberOfParticlesLessTwo;particleIndex>0;--particleIndex,--particleIndexPlusOne)
    {
        powExponential[particleIndex] = 1./particleIndexPlusOne;
    }
    powExponential[0] = 1.;
        
    // Run energy fraction sampling
    vector<double> energyFraction(m_numberOfParticles);
    energyFraction[numberOfParticlesLessOne] = 1.;
    particleIndexPlusOne = numberOfParticlesLessOne;
    for(int particleIndex=numberOfParticlesLessTwo;particleIndex>=0;--particleIndex,--particleIndexPlusOne)
    {
        // TODO: avoid pow(?,1)!
        // TODO: debug
        energyFraction[particleIndex] = pow(m_randomGenerator->getUniformRandom(),powExponential[particleIndex])*energyFraction[particleIndexPlusOne];
    }
    
    // Run particle momenta and energy calculation
    vector<double> particleMomenta(m_numberOfParticles);
    double integrandFunction;
    bool functionCalculationExitCode;
    double particleKineticEnergy(energyFraction[0]*phaseSpaceKineticEnergy);// TODO: factorize these operation into a function..
    double particleSquaredMomenta(particleKineticEnergy*(particleKineticEnergy+particleMassesTimesTwo[0]));
    double particleEnergyFactor(sqrt(particleSquaredMasses[0]+particleSquaredMomenta));
    particleMomenta[0] = sqrt(particleSquaredMomenta);
    unsigned int particleIndexLessOne(0);
    for(unsigned int particleIndex=1;particleIndex<numberOfParticlesLessOne;++particleIndex,++particleIndexLessOne)
    {
        particleKineticEnergy = (energyFraction[particleIndex] - energyFraction[particleIndexLessOne])*phaseSpaceKineticEnergy;
        particleSquaredMomenta = particleKineticEnergy*(particleKineticEnergy+particleMassesTimesTwo[particleIndex]);
        particleEnergyFactor *= sqrt(particleSquaredMasses[particleIndex]+particleSquaredMomenta);
        particleMomenta[particleIndex] = sqrt(particleSquaredMomenta);
    }
    particleKineticEnergy = (1. - energyFraction[numberOfParticlesLessTwo])*phaseSpaceKineticEnergy;
    particleSquaredMomenta = particleKineticEnergy*(particleKineticEnergy+particleMassesTimesTwo[numberOfParticlesLessOne]);
    particleEnergyFactor *= sqrt(particleSquaredMasses[numberOfParticlesLessOne]+particleSquaredMomenta);
    particleMomenta[numberOfParticlesLessOne] = sqrt(particleSquaredMomenta);
    
    // Compute integrand function
    try
    {                
        // Try to compute the integrand function with current method
        functionCalculationExitCode = (*m_phaseSpaceIntegrandFunction)(particleMomenta,
                                                                       integrandFunction);
        // Check calcuation status
        if(functionCalculationExitCode == true)
        {
            // Error during calculation, try changing calculation strategy
            // NOTE: no switch from exact to approximated is currently foreseen
            if(m_phaseSpaceIntegrandFunction == m_phaseSpaceApproximatedIntegrandFunction)
            {
                // Switch from approximated to exact calculation method
                m_phaseSpaceIntegrandFunction = m_phaseSpaceExactIntegrandFunction;
                
                // Try again the calculation
                functionCalculationExitCode = (*m_phaseSpaceIntegrandFunction)(particleMomenta,
                                                                               integrandFunction);
                
                // Check again calculation status
                if(functionCalculationExitCode == true)
                {
                    // Error also with alternative method, throw exception
                    throw HadronizationException("Error during phase space configuration sampling: error in integrand function calculation after method switch from approximated to exact",
                                                 __FUNCTION__,824);// TODO: switch to 4 digit exit code..
                }
            }
        }
        
        // Integrand function can assume negative value due to numerical error, since
        // in these case the function value should be zero. 
        if(integrandFunction<0.)
        {
            integrandFunction = 0.;
        }
        else
        {
            integrandFunction *= particleEnergyFactor;
        }
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
    // TODO: avoid factorial calculation when possible (can be stored as data member and reused if numbr of partcles do not change)
    // TODO: can be avoided at all? Check matching term in particle sampling
    // TODO: (factorial) move outside and unify with the 2 body case
    // TODO: check on integrand function sign? (e.g., avoid -3.e-15)
    // TODO: Optimization (mult, div, pow...)            
    // TODO: check the casting below
    // TODO: check denominator value
    return (integrandFunction*
            pow(fourPi,static_cast<int>(m_numberOfParticles))*
            pow(phaseSpaceKineticEnergy,static_cast<int>(numberOfParticlesLessOne))/
            factorial(numberOfParticlesLessOne));
}
