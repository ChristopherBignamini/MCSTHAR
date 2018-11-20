#include "../Include/PhaseSpaceIntegrator.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/Constants.h"
#include "../../../Utilities/Include/mathFunctions.h"
#include "../../../Utilities/Include/logMessage.h"
#include <sstream>
#include <limits>
#ifdef _USEROOTRANDOM
    #include "../../../Utilities/Include/ROOTRandomNumberGenerator.h"
#else
    #include "../../../Utilities/Include/GSLRandomNumberGenerator.h"// TODO: use preprocessor directives
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif

// TODO: debug
#include <cassert>
#include <iomanip>
#include <iostream>

using namespace std;

const unsigned int PhaseSpaceIntegrator::m_exactIntegrandCalculationMaxMultiplicity(10);

PhaseSpaceIntegrator::PhaseSpaceIntegrator(const unsigned int i_maxPhaseSpaceSamplingNumber,
                                           const double i_phaseSpaceIntegrationErrorThreshold,
                                           RandomNumberGenerator& io_randomGenerator)
                                          :m_isIntegratorReady(false)
                                          ,m_isIntegralAvailable(false)
                                          ,m_numberOfParticles(0)
                                          ,m_maxSamplingNumber(i_maxPhaseSpaceSamplingNumber)
                                          ,m_integrationSquaredErrorThreshold(i_phaseSpaceIntegrationErrorThreshold*
                                                                              i_phaseSpaceIntegrationErrorThreshold)
                                          ,m_minimumStabilityCheckSamplingNumber(100)
{
    if(i_maxPhaseSpaceSamplingNumber==0)
    {
        throw HadronizationException("Error in phase space integration setup, maximum number of sampling attempt set to zero",
                                     __FUNCTION__,
                                     826);
    }

    // Check provided number of minimum and maximum samplings
    if(m_maxSamplingNumber<m_minimumStabilityCheckSamplingNumber)
    {
        // Throw exception for bad sampling number configuration provided
        throw HadronizationException(
                                     "Error in phase space integration setup, Minimum number of phase space integration samplings larger than the provided maximum one",
                                     __FUNCTION__,826);
    }

    if(i_phaseSpaceIntegrationErrorThreshold<0.)
    {
        throw HadronizationException("Error in phase space integration setup, negative error integration relative error threshold provided",
                                     __FUNCTION__,
                                     826);
    }
    
    
#ifdef _OPENMP
    
    // Retrieve number of threads
    #pragma omp parallel// TODO: check for best practices
    {
        m_numberOfThreads = omp_get_num_threads();
    }
    
    // Compute single thread minimum sampling number
    if(m_minimumStabilityCheckSamplingNumber<m_numberOfThreads)
    {
        // If number of minimum sampling number is smaller than number of thread reset it to number of threads itself
        m_minimumStabilityCheckSamplingNumber = m_numberOfThreads;
        m_localMinStabilityCheckSamplingNumber = 1;
        
        // Log warning for bad sampling number configuration
        stringstream numberOfSamplingsStr;
        numberOfSamplingsStr<<m_minimumStabilityCheckSamplingNumber;
        logMessage("Minimum number of phase space integration samplings reset to " +
                   numberOfSamplingsStr.str() +
                   "(one per thread).",
                   WARNING);
    }
    else
    {
        // Check if minimum sampling number can be evenly distributed among threads, otherwise recompute it 
        const unsigned int samplingNumberPerThread(m_minimumStabilityCheckSamplingNumber/m_numberOfThreads);
        if((samplingNumberPerThread*m_numberOfThreads)!=m_minimumStabilityCheckSamplingNumber)
        {
            m_localMinStabilityCheckSamplingNumber = samplingNumberPerThread + 1;
            m_minimumStabilityCheckSamplingNumber = m_localMinStabilityCheckSamplingNumber*m_numberOfThreads;
            
            // Log warning for unbalanced sampling number configuration reset
            stringstream numberOfSamplingsStr;
            stringstream numberOfSamplingsPerThreadStr;
            numberOfSamplingsStr<<m_minimumStabilityCheckSamplingNumber;
            numberOfSamplingsPerThreadStr<<m_localMinStabilityCheckSamplingNumber;
            logMessage("Minimum number of phase space integration samplings reset to " +
                       numberOfSamplingsStr.str() +
                       " (" +
                       numberOfSamplingsPerThreadStr.str() +
                       " per thread).",
                       WARNING);

        }
        else
        {
            m_localMinStabilityCheckSamplingNumber = samplingNumberPerThread;        
        }
    }
    
    // Compute single thread maximum sampling number
    if(m_maxSamplingNumber<m_numberOfThreads)
    {
        // If number of maximum sampling number is smaller than number of thread reset it to number of threads itself
        m_maxSamplingNumber = m_numberOfThreads;
        m_localMaxSamplingNumber = 1;
        
        // Log warning for bad sampling number configuration
        stringstream numberOfSamplingsStr;
        numberOfSamplingsStr<<m_maxSamplingNumber;
        logMessage("Maximum number of phase space integration samplings reset to " +
                   numberOfSamplingsStr.str() +
                   "(one per thread). Bad thread configuration provided.",
                   WARNING);// TODO: throw exception in this case
    }
    else
    {
        // Check if maximum sampling number can be evenly distributed among threads, otherwise recompute it
        const unsigned int samplingNumberPerThread(m_maxSamplingNumber/m_numberOfThreads);
        if((samplingNumberPerThread*m_numberOfThreads)!=m_maxSamplingNumber)
        {
            m_localMaxSamplingNumber = samplingNumberPerThread+1;
            m_maxSamplingNumber = m_localMaxSamplingNumber*m_numberOfThreads;
            
            // Log warning for unbalanced sampling number configuration reset
            stringstream numberOfSamplingsStr;
            stringstream numberOfSamplingsPerThreadStr;
            numberOfSamplingsStr<<m_maxSamplingNumber;
            numberOfSamplingsPerThreadStr<<m_localMaxSamplingNumber;
            logMessage("Maximum number of phase space integration samplings reset to " +
                       numberOfSamplingsStr.str() +
                       " (" +
                       numberOfSamplingsPerThreadStr.str() +
                       " per thread).",
                       WARNING);
        }
        else
        {
            m_localMaxSamplingNumber = samplingNumberPerThread;
        }
    }
    
    // Build single thread random number generator and integrand functions
    m_phaseSpaceIntegrandFunctions.resize(m_numberOfThreads,NULL);
    m_phaseSpaceExactIntegrandFunctions.resize(m_numberOfThreads,NULL);
    m_phaseSpaceApproximatedIntegrandFunctions.resize(m_numberOfThreads,NULL);
    m_randomGenerators.resize(m_numberOfThreads,NULL);
    m_randomGenerators[0] = &io_randomGenerator;// TODO: only for debug??
    for(unsigned int threadNumber=0;threadNumber<m_numberOfThreads;++threadNumber)
    {
        m_phaseSpaceExactIntegrandFunctions[threadNumber] = new PhaseSpaceExactIntegrandFunction;
        m_phaseSpaceApproximatedIntegrandFunctions[threadNumber] = new PhaseSpaceApproximatedIntegrandFunction;
        if(threadNumber>0)
        {
            #ifdef _USEROOTRANDOM // TODO: move to final code (remove root)
                m_randomGenerators[threadNumber] = new ROOTRandomNumberGenerator;// TODO: use preprocessor directives
            #else   
                m_randomGenerators[threadNumber] = new GSLRandomNumberGenerator;// TODO: use preprocessor directives
            #endif
            
            m_randomGenerators[threadNumber]->setSeed(threadNumber);// TODO: check for best practices
        }
    }
    // TODO: check for best practices, must the seed be reassigne at the end of calculation?
    
#else
    
    m_phaseSpaceIntegrandFunction = NULL;
    m_phaseSpaceExactIntegrandFunction = new PhaseSpaceExactIntegrandFunction;
    m_phaseSpaceApproximatedIntegrandFunction = new PhaseSpaceApproximatedIntegrandFunction;
    m_randomGenerator = &io_randomGenerator;
    
#endif
    
    // TODO: add check and initialization status
    m_isIntegratorReady = true;
}

PhaseSpaceIntegrator::~PhaseSpaceIntegrator(void)
{
    
#ifdef _OPENMP

    for(unsigned int threadNumber=0;threadNumber<m_numberOfThreads;++threadNumber)
    {
        if(m_phaseSpaceExactIntegrandFunctions[threadNumber]!=NULL)
        {
            delete m_phaseSpaceExactIntegrandFunctions[threadNumber];
            m_phaseSpaceExactIntegrandFunctions[threadNumber] = NULL;
        }

        if(m_phaseSpaceApproximatedIntegrandFunctions[threadNumber]!=NULL)
        {
            delete m_phaseSpaceApproximatedIntegrandFunctions[threadNumber];
            m_phaseSpaceApproximatedIntegrandFunctions[threadNumber] = NULL;
        }

        if(threadNumber>0)// TODO: debug
        {
            if(m_randomGenerators[threadNumber]!=NULL)
            {
                delete m_randomGenerators[threadNumber];
                m_randomGenerators[threadNumber] = NULL;
            }
        }
    }
    
#else
    
    delete m_phaseSpaceExactIntegrandFunction;
    delete m_phaseSpaceApproximatedIntegrandFunction;
    
#endif
    
}

void PhaseSpaceIntegrator::run(const double i_availableEnergy,const vector<double>& i_particleMasses)
{
    m_isIntegralAvailable = false;
    if(m_isIntegratorReady)
    {
        try
        {
            // Check number of input particle masses
            m_numberOfParticles = i_particleMasses.size();
            
            if(m_numberOfParticles<2)
            {
                // Phase space integral not defined, throw exception
                throw HadronizationException("Error during phase space integration, less than 2 particle masses provided",
                                             __FUNCTION__,822);
            }
            else if(m_numberOfParticles==2)
            {
                // Two-body phase space, analytical integral available
                runTwoBodyCalculation(i_availableEnergy,i_particleMasses);
            }
            else
            {

#ifdef _OPENMP
                
                // Monte carlo integration required, choose integrand function calculation method
                if(m_numberOfParticles<=m_exactIntegrandCalculationMaxMultiplicity)
                {
                    for(unsigned int threadNumber=0;threadNumber<m_numberOfThreads;++threadNumber)
                    {
                        m_phaseSpaceIntegrandFunctions[threadNumber] = m_phaseSpaceExactIntegrandFunctions[threadNumber];
                    }
                }
                else
                {
                    for(unsigned int threadNumber=0;threadNumber<m_numberOfThreads;++threadNumber)
                    {
                        m_phaseSpaceIntegrandFunctions[threadNumber] = m_phaseSpaceApproximatedIntegrandFunctions[threadNumber];
                    }
                }
                
#else
                
                // Monte carlo integration required, choose integrand function calculation method
                if(m_numberOfParticles<=m_exactIntegrandCalculationMaxMultiplicity)
                {
                    m_phaseSpaceIntegrandFunction = m_phaseSpaceExactIntegrandFunction;
                }
                else
                {
                    m_phaseSpaceIntegrandFunction = m_phaseSpaceApproximatedIntegrandFunction;
                }
                
#endif
                // Run Monte Carlo integration
                runNBodyCalculation(i_availableEnergy,i_particleMasses);
            }
            
            // TODO: Include missing pi factor
            // TODO: can we avoid the calculation of this factor (check presence of matching terms in particle sampling)
        }
        catch(HadronizationException& ex)
        {
            throw ex;
        }
        m_isIntegralAvailable = true;
    }
    else
    {
        // Phase space integrator not correctly initialized, throw exception
        throw HadronizationException("Error during phase space integration, integrator not correctly initialized",
                                     __FUNCTION__,824);
    }
    return;
}

void PhaseSpaceIntegrator::runTwoBodyCalculation(const double i_availableEnergy,const vector<double>& i_particleMasses)
{
    // TODO: factorize with other methods
    const double phaseSpaceKineticEnergy(i_availableEnergy-i_particleMasses[0]-i_particleMasses[1]);
    
    // TODO: factorize with phase space generation class
    // Check phase space availability
    if(phaseSpaceKineticEnergy<=0.)
    {
        throw HadronizationException("Error during 2-body phase space integral calculation, no phase space availability for the provided masses",
                                     __FUNCTION__,823);
    }
    
    const double firstParticleSquaredMass(i_particleMasses[0]*i_particleMasses[0]);
    const double secondParticleSquaredMass(i_particleMasses[1]*i_particleMasses[1]);
    const double particleSquaredMassDifference(firstParticleSquaredMass-secondParticleSquaredMass);
    const double squaredAvailableEnergy(i_availableEnergy*i_availableEnergy);
    const double particleSquaredMomenta(0.25*(squaredAvailableEnergy -
                                              2.*(firstParticleSquaredMass+secondParticleSquaredMass) +
                                              particleSquaredMassDifference*particleSquaredMassDifference/squaredAvailableEnergy));
            
    m_phaseSpaceIntegrationData.phaseSpaceIntegral =
        4.*pi*sqrt(particleSquaredMomenta)*
        sqrt(particleSquaredMomenta+firstParticleSquaredMass)*
        sqrt(particleSquaredMomenta+secondParticleSquaredMass)/i_availableEnergy;    
    m_phaseSpaceIntegrationData.phaseSpaceIntegralError = 0.;
    m_phaseSpaceIntegrationData.isErrorUnderThreshold = true;
    m_phaseSpaceIntegrationData.numberOfSamplings = 0;
    
    return;
}

void PhaseSpaceIntegrator::runNBodyCalculation(const double i_availableEnergy,const vector<double>& i_particleMasses)
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

    // Check phase space availability
    if(phaseSpaceKineticEnergy<=0.)
    {        
        // No phase space available, throw exception
        string errorMessage("Error during N-body phase space integral calculation, no phase space availability for the provided masses. Mass values (GeV):\n");
        for(unsigned int massIndex=0;massIndex<m_numberOfParticles;++massIndex)
        {
            stringstream massStr;
            massStr<<setprecision(17);// TODO: switch to appropriate double precision
            massStr<<i_particleMasses[massIndex];
            errorMessage += "\n" + massStr.str();
        }
        stringstream massStr;
        massStr<<i_availableEnergy;
        errorMessage += "\n\nAvailable energy (GeV): " + massStr.str();
        
        throw HadronizationException(errorMessage,__FUNCTION__,823);
    }
    
    // Reset integration data
    // TODO: use structure method?
    m_phaseSpaceIntegrationData.phaseSpaceIntegral = 0.;
    m_phaseSpaceIntegrationData.phaseSpaceIntegralError = 0.;
    m_phaseSpaceIntegrationData.isErrorUnderThreshold = false;
    m_phaseSpaceIntegrationData.numberOfSamplings = 0;
    
    // Compute particle number (less one) factorial
    // TODO: factorize with factorial function
    // TODO: avoid factorial calculation when possible (can be stored as data member and reused if numbr of partcles do not change)
    // TODO: can be avoided at all? Check matching term in particle sampling
    // TODO: move outside and unify with the 2 body case
//    unsigned int particleNumberFactorial(1);
//    unsigned int particleCounter(1);
//    for(unsigned int particleIndex=0;particleIndex<m_numberOfParticles;++particleIndex,++particleCounter)
//    {
//        particleNumberFactorial *= particleCounter;
//    }
    const unsigned int numberOfParticlesLessOne(m_numberOfParticles-1);
    double particleNumberFactorial(factorial(numberOfParticlesLessOne));
    
    // Run Monte-Carlo integration
    const unsigned int numberOfParticlesLessTwo(numberOfParticlesLessOne-1);
    vector<double> powExponential(numberOfParticlesLessOne);// TODO: store somewhere precalculated pow?
    unsigned int particleIndexPlusOne(numberOfParticlesLessOne);
    for(unsigned int particleIndex=numberOfParticlesLessTwo;particleIndex>0;--particleIndex,--particleIndexPlusOne)
    {
        powExponential[particleIndex] = 1./particleIndexPlusOne;
    }
    powExponential[0] = 1.;
        
#ifdef _OPENMP
    
    // TODO: debug! unsigned long int or unsigned long long int? Is there a better calculation method to avoid huge numbers?
    unsigned long long int samplingNumber(0);
    
    // Single threads exception/warning collection
    vector<HadronizationException> singleThreadException(m_numberOfThreads);
    bool* singleThreadErrorStatus = new bool[m_numberOfThreads];
    
    // Global error status
    bool globalErrorStatus(false);
    
    #pragma omp parallel
    {
        // Retrieve current thread id
        const unsigned int threadId(omp_get_thread_num());
        
        // Single thread error status
        bool& errorStatus(singleThreadErrorStatus[threadId]);
        errorStatus = false;
        
        // Retrieve pointer to PhaseSpaceIntegrandFunction functor class for current thread
        // TODO: in case of OPENMP the vector of base pointers to the functor is no more required
        PhaseSpaceIntegrandFunction* phaseSpaceIntegrandFunction(m_phaseSpaceIntegrandFunctions[threadId]);

        // Retrieve current thread random number generator
        RandomNumberGenerator* randomGenerator(m_randomGenerators[threadId]);
        
        vector<double> energyFraction(m_numberOfParticles);
        energyFraction[numberOfParticlesLessOne] = 1.;
        vector<double> particleMomenta(m_numberOfParticles);

        // TODO: debug! unsigned long int or unsigned long long int? Is there a better calculation method to avoid huge numbers?
        unsigned long long int localSamplingNumber(0);
        int nextParticleIndex;
        unsigned int stabilityCheckSamplingNumber(m_localMinStabilityCheckSamplingNumber);
        unsigned int stabilityCheckSamplingNumberUpdate(m_localMinStabilityCheckSamplingNumber);
        unsigned int stabilityCheckSamplingNumberIncrement(0);
        // TODO: what about storing random number for next channels???
        bool functionCalculationExitCode;
        double localIntegrandFunctionSum;
        double localSquaredIntegrandFunctionSum;
        double integrandFunction;
        while(true)
        {

            localIntegrandFunctionSum  = 0.;
            localSquaredIntegrandFunctionSum = 0.;
            for(;localSamplingNumber<stabilityCheckSamplingNumber;++localSamplingNumber)
            {                
                // Run energy fraction sampling
                nextParticleIndex = numberOfParticlesLessOne;
                for(int particleIndex=numberOfParticlesLessTwo;particleIndex>=0;--particleIndex,--nextParticleIndex)
                {
                    // TODO: avoid pow(?,1)!
                    energyFraction[particleIndex] = pow(randomGenerator->getUniformRandom(),powExponential[particleIndex])*
                                                    energyFraction[nextParticleIndex];
                }
                // Run particle momenta and energy calculation
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
                    functionCalculationExitCode = (*phaseSpaceIntegrandFunction)(particleMomenta,
                                                                                 integrandFunction);// TODO: what about
                                                                                                    // numeric_limits<double>::epsilon() level
                                                                                                    // negative values? (should be zero?) Are zero
                                                                                                    // values possible??
                    // Check calcuation status
                    if(functionCalculationExitCode == true)// TODO: check if better strategy for calculation method switch exist
                    {
                        // Error during calculation, try changing calculation strategy
                        // NOTE: no switch from exact to approximated is currently foreseen
                        // TODO: comment the "if" below leaving the above comment..(also for serial case)
                        if(phaseSpaceIntegrandFunction == m_phaseSpaceApproximatedIntegrandFunctions[threadId])
                        {
                            // Switch from approximated to exact calculation method
                            phaseSpaceIntegrandFunction = m_phaseSpaceExactIntegrandFunctions[threadId];
                            
                            // Try again the calculation
                            functionCalculationExitCode = (*phaseSpaceIntegrandFunction)(particleMomenta,
                                                                                         integrandFunction);
                            
                            // Check again calculation status
                            if(functionCalculationExitCode == false)
                            {
                                // Calculation correctly performed, restore original calculation method
                                phaseSpaceIntegrandFunction = m_phaseSpaceApproximatedIntegrandFunctions[threadId];
                            }
                            else
                            {
                                // Error also with alternative method, throw exception
                                throw HadronizationException("Error during phase space integration: error in integrand function calculation after method switch from approximated to exact",
                                                             __FUNCTION__,828);// TODO: switch to 4 digit exit code..
                            }
                        }
                    }
                    
                    assert(localIntegrandFunctionSum>-numeric_limits<double>::epsilon());
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
                    // Retrieve thrown exception message
                    string errorMessage(ex.getErrorMessage());
                    
                    // Add current sampling number and thread id (//TODO: factorize into method?)
                    stringstream samplingNumberStr;
                    stringstream threadIdStr;
                    samplingNumberStr<<localSamplingNumber;
                    threadIdStr<<threadId;
                    
                    
                    // Update error message with sampling number
                    errorMessage += "\nLocal sampling number: " + samplingNumberStr.str();

                    // Update error message with trhead id
                    errorMessage += "\nThread id: " + threadIdStr.str();
                    
                    ex.setErrorMessage(errorMessage);
                    
                    // Store current thread exception
                    singleThreadException[threadId] = ex;
                    
                    // Update error status for current thread
                    errorStatus = true;
                    
                    // Break loop
                    break;
                }
                
                // TODO: check on integrand function sign? (e.g., avoid -3.e-15)
                // TODO: Optimization (mult, div, pow...)
                // Update local phase space weight and squared weight sum
                localIntegrandFunctionSum += integrandFunction;
                localSquaredIntegrandFunctionSum += integrandFunction*integrandFunction;
            }
            
            // Update global phase space weight and squared weight sum
            #pragma omp critical
            {
                // Check current thread error status
                if(errorStatus==false)
                {
                    // No error present, update integral
                    m_phaseSpaceIntegrationData.phaseSpaceIntegral += localIntegrandFunctionSum;
                    m_phaseSpaceIntegrationData.phaseSpaceIntegralError += localSquaredIntegrandFunctionSum;
                }
                else
                {
                    // Error found by current thread, update global error statu
                    globalErrorStatus = true;
                }
            }
            
            //Wait for all thread global sums update
            #pragma omp barrier

            // Check global error status
            if(globalErrorStatus==true)
            {
                // At least one thread reported an error, break loop
                break;
            }
            
            // Check error status (using the first available thread)
            #pragma omp single
            {
                // Update total sampling number
                samplingNumber = localSamplingNumber*m_numberOfThreads;
                
                // TODO: localSamplingNumber*phaseSpaceSquaredWeightSum can be very big???
                if(m_phaseSpaceIntegrationData.phaseSpaceIntegral>0.)
                {
                    // TODO: check unsigned int overflow everywhere!
                    const double oneOnWeightSumSquared(1./(m_phaseSpaceIntegrationData.phaseSpaceIntegral*
                                                           m_phaseSpaceIntegrationData.phaseSpaceIntegral));
                    double squaredRelativeError((samplingNumber*m_phaseSpaceIntegrationData.phaseSpaceIntegralError*oneOnWeightSumSquared - 1.)/
                                                (samplingNumber-1));
                    // Check integration error status
                    // TODO: check on relative or absolute error?
                    // TODO: check denominator
                    // TODO: add numerical check for /0, etc...
                    // TODO: avoid eps and use constant as for channel sampling?
                    // Check for numerical squared error non positiveness
                    if(std::abs(squaredRelativeError)<=numeric_limits<double>::epsilon())
                    {
                        // Squared error is approximatively zero, set it to exactly zero
                        squaredRelativeError = 0.;
                    }
                    else if(squaredRelativeError<-numeric_limits<double>::epsilon())
                    {
                        // TODO: multi threading exception handling support needed
                        // TODO: add global sampling number
                        // Squared error is negative, create exception
                        singleThreadException[threadId] =
                            HadronizationException("Error during phase space integration: negative relative squared error",
                                                     __FUNCTION__,828);
                        
                        // Update global and single thread error status 
                        errorStatus = true;
                        globalErrorStatus = true;
                    }
                    
                    // Check error status
                    if(errorStatus==false)
                    {
                        // TODO: optimize
                        if(squaredRelativeError<m_integrationSquaredErrorThreshold)
                        {
                            // Integration stability reached
                            m_phaseSpaceIntegrationData.isErrorUnderThreshold = true;
                        }
                    }
                }
            }

            // Check global error status
            if(globalErrorStatus==true)
            {
                // Error detected in previous single section
                break;
            }

            // Check if the maximum number of sampling number or stability has been reached
            if((samplingNumber==m_maxSamplingNumber) ||
               m_phaseSpaceIntegrationData.isErrorUnderThreshold)
            {
                break;
            }
            
            // TODO: check if better conter solution exists
            // Compute next sampling number value for stability check
            if(localSamplingNumber==stabilityCheckSamplingNumberUpdate)
            {
                stabilityCheckSamplingNumberUpdate = (localSamplingNumber*10);
                stabilityCheckSamplingNumberIncrement = stabilityCheckSamplingNumberUpdate/2;
                stabilityCheckSamplingNumber = stabilityCheckSamplingNumberIncrement;
            }
            else
            {
                stabilityCheckSamplingNumber += stabilityCheckSamplingNumberIncrement;
            }
            if(stabilityCheckSamplingNumber>m_localMaxSamplingNumber)
            {
                stabilityCheckSamplingNumber=m_localMaxSamplingNumber;
            }
        }
    
    }
        
    // Check if errors has been reported during calculation
    if(globalErrorStatus==true)
    {
        // At least one thread reported an error, build and throw exception
        HadronizationException globalException;
        for(unsigned int threadId=0;threadId<m_numberOfThreads;++threadId)
        {
            // Check if current thread reported an error
            if(singleThreadErrorStatus[threadId]==true)
            {
                globalException.merge(singleThreadException[threadId]);
            }
        }
        // Free memory
        delete [] singleThreadErrorStatus;
        
        // Set exception return value
        globalException.setReturnValue(827);
        throw globalException;
    }

    // Free memory // TODO: factorize with code above
    delete [] singleThreadErrorStatus;

#else// TODO: check best pratice! (code is duplicated)

    vector<double> energyFraction(m_numberOfParticles);
    energyFraction[numberOfParticlesLessOne] = 1.;
    vector<double> particleMomenta(m_numberOfParticles);
    // TODO: debug! unsigned long int or unsigned long long int? Is there a better calculation method to avoid huge numbers?
    unsigned int stabilityCheckSamplingNumber(m_minimumStabilityCheckSamplingNumber);
    unsigned int stabilityCheckSamplingNumberUpdate(m_minimumStabilityCheckSamplingNumber);
    unsigned int stabilityCheckSamplingNumberIncrement(0);
    // TODO: what about storing random number for next channels???
    double integrandFunction;
    bool functionCalculationExitCode;
    unsigned long long int samplingNumber(0);
    while(true)
    {        
        for(;samplingNumber<stabilityCheckSamplingNumber;++samplingNumber)
        {            
            // Run energy fraction sampling
            particleIndexPlusOne = numberOfParticlesLessOne;
            for(int particleIndex=numberOfParticlesLessTwo;particleIndex>=0;--particleIndex,--particleIndexPlusOne)
            {
                // TODO: avoid pow(?,1)!
                // TODO: debug
                energyFraction[particleIndex] = pow(m_randomGenerator->getUniformRandom(),powExponential[particleIndex])*energyFraction[particleIndexPlusOne];
            }
            
            // Run particle momenta and energy calculation
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
                        if(functionCalculationExitCode == false)
                        {
                            // Calculation correctly performed, restore original calculation method
                            m_phaseSpaceIntegrandFunction = m_phaseSpaceApproximatedIntegrandFunction;
                        }
                        else
                        {
                            // Error also with alternative method, throw exception
                            throw HadronizationException("Error during phase space integration: error in integrand function calculation after method switch from approximated to exact",
                                                         __FUNCTION__,828);// TODO: switch to 4 digit exit code..
                        }
                    }
                }
                
                assert(m_phaseSpaceIntegrationData.phaseSpaceIntegral>-numeric_limits<double>::epsilon());
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
            
            // TODO: check on integrand function sign? (e.g., avoid -3.e-15)
            // TODO: Optimization (mult, div, pow...)
            // Update phase space weight and squared weight sum
            m_phaseSpaceIntegrationData.phaseSpaceIntegral += integrandFunction;
            m_phaseSpaceIntegrationData.phaseSpaceIntegralError += integrandFunction*integrandFunction;
        }
        
        // Check error status
        // TODO: samplingNumber*phaseSpaceSquaredWeightSum can be very big???
        if(m_phaseSpaceIntegrationData.phaseSpaceIntegral>0.)
        {
            // TODO: check unsigned int overflow everywhere!
            const double oneOnWeightSumSquared(1./(m_phaseSpaceIntegrationData.phaseSpaceIntegral*
                                                   m_phaseSpaceIntegrationData.phaseSpaceIntegral));
            double squaredRelativeError(
                                        (samplingNumber*m_phaseSpaceIntegrationData.phaseSpaceIntegralError*oneOnWeightSumSquared - 1.)/(samplingNumber-1));
            // Check integration error status
            // TODO: check on relative or absolute error?
            // TODO: check denominator
            // TODO: add numerical check for /0, etc...
            // TODO: avoid eps and use constant as for channel sampling?
            // Check for numerical squared error non positiveness
            if(std::abs(squaredRelativeError)<=numeric_limits<double>::epsilon())
            {
                // Squared error is approximatively zero, set it to exactly zero
                squaredRelativeError = 0.;
            }
            else if(squaredRelativeError<-numeric_limits<double>::epsilon())
            {
                // Squared error is negative, throw exception
                throw HadronizationException("Error during phase space integration: negative relative squared error",
                                             __FUNCTION__,828);
            }
            
            
            // TODO: optimize
            if(squaredRelativeError<m_integrationSquaredErrorThreshold)
            {
                // Integration stability reached
                m_phaseSpaceIntegrationData.isErrorUnderThreshold = true;
                break;
            }
        }
        
        // Check if the maximum number of sampling number has been reached
        if(samplingNumber==m_maxSamplingNumber)
        {
            break;
        }

        
        // TODO: check if better conter solution exists
        // Compute next sampling number value for stability check
        if(samplingNumber==stabilityCheckSamplingNumberUpdate)
        {
            stabilityCheckSamplingNumberUpdate = (samplingNumber*10);
            stabilityCheckSamplingNumberIncrement = stabilityCheckSamplingNumberUpdate/2;
            stabilityCheckSamplingNumber = stabilityCheckSamplingNumberIncrement;
        }
        else
        {
            stabilityCheckSamplingNumber += stabilityCheckSamplingNumberIncrement;
        }
        if(stabilityCheckSamplingNumber>m_maxSamplingNumber)
        {
            stabilityCheckSamplingNumber=m_maxSamplingNumber;
        }        
    }

#endif

    // Compute integral value and error
    if(m_phaseSpaceIntegrationData.isErrorUnderThreshold==false)
    {        
        // Non positive phase space integral, some error occured
        if(m_phaseSpaceIntegrationData.phaseSpaceIntegral<=0.)
        {
            stringstream numberOfSamplingStr;
            numberOfSamplingStr<<samplingNumber;
            // TODO: add random seed in error message!
            const string errorMessage("Non positive phase space integral value after " +
                                      numberOfSamplingStr.str() +
                                      " sampling attempts");
            throw HadronizationException(errorMessage,__FUNCTION__,825);
        }
    }
    
    // TODO: check the casting below
    // TODO: check denominator value
    // TODO: Error->Variance
    // TODO: Threshold on absolute or relative error?
    // TODO: /N -> /(N-1)
    const double kineticFactor(pow(fourPi,static_cast<int>(m_numberOfParticles))*
                               pow(phaseSpaceKineticEnergy,static_cast<int>(numberOfParticlesLessOne))/particleNumberFactorial);
    
    // Check for numerical squared error non positiveness
    double scaledSquaredError((samplingNumber*m_phaseSpaceIntegrationData.phaseSpaceIntegralError -
                               m_phaseSpaceIntegrationData.phaseSpaceIntegral*m_phaseSpaceIntegrationData.phaseSpaceIntegral)/
                              ((1.*samplingNumber)*samplingNumber*(samplingNumber-1)));
    if(scaledSquaredError<0.)
    {
        if(scaledSquaredError>=-numeric_limits<double>::epsilon())
        {
            // Squared error is approximatively zero, set it to exactly zero
            scaledSquaredError = 0.;
        }
        else if(scaledSquaredError<-numeric_limits<double>::epsilon())
        {
            // Squared error is negative, throw exception
            throw HadronizationException("Error during phase space integration: negative squared error",
                                         __FUNCTION__,829);
        }
    }
    m_phaseSpaceIntegrationData.phaseSpaceIntegralError = kineticFactor*sqrt(scaledSquaredError);
    m_phaseSpaceIntegrationData.phaseSpaceIntegral *= kineticFactor/samplingNumber;
    m_phaseSpaceIntegrationData.numberOfSamplings = samplingNumber;
    
    return;
}

const PhaseSpaceIntegrationData& PhaseSpaceIntegrator::getIntegrationData(void) const
{
    if(m_isIntegralAvailable)
    {
        return m_phaseSpaceIntegrationData;
    }
    else
    {
        throw HadronizationException("Error during phase space integration data retrieval, integration data not available",
                                     __FUNCTION__,821);
    }
}