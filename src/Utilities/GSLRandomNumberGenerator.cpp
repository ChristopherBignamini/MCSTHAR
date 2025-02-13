#include "../Include/GSLRandomNumberGenerator.h"
#include "../Include/HadronizationException.h"
#include "gsl/gsl_randist.h"

GSLRandomNumberGenerator::GSLRandomNumberGenerator(void)
{
    buildGenerator();
}

GSLRandomNumberGenerator::GSLRandomNumberGenerator(const unsigned int i_seed)
{
    buildGenerator();
    setSeed(i_seed);
}

GSLRandomNumberGenerator::GSLRandomNumberGenerator(const GSLRandomNumberGenerator& i_randomNumberGenerator)
{
    copyGenerator(i_randomNumberGenerator);
}

GSLRandomNumberGenerator::~GSLRandomNumberGenerator(void)
{
    if(m_gslRandomGenerator!=NULL)
    {
        gsl_rng_free(m_gslRandomGenerator);
    }
}

void GSLRandomNumberGenerator::setSeed(unsigned int i_seed)
{
    gsl_rng_set(m_gslRandomGenerator,i_seed);
}

double GSLRandomNumberGenerator::getUniformRandom(void)
{
    return gsl_rng_uniform(m_gslRandomGenerator);
}

double GSLRandomNumberGenerator::getPoissonRandom(const double i_poissonMean)
{
    return gsl_ran_poisson(m_gslRandomGenerator,i_poissonMean);
}

void GSLRandomNumberGenerator::buildGenerator(void)
{
    m_gslRandomGenerator = gsl_rng_alloc(gsl_rng_ranlux389);
    
    // Check random generator creation status
    if(m_gslRandomGenerator==NULL)
    {
        // Error during random generator creation, throw exception
        throw HadronizationException("Error during GSL random number generator creation, insufficient memory availability",
                                     __FUNCTION__,721);
    }
    
    return;
}

void GSLRandomNumberGenerator::writeGeneratorStatus(const string& i_generatorStatusFile)
{
    // Create/open C stream file in binary mode
    FILE* generatorStatusStream;
    generatorStatusStream = fopen(i_generatorStatusFile.c_str(), "wb");
    
    // Write status on stream
    const int writeExitCode(gsl_rng_fwrite(generatorStatusStream,m_gslRandomGenerator));
    
    // Close stream
    fclose(generatorStatusStream);
    
    // Check exit code
    if(writeExitCode!=0)
    {
        // Error during generator status writing, throw exception
        throw HadronizationException("Error during GSL random number generator status writing on file",
                                     __FUNCTION__,722);
    }
    
    return;
}

void GSLRandomNumberGenerator::writeGeneratorStatus(const GSLRandomNumberGeneratorStatus& i_generatorStatus,
                                                    const string& i_generatorStatusFile)
{
    // Create/open C stream file in binary mode
    FILE* generatorStatusStream;
    generatorStatusStream = fopen(i_generatorStatusFile.c_str(), "rb");
    
    // Write status on file
    const size_t writtenDataSize(fwrite(i_generatorStatus.generatorStatus,
                                        i_generatorStatus.generatorStatusSize,
                                        1,
                                        generatorStatusStream));
    
    // Close stream
    fclose(generatorStatusStream);
    
    // Check written data size for writing error
    if(writtenDataSize != i_generatorStatus.generatorStatusSize)
    {
        // Error during generator status writing, throw exception
        throw HadronizationException("Error during GSL random number generator status writing on file",
                                     __FUNCTION__,724);
    }
    
    return;
}

void GSLRandomNumberGenerator::setGeneratorStatusFromFile(const string& i_generatorStatusFile)
{
    // Create/open C stream file in binary mode
    FILE* generatorStatusStream;
    generatorStatusStream = fopen(i_generatorStatusFile.c_str(), "rb");
    
    // Set current generator with the stored status
    const int readAndSetExitCode(gsl_rng_fread(generatorStatusStream,m_gslRandomGenerator));
    
    // Close stream
    fclose(generatorStatusStream);
    
    // Check exit code
    if(readAndSetExitCode!=0)
    {
        // Error during generator status reading and setting, throw exception
        throw HadronizationException("Error during GSL random number generator status setting from file",
                                     __FUNCTION__,723);
        
    }
    
    return;
}

GSLRandomNumberGeneratorStatus GSLRandomNumberGenerator::getGeneratorStatus(void) const
{
    GSLRandomNumberGeneratorStatus o_generatorStatus;
    
    // Retrieve generator status
    o_generatorStatus.generatorStatus = gsl_rng_state(m_gslRandomGenerator);
    
    // Retrieve generator status size
    o_generatorStatus.generatorStatusSize = gsl_rng_size(m_gslRandomGenerator);
    
    return o_generatorStatus;
}

GSLRandomNumberGenerator& GSLRandomNumberGenerator::operator=(const GSLRandomNumberGenerator& i_randomNumberGenerator)
{
    copyGenerator(i_randomNumberGenerator);
    
    return *this;
}

void GSLRandomNumberGenerator::copyGenerator(const GSLRandomNumberGenerator& i_randomNumberGenerator)
{
    // Copy GSL random generator
    m_gslRandomGenerator = gsl_rng_clone(i_randomNumberGenerator.m_gslRandomGenerator);
    
    return;
}
