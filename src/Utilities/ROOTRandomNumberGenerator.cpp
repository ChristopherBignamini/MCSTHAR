#include "../Include/ROOTRandomNumberGenerator.h"
#include "../Include/HadronizationException.h"

ROOTRandomNumberGenerator::ROOTRandomNumberGenerator(void)
{
}

ROOTRandomNumberGenerator::ROOTRandomNumberGenerator(const unsigned int i_seed)
{
    setSeed(i_seed);
}

ROOTRandomNumberGenerator::~ROOTRandomNumberGenerator(void)
{
}

void ROOTRandomNumberGenerator::setSeed(unsigned int i_seed)
{
    m_rootRandomGenerator.SetSeed2(i_seed,4);
}

double ROOTRandomNumberGenerator::getUniformRandom(void)
{
    return m_rootRandomGenerator.Uniform(1.);
}

double ROOTRandomNumberGenerator::getPoissonRandom(const double i_poissonMean)
{
    return m_rootRandomGenerator.Poisson(i_poissonMean);
}

void ROOTRandomNumberGenerator::writeGeneratorStatus(const string& i_generatorStatusFile)
{
    // TODO: implement
    return;
}

void ROOTRandomNumberGenerator::setGeneratorStatusFromFile(const string& i_generatorStatusFile)
{
    // TODO: implement
    return;
}
