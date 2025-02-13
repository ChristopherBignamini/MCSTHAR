#include "../Include/RandomNumberGenerator.h"

RandomNumberGenerator::RandomNumberGenerator(void)
{
}

RandomNumberGenerator::RandomNumberGenerator(const unsigned int i_seed)// TODO: to be removed
{
}

RandomNumberGenerator::~RandomNumberGenerator(void)// TODO: to be removed
{
}

double RandomNumberGenerator::getUniformRandom(const int i_rangeMin, const int i_rangeMax)
{
    return i_rangeMin + (i_rangeMax - i_rangeMin)*getUniformRandom();
}
