#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

#include <string>
using namespace std;

/**
* @brief Random number generator interface class
*
* @author Christopher Bignamini
*/
class RandomNumberGenerator
{
    public:
    
        /**
        * @brief Constructor
        *
        * The usage of this constructor leads to the initialization
        * of the internal random number generator with its default
        * seed 
        *
        * @throw HadronizationException in case of error during generator building
        */
        RandomNumberGenerator(void);

        /**
        * @brief Constructor
        *
        * @param i_seed Random number generator initialization seed
        */
        RandomNumberGenerator(unsigned int i_seed);

        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Destructor
        */
        virtual ~RandomNumberGenerator(void);

        /**
        * @brief Random number generator seed set method
        *
        * @param i_seed Random number generator initialization seed
        */
        virtual void setSeed(unsigned int i_seed) = 0;
    
        /**
        * @brief Uniform distribution random number generation method 
        *
        * @return Random number uniformingly extracted in [0,1)
        */
        virtual double getUniformRandom(void) = 0;

        /**
        * @brief Uniform distribution random number generation method
        *
        * @param i_rangeMin Random number range minimum value
        * @param i_rangeMax Random number range maximum value
        * @return Random number uniformingly extracted in [i_rangeMin,i_rangeMax)
        */
        double getUniformRandom(int i_rangeMin,int i_rangeMax);
    
        /**
        * @brief Poisson distribution random number generation method
        *
        * @param i_poissonMean Poisson distribution mean
        * @return Random number extracted according to a Poisson distribution with the provided mean value
        */
        virtual double getPoissonRandom(const double i_poissonMean) = 0;
    
        /**
        * @brief Random number generator status writing method
        *
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status writing
        */
        virtual void writeGeneratorStatus(const string& i_generatorStatusFile) = 0;


        /**
        * @brief Random number generator status setting from file method
        *
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status loading/setting
        */
        virtual void setGeneratorStatusFromFile(const string& i_generatorStatusFile) = 0;
        
//    protected:
//
//        /**
//        * @brief Internal generator build method
//        *
//        * @throw HadronizationException in case of error during generator building
//        */
//        virtual void buildGenerator(void) = 0;
//
//        /**
//        * @brief Internal generator copy method
//        *
//        * @param i_randomNumberGenerator Random number generator to be copied
//        */
//        virtual void copyGenerator(const RandomNumberGenerator& i_randomNumberGenerator) = 0;        
};

#endif