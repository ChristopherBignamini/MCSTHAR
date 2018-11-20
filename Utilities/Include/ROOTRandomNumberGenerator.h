#ifndef ROOTRANDOMNUMBERGENERATOR_H
#define ROOTRANDOMNUMBERGENERATOR_H

#include "RandomNumberGenerator.h"
#include "TRandom1.h"

using namespace std;

/**
* @brief ROOT Random number generator // TODO: class to be completed or removed
*
* @author Christopher Bignamini
*/
class ROOTRandomNumberGenerator : public RandomNumberGenerator
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
        ROOTRandomNumberGenerator(void);
        
        /**
        * @brief Constructor
        *
        * @param i_seed Random number generator initialization seed
        * @throw HadronizationException in case of error during generator building
        */
        ROOTRandomNumberGenerator(unsigned int i_seed);
            
        /**
        * @brief Destructor
        */
        ~ROOTRandomNumberGenerator(void);
    
        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Random number generator seed set method
        *
        * @param i_seed Random number generator initialization seed
        */
        void setSeed(unsigned int i_seed);
        
        /**
        * @brief Uniform distribution random number generation method
        *
        * @return Random number uniformingly extracted in [0,1)
        */
        double getUniformRandom(void);
        
        /**
        * @brief Poisson distribution random number generation method
        *
        * @param i_poissonMean Poisson distribution mean
        * @return Random number extracted according to a Poisson distribution with the provided mean value
        */
        double getPoissonRandom(const double i_poissonMean);
    
        // TODO: add check for mean positiveness
        /**
        * @brief Random number generator status writing method
        *
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status writing
        */
        void writeGeneratorStatus(const string& i_generatorStatusFile);
        
        /**
        * @brief Random number generator status setting from file method
        *
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status loading/setting
        */
        void setGeneratorStatusFromFile(const string& i_generatorStatusFile);
                    
    private:
                
        /**
        * Random number generator engine
        */
        TRandom1 m_rootRandomGenerator;
};

#endif