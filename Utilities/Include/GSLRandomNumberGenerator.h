#ifndef GSLRANDOMNUMBERGENERATOR_H
#define GSLRANDOMNUMBERGENERATOR_H

#include "RandomNumberGenerator.h"
#include "gsl/gsl_rng.h"

using namespace std;

/**
* @brief GSL Random number generator storage structure
*
* @author Christopher Bignamini
*/
struct GSLRandomNumberGeneratorStatus
{
    /**
    * Generator status data
    */
    void* generatorStatus;
    
    /**
    * Generator status data size
    */
    size_t generatorStatusSize;
};

// TODO: Implement usage of GSL environment variables for generator handling
/**
* @brief GSL Random number generator
*
* @author Christopher Bignamini
*/
class GSLRandomNumberGenerator : public RandomNumberGenerator
{
    public:
    
        /**
        * @brief Constructor
        *
        * The usage of this constructor leads to the initialization
        * of the internal random number generator with its default
        * seed
        *
        * @throw HadronizationException in case of error during generator building (e.g., no memory available (GSL error))
        */
        GSLRandomNumberGenerator(void);
        
        /**
        * @brief Constructor
        *
        * @param i_seed Random number generator initialization seed
        * @throw HadronizationException in case of error during generator building (e.g., no memory available (GSL error))
        */
        GSLRandomNumberGenerator(unsigned int i_seed);
        
        /**
        * @brief Copy constructor
        *
        * @param i_randomNumberGenerator Random number generator to be copied
        */
        GSLRandomNumberGenerator(const GSLRandomNumberGenerator& i_randomNumberGenerator);
        
        /**
        * @brief Destructor
        */
        ~GSLRandomNumberGenerator(void);
        
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
        * NOTE from GSL manual: "Since the data is written in the native
        *                        binary format it may not be portable between
        *                        different architectures."
        *
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status writing (GSL error)
        */
        void writeGeneratorStatus(const string& i_generatorStatusFile);
        
        /**
        * @brief Random number generator status writing method
        *
        * NOTE from GSL manual: "Since the data is written in the native
        *                        binary format it may not be portable between
        *                        different architectures."
        *
        * @param i_generatorStatus Random number generator status
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status writing
        */
        void writeGeneratorStatus(const GSLRandomNumberGeneratorStatus& i_generatorStatus,
                                  const string& i_generatorStatusFile);
        
        /**
        * @brief Random number generator status setting from file method
        *
        * NOTE from GSL manual: "The data is assumed to have been written in the native binary
        *                        format on the same architecture."
        *
        * WARNING: the generator type is not stored in the status file, it is therefore user
        *          responsibility of the user to use the stored status to set a generator
        *          of the same type of the one the status data refers to. (At present only the
        *          gsl_rng_ranlux389 generator is available).
        *
        *
        * @param i_generatorStatusFile Status storage file
        * @throw HadronizationException in case of error during status loading/setting (GSL error)
        */
        void setGeneratorStatusFromFile(const string& i_generatorStatusFile);
        
        /**
        * @brief Random number generator status getter
        *
        * This method can be used at debug time to store status of the generator
        * in order to reproduce the conditions which caused a given error.
        *
        * @return Current status of the random number generator
        */
        GSLRandomNumberGeneratorStatus getGeneratorStatus(void) const;
        
        /**
        * @brief Overloaded assignement operator
        *
        * @param i_randomNumberGenerator Random number generator to be used for assignement
        * @return Assigned random number generator
        */
        GSLRandomNumberGenerator& operator=(const GSLRandomNumberGenerator& i_randomNumberGenerator);
    
    private:
    
        /**
        * @brief Internal generator build method
        *
        * @throw HadronizationException in case of error during generator building (e.g., no memory available (GSL error))
        */
        void buildGenerator(void);
        
        /**
        * @brief Internal generator copy method
        *
        * @param i_randomNumberGenerator Random number generator to be copied
        */
        void copyGenerator(const GSLRandomNumberGenerator& i_randomNumberGenerator);
        
        /**
        * Random number generator engine
        */
        gsl_rng* m_gslRandomGenerator;
};

#endif