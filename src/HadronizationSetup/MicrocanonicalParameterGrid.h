#ifndef MICROCANONICALPARAMETERGRID_H
#define MICROCANONICALPARAMETERGRID_H

namespace MCSTHAR
{
    namespace Setup
    {

        /**
        * @brief Microcanonical parameter grid storage structure
        *
        * This structure is used to store the microcanonical parameter
        * (energy density and strangeness suppression parameter) grid
        * structure 
        *
        * @author Christopher Bignamini
        */
        // TODO: to be factorized with xml strorage stuctures
        struct MicrocanonicalParameterGrid
        {
            /**
            * @brief Constructor
            *
            * @param i_minEnergyDensity Minimum energy density value
            * @param i_maxEnergyDensity Maximum energy density value
            * @param i_numberOfEnergyDensityValues Number of energy density values
            * @param i_minGammaS Minimum strangeness suppression parameter value
            * @param i_maxGammaS Maximum strangeness suppression parameter value
            * @param i_numberOfGammaSValues Number of strangeness suppression parameter values
            */
            // TODO: add grid structure check and parameter automatic update (e.g. totalNumberOfGridPoints)
            MicrocanonicalParameterGrid(const double i_minEnergyDensity,
                                        const double i_maxEnergyDensity,
                                        const unsigned int i_numberOfEnergyDensityValues,
                                        const double i_minGammaS,
                                        const double i_maxGammaS,
                                        const unsigned int i_numberOfGammaSValues)
                                       :minEnergyDensity(i_minEnergyDensity)
                                       ,maxEnergyDensity(i_maxEnergyDensity)
                                       ,numberOfEnergyDensityValues(i_numberOfEnergyDensityValues)
                                       ,minGammaS(i_minGammaS)
                                       ,maxGammaS(i_maxGammaS)
                                       ,numberOfGammaSValues(i_numberOfGammaSValues)
                                       ,totalNumberOfGridPoints(i_numberOfEnergyDensityValues*i_numberOfGammaSValues)
            {
            }
            
            /**
            * Minimum energy density value
            */
            double minEnergyDensity;
            
            /**
            * Maximum energy density value
            */
            double maxEnergyDensity;
            
            /**
            * Number of energy density values
            */
            unsigned int numberOfEnergyDensityValues;
            
            /**
            * Minimum strangeness suppression parameter value
            */
            double minGammaS;
            
            /**
            * Maximum strangeness suppression parameter value
            */
            double maxGammaS;
            
            /**
            * Number of strangeness suppression parameter values
            */
            unsigned int numberOfGammaSValues;
            
            /**
            * Microcanonical parameter total number of grid points
            */
            unsigned int totalNumberOfGridPoints;
            
        };

    }
}
#endif
