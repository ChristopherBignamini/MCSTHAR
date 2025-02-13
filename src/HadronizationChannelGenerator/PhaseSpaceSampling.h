#ifndef PHASESPACESAMPLING_H
#define PHASESPACESAMPLING_H

#include <vector>
#include "../../HadronizationHadronSet/Include/HadronData.h"
#include "../../HadronizationObjects/Include/Cluster.h"
#include "../../../Utilities/Include/RandomNumberGenerator.h"
#include "TLorentzVector.h"

using namespace std;

/**
* @brief Phase space configuration sampling class
*
* Phase space sampling class used to generate the kinematical configuration
* of a given cluster hadronization channel. The implemented generation strategy
* is based on a two-body decay chain (see for example Sec. 43.4 of 
* http://pdg.lbl.gov/2012/reviews/rpp2012-rev-kinematics.pdf). 
* The phase space configuration accessible by means of the getPhaseSpace method 
* below refers to the same reference frame of the provided cluster momentum and 
* all hadron momenta are on mass shell.
*
* @author Christopher Bignamini
*/
class PhaseSpaceSampling
{

    public:
        
        /**
        * @brief Constructor
        * 
        * This costructor does not require an input random number generator, needed by the internal
        * algorithms. It must be specified using the setRandomNumberGenerator method.        
        *
        */
        PhaseSpaceSampling(void);
    
        /**
        * @brief Constructor
        *
        * @param io_randomGenerator Random number generator used for phase configuration sampling
        * @throw HadronizationException in case of random number generator unavailability
        */
        PhaseSpaceSampling(RandomNumberGenerator& io_randomGenerator);
        
        /**
        * @brief Destructor
        *
        */
        ~PhaseSpaceSampling(void);
    
        // Default copy constructor and assignement operator are being used
    
        /**
        * @brief Phase space configuration sampling method
        *
        * Phase space configuration generation method: after a (successful) execution
        * of this method a phase space configuration for the provided input is available 
        * as well as a sampling weight. 
        *
        * @param i_clusterMomentum Hadronizing cluster 4-momentum
        * @param i_hadronMasses Hadronization channel hadron masses
        * @throw HadronizationException in case of random number generator unavailability, 
        * no energy availability for the considered mass configuration, less than two hadron 
        * masses provided or error during sampling algorithm execution
        */
        void run(const TLorentzVector& i_clusterMomentum, const vector<double>& i_hadronMasses);
    
        /**
        * @brief Phase space configuration sampling method (cluster and hadron data input)
        *
        * Phase space configuration generation method: analogous to the previous one, 
        * the only difference being the input abstraction level, in this case the 
        * required kinematical information are retrieved by the method itself.
        * 
        * @param i_cluster Hadronizing cluster
        * @param i_hadronData Hadronization channel hadron data
        * @throw HadronizationException in case of random number generator unavailability, 
        * no energy availability for the considered hadronization channel, less than two hadron
        * masses provided or error during sampling algorithm execution
        */
        void run(const Cluster& i_cluster, const vector<HadronData>& i_hadronData);
        
        /**
        * @brief Random number generator set method
        *
        * @param io_randomGenerator Random number generator
        * @throw HadronizationException in case of random number generator unavailability
        */
        void setRandomNumberGenerator(RandomNumberGenerator& io_randomGenerator);
    
        /**
        * @brief Phase space configuration get method
        *
        * This method returns the hadron momenta corresponding to the generated phase space 
        * configuration. The momenta contained in the returned vector are ordered accordingly
        * to the provided hadron masses.
        *
        * @return Hadron momenta
        * @throw HadronizationException if the phase space configuration is not available
        */
        const vector<TLorentzVector>& getPhaseSpace(void) const;
    
        /**
        * @brief Phase space configuration weight get method
        *
        * This method returns the sampling weight corresponding to the generated phase space
        * configuration.
        *
        * @return Phase space sampling weight
        * @throw HadronizationException if the phase space configuration is not available
        */
        double getSamplingWeight(void) const;
    
    private:
    
        /**
        * @brief Two-body split method
        *
        * This method is used to generate the 2-body splits in the implemented decay chain
        * simulation. The returned momentum, corresponding to the input hadron mass, is already 
        * in the correct reference frame, namely in the laboratory frame. All the information 
        * needed for the decay simulation not provided in input are retrieved from class data members.
        * 
        * @param i_hadronMass New momentum hadron mass
        * @param i_residualSquaredMass Multiple two-body decay decay residual squared invariant mass
        * @param i_cosTheta New momentum cos(theta) value, with theta = elevation angle
        * @param i_phi New momentum azimuth angle
        * @throw HadronizationException in case of negative hadron momentum module generation 
        * @return Generated 4-momentum for the provided hadron mass
        */
        TLorentzVector performTwoBodySplit(double i_hadronMass,double i_residualSquaredMass,double i_cosTheta,double i_phi);
	       
        /**
        * @brief Residual squared mass generation method
        *
        * This method in used in the implemented decay chain simulation for the generation
        * of the two-body split residual invariant masses
        *
        * @param i_hadronIndex Hadron index corresponding to the generated momentum
        * @param i_hadronMasses Hadron masses for the generated phase space
        * @throw HadronizationException in case of empty residual mass sampling interval  
        */
        double generateSquaredResidualMass(unsigned int i_hadronIndex,const vector<double>& i_hadronMasses);
    
        /**
        * Phase space configuration availability flag
        */
        bool m_isPhaseSpaceAvailable;
    
        /**
        * Generated hadron momenta
        */
        vector<TLorentzVector> m_hadronMomenta;

        /**
        * Sampling weight
        */
        double m_weight;

        /**
        * Random number generator (externally provided)
        */
        RandomNumberGenerator* m_randomGenerator;// TODO: using reference could improve performances

        /**
        * Phase space sampling jacobian factor
        */
        double m_jacobianFactor;

        /**
        * Variable residual momentum for decay chain
        */
        TLorentzVector m_residualMomentum;
    
};

#endif