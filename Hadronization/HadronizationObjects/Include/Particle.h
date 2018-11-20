#ifndef PARTICLE_H
#define PARTICLE_H

#include "HadronizationObject.h"
#include "../../HadronizationHadronSet/Include/HadronData.h"

/**
* @brief Hadron representation class
*
* This class is derived from the HadronizationObject class and
* is used to represent the hadrons produced during the hadronization
* process. See HadronizationObject class description for on-shellness 
* condition.
*
* @author Christopher Bignamini
*/
class Particle : public HadronizationObject
{
    	    
    public:
        
        /**
        * @brief Constructor
        *
        */
        Particle(void);
    
        /**
        * @brief Constructor (null 4-momentum and mother index)
        *
        * @param i_hadronData Particle data
        * @param i_parentIndex Particle parent index in global event record
        * @param i_particleIndex Particle index in global event record
        */
        Particle(const HadronData& i_hadronData,unsigned int i_parentIndex,unsigned int i_particleIndex = 0);
        
        /**
        * @brief Destructor
        *
        */
        ~Particle(void);
        
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Particle parent index setter
        *
        * @param i_parentIndex Particle parent index in event record
        */
        inline void setParentIndex(const unsigned int i_parentIndex) { m_parentIndex = i_parentIndex; }
    
        /**
        * @brief Particle data setter
        *
        * This method is used to set the particle "hadronic" data. 4-momentum, parent 
        * index and object index are set to zero.
        *
        * @param i_hadronData Particle data
        */
        void setData(const HadronData& i_hadronData);

        /**
        * @brief Particle data setter
        *
        * This method is used to set the particle "hadronic" data and momentum. Parent 
        * index and object index are set to zero.
        *
        * @param i_hadronData Particle data
        * @param i_momentum Particle momentum 
        */
        void setDataAndMomentum(const HadronData& i_hadronData,const TLorentzVector& i_momentum);

        /**
        * @brief Particle spin multiplicity setter
        *
        * @param i_spinMultiplicity Particle spin multiplicity
        */
        inline void setSpinMultiplicity(const unsigned int i_spinMultiplicity) { m_spinMultiplicity = i_spinMultiplicity; }
                
        /**
        * @brief Particle spin multiplicity getter
        *
        * @return Particle spin multiplicity
        */
        inline unsigned int getSpinMultiplicity(void) const { return m_spinMultiplicity; }

        /**
        * @brief Particle wave function strange (ssbar) wave function component value
        *
        * @return Particle wave function strange (ssbar) wave function component value
        */
        inline double getSSBarComponent(void) const { return m_ssbarComponentValue; }
    
        /**
        * @brief Particle parent index getter
        *
        * @return Particle parent index in event record
        */
        inline unsigned int getParentIndex(void) const { return m_parentIndex; }

    
    private:
    
        /** 
        * Particle spin multiplicity
        */
        unsigned int m_spinMultiplicity;
    
        /**
        * Strange (ssbar) wave function component value (square modulus)
        */
        double m_ssbarComponentValue;
    
        /**
        * Particle parent index in event record
        */
        unsigned int m_parentIndex;
};

#endif