#ifndef HADRONIZATIONOBJECT_H
#define HADRONIZATIONOBJECT_H

#include "TLorentzVector.h"

using namespace std;

// TODO: shellness not imposed

/**
* @brief Hadronization object (Cluster and Particle) base class
*
* Base class for the hadronization related objects, namely clusters
* and particles. No on mass shellness condition is imposed, except for 
* signaled cases, being mass and momentum independently setted.
*
* @author Christopher Bignamini
*/ 
class HadronizationObject
{
        
    public:
            
        /**
        * @brief Constructor
        *
        */
        HadronizationObject(void);

        /**
        * @brief Constructor
        *
        * Using this constructor the created object is on mass shell
        *
        * @param i_idCode Hadronization object identification code
        * @param i_strangeCharge Hadronization object strange charge
        * @param i_charmCharge Hadronization object charm charge
        * @param i_bottomCharge Hadronization object bottom charge
        * @param i_electricCharge Hadronization object electric charge
        * @param i_baryonicCharge Hadronization object baryonic charge
        * @param i_p Hadronization object 4-momentum
        * @param i_index Hadronization object index in event record
        */
        HadronizationObject(int i_idCode,
                            int i_strangeCharge,
                            int i_charmCharge,
                            int i_bottomCharge,
                            double i_electricCharge,
                            double i_baryonicCharge,
                            const TLorentzVector& i_p,
                            unsigned int i_index = 0);

        /**
        * @brief Constructor
        *
        * @param i_idCode Hadronization object identification code
        * @param i_strangeCharge Hadronization object strange charge
        * @param i_charmCharge Hadronization object charm charge
        * @param i_bottomCharge Hadronization object bottom charge
        * @param i_electricCharge Hadronization object electric charge
        * @param i_baryonicCharge Hadronization object baryonic charge
        * @param i_mass Hadronization object mass
        * @param i_index Hadronization object index in event record
        */
        HadronizationObject(int i_idCode,
                            int i_strangeCharge,
                            int i_charmCharge,
                            int i_bottomCharge,
                            double i_electricCharge,
                            double i_baryonicCharge,
                            double i_mass,
                            unsigned int i_index = 0);
    
        /**
        * @brief Destructor
        *
        */
        virtual ~HadronizationObject(void);
        
        // Default copy constructor and overloaded assignement operator are being used
        
        /**
        * @brief Object identification code setter
        *
        * @param i_idCode Hadronization object identification code
        */
        inline void setIdCode(const int i_idCode) { m_idCode = i_idCode; }

        /**
        * @brief Hadronization object index in event record setter
        *
        * @param i_objectIndex Hadronization object index in event record
        */
        inline void setIndex(const unsigned int i_objectIndex) { m_objectIndex = i_objectIndex; }
    
        /**
        * @brief Object strange charge setter
        *
        * @param i_strangeCharge Hadronization object strange charge
        */
        inline void setStrangeCharge(const int i_strangeCharge) { m_strangeCharge = i_strangeCharge; }
    
        /**
        * @brief Object charm charge setter
        *
        * @param i_charmCharge Hadronization object charm charge
        */
        virtual inline void setCharmCharge(const int i_charmCharge) { m_charmCharge = i_charmCharge; }        

        /**
        * @brief Object bottom charge setter
        *
        * @param i_bottomCharge Hadronization object bottom charge
        */
        virtual inline void setBottomCharge(const int i_bottomCharge) { m_bottomCharge = i_bottomCharge; }
    
        /**
        * @brief Object electric charge setter
        *
        * @param i_electricCharge Hadronization object electric charge
        */
        inline void setElectricCharge(const double i_electricCharge) { m_electricCharge = i_electricCharge; }

        /**
        * @brief Object baryonic charge setter
        *
        * @param i_baryonicCharge Hadronization object baryonic charge
        */
        inline void setBaryonicCharge(const double i_baryonicCharge) { m_baryonicCharge = i_baryonicCharge; }
        
        /**
        * @brief Object mass setter
        *
        * @param i_mass Hadronization object mass
        */
        inline void setMass(const double i_mass) { m_mass = i_mass; }

        /**
        * @brief Object 4-momentum setter
        *
        * @param i_p Hadronization object 4-momentum
        */
        inline void setP(const TLorentzVector &i_p) { m_p = i_p; }
    
        /**
        * @brief Object identification code getter
        *
        * @return Object identification code
        */
        inline int getIdCode(void) const { return m_idCode; }

        /**
        * @brief Hadronization object index in event record getter
        *
        * @return Hadronization object index in event record
        */
        inline unsigned int getIndex(void) const { return m_objectIndex; }
        
        /**
        * @brief Object strange charge getter
        *
        * @return Object strange charge value
        */
        inline int getStrangeCharge(void) const { return m_strangeCharge; }
        
        /**
        * @brief Object charm charge getter
        *
        * @return Object charm charge value
        */
        inline int getCharmCharge(void) const { return m_charmCharge; }
        
        /**
        * @brief Object bottom charge getter
        *
        * @return Object bottom charge value
        */
        inline int getBottomCharge(void) const { return m_bottomCharge; }
        
        /**
        * @brief Object electric charge getter
        *
        * @return Object electric charge value
        */
        inline double getElectricCharge(void) const { return m_electricCharge; }
        
        /**
        * @brief Object baryonic charge getter
        *
        * @return Object baryonic charge value
        */
        inline double getBaryonicCharge(void) const { return m_baryonicCharge; }
        
        /**
        * @brief Object mass getter
        *
        * @return Object mass value
        */
        inline double getMass(void) const { return m_mass; }
        
        /**
        * @brief Object 4-momentum getter
        *
        * @return Object 4-momentum
        */
        inline const TLorentzVector& getP(void) const { return m_p; }

    
    protected:
    
        /**
        * Object identification code
        */
        int m_idCode;

        /**
        * Object index in event record
        */
        unsigned int m_objectIndex;
    
        /**
        * Object strange charge
        */
        int m_strangeCharge;
        
        /**
        * Object charm charge
        */
        int m_charmCharge;
    
        /**
        * Object bottom charge
        */
        int m_bottomCharge;
    
        /** 
        * Object electric charge
        */
        double m_electricCharge;
    
        /** 
        * Object baryonic charge
        */
        double m_baryonicCharge;
    
        /**
        * Object mass
        */
        double m_mass;
    
        /**
        * Object 4-momentum
        */
        TLorentzVector m_p;
    
};

#endif