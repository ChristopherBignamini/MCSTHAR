#ifndef HADRONDATA_H
#define HADRONDATA_H

/**
* @brief Hadron data class
*
* This class is used to store the hadron
* physical (+ ID) data
*
* @author Christopher Bignamini
*/
class HadronData
{
    public:

        /**
        * Constructor
        *
        */
        HadronData(void);
    
        /**
        * Constructor
        *
        * @param i_idCode Hadron identification code
        * @param i_strangeCharge Hadron strange charge
        * @param i_charmCharge Hadron charm charge
        * @param i_bottomCharge Hadron bottom charge
        * @param i_electricCharge Hadron electric charge
        * @param i_baryonicCharge Hadron baryonic charge
        * @param i_mass Hadron mass
        * @param i_spinMultiplicity Hadron spin multiplicity
        * @param i_strangeQuarkNumber Number of strange and antistrange constituent quarks
        * @param i_ssbarComponentValue Hadron strange (ssbar) wave function component value (square modulus)
        * @throw HadronizationException if non positive mass or ssbar wave function component not in range [0,1] is provided
        */
        HadronData(int i_idCode,
                   int i_strangeCharge,
                   int i_charmCharge,
                   int i_bottomCharge,
                   double i_electricCharge,
                   double i_baryonicCharge,
                   double i_mass,
                   unsigned int i_spinMultiplicity,
                   unsigned int i_strangeQuarkNumber,
                   double i_ssbarComponentValue = 0.);
    
        /**
        * Destructor
        */
        ~HadronData(void);

        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Hadron id code setter
        *
        * @param i_idCode Hadron id code
        */
        inline void setIdCode(const int i_idCode) { m_idCode = i_idCode; }
        
        /**
        * @brief Hadron strange charge setter
        *
        * @param i_strangeCharge Hadron strange charge
        */
        inline void setStrangeCharge(const int i_strangeCharge) { m_strangeCharge = i_strangeCharge; }
        
        /**
        * @brief Hadron charm charge setter
        *
        * @param i_charmCharge Hadron charm charge
        */
        inline void setCharmCharge(const int i_charmCharge) { m_charmCharge = i_charmCharge; }
                
        /**
        * @brief Hadron bottom charge setter
        *
        * @param i_bottomCharge Hadron bottom charge
        */
        inline void setBottomCharge(const int i_bottomCharge) { m_bottomCharge = i_bottomCharge; }
        
        /**
        * @brief Hadron electric charge setter
        *
        * @param i_electricCharge Hadron electric charge
        */
        inline void setElectricCharge(const double i_electricCharge) { m_electricCharge = i_electricCharge; }
        
        /**
        * @brief Hadron baryonic charge setter
        *
        * @param i_baryonicCharge Hadron baryonic charge
        */
        inline void setBaryonicCharge(const double i_baryonicCharge) { m_baryonicCharge = i_baryonicCharge; }

        /**
        * @brief Hadron mass setter
        *
        * @return i_mass hadron mass
        * @throw HadronizationException if non positive mass is provided
        */
        void setMass(double i_mass);
    
        /**
        * @brief Hadron spin multiplicity setter
        *
        * @return i_spin Hadron spin multiplicity
        */
        inline void setSpinMultiplicity(const unsigned int i_spinMultiplicity) { m_spinMultiplicity = i_spinMultiplicity; }
                
        /**
        * @brief Hadron spin multiplicity setter
        *
        * @param i_strangeQuarkNumber number of strange and antistrange contituent quarks
        */
        inline void setNumberOfStrangeQuarks(const unsigned int i_strangeQuarkNumber) { m_strangeQuarkNumber = i_strangeQuarkNumber; }

        /**
        * @brief Hadron wave function ssbar component square modulus setter
        *
        * @param i_ssbarComponentValue Wave function ssbar component square modulus
        * @throw HadronizationException if a wave function ssbar component square modulus not in [0,1] is provided
        */
        void setSSBARComponentValue(double i_ssbarComponentValue);
    
        /**
        * @brief Hadron id code getter
        *
        * @return Hadron id code
        */
        inline const int getIdCode(void) const { return m_idCode; }
        
        /**
        * @brief Hadron strange charge getter
        *
        * @return Hadron strange charge
        */
        inline int getStrangeCharge(void) const { return m_strangeCharge; }
        
        /**
        * @brief Hadron charm charge getter
        *
        * @return Hadron charm charge
        */
        inline int getCharmCharge(void) const { return m_charmCharge; }
        
        /**
        * @brief Hadron bottom charge getter
        *
        * @return Hadron bottom charge
        */
        inline int getBottomCharge(void) const { return m_bottomCharge; }
        
        /**
        * @brief Hadron electric charge getter
        *
        * @return Hadron electric charge
        */
        inline double getElectricCharge(void) const { return m_electricCharge; }
        
        /**
        * @brief Hadron baryonic charge getter
        *
        * @return Hadron baryonic charge
        */
        inline double getBaryonicCharge(void) const { return m_baryonicCharge; }
        
        /**
        * @brief Hadron mass getter
        *
        * @return Hadron mass
        */
        inline double getMass(void) const { return m_mass; }
        
        /**
        * @brief Hadron spin multiplicity getter
        *
        * @return Hadron spin multiplicity
        */
        inline unsigned int getSpinMultiplicity(void) const { return m_spinMultiplicity; }

        /**
        * @brief Number of strange and antistrange constituent quark getter
        *
        * @return Number of strange or antistrange contituent quarks
        */
        inline unsigned int getNumbeOfStrangeQuarks(void) const { return m_strangeQuarkNumber; }
        
        /**
        * @brief Hadron wave function ssbar component square modulus getter
        *
        * @return Hadron wave function ssbar component square modulus
        */
        inline double getSSBARComponentValue(void) const { return m_ssbarComponentValue; }
    
        /**
        * @brief 1 - hadron wave function ssbar component square modulus getter
        *
        * @return 1 - hadron wave function ssbar component square modulus
        */
        inline double getOneMinusSSBARComponentValue(void) const { return m_oneMinusSSBARComponentValue; }
    
    private:
    
        /**
        * Hadron identification code
        */
        int m_idCode;
    
        /**
        * Hadron strange charge
        */
        int m_strangeCharge;
        
        /**
        * Hadron charm charge
        */
        int m_charmCharge;
        
        /**
        * Hadron bottom charge
        */
        int m_bottomCharge;
        
        /**
        * Hadron electric charge
        */
        double m_electricCharge;
        
        /**
        * Hadron baryonic charge
        */
        double m_baryonicCharge;
        
        /**
        * Hadron mass
        */
        double m_mass;
        
        /**
        * Hadron spin multiplicity
        */
        unsigned int m_spinMultiplicity;
        
        /**
        * NUmber of strange and antistrange quarks
        */
        unsigned int m_strangeQuarkNumber;

        /**
        * Strange (ssbar) wave function component value (square modulus)
        */
        double m_ssbarComponentValue;
        
        /**
        * 1 - Strange (ssbar) wave function component value (square modulus) 
        */
        double m_oneMinusSSBARComponentValue;
    
};

#endif