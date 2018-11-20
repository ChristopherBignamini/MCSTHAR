#ifndef LORENTZVECTOR_H
#define LORENTZVECTOR_H

#include "Vector3D.h"

// TODO: check physics content of the present class (boost, etc...)
/**
* @brief Lorentz vector (metric +,-,-,-)
*
* @author Christopher Bignamini
*/
class LorentzVector
{
    public:
    
        /**
        * @brief Constructor
        */
        LorentzVector(void);

        /**
        * @brief Constructor
        *
        * @param i_t Time component value
        * @param i_s 3D space vector
        */
        LorentzVector(double i_t,
                      const Vector3D& i_s);

        /**
        * @brief Constructor
        *
        * @param i_t Time component value
        * @param i_x X component value
        * @param i_y Y component value
        * @param i_z Z component value
        */
        LorentzVector(double i_t,
                      double i_x,
                      double i_y,
                      double i_z);
    
        /**
        * @brief Destructor
        */
        ~LorentzVector(void);
    
        // TODO: default copy constructor and assignement operator are being used
    
        /**
        * @brief Time component value get method
        *
        * @return Time component value
        */
        inline double getT(void) const { return m_t; }

        /**
        * @brief Space 3D vector get method
        *
        * @return Space 3D vector
        */
        inline const Vector3D& getS(void) const { return m_s; }
        
        /**
        * @brief X component value get method
        *
        * @return X component value
        */
        inline double getX(void) const { return m_s.getX(); }

        /**
        * @brief Y component value get method
        *
        * @return Y component value
        */
        inline double getY(void) const { return m_s.getY(); }

        /**
        * @brief Z component value get method
        *
        * @return Z component value
        */
        inline double getZ(void) const { return m_s.getZ(); }

        /**
        * @brief Time component value set method
        *
        * @param i_t Time component value
        */
        inline void setT(const double i_t) { m_t = i_t; }
                
        /**
        * @brief Space components set method
        *
        * @param i_s Space 3D vector
        */
        inline void setS(const Vector3D& i_s) { m_s = i_s; }
    
        /**
        * @brief X component value set method
        *
        * @param i_x X component value
        */
        inline void setX(const double i_x) { m_s.setX(i_x); }
        
        /**
        * @brief Y component value set method
        *
        * @param i_y Y component value
        */
        inline void setY(const double i_y) { m_s.setY(i_y); }
        
        /**
        * @brief Z component value set method
        *
        * @param i_z Z component value
        */
        inline void setZ(const double i_z) { m_s.setZ(i_z); }
            
        /**
        * @brief Lorentz vector squared module get method
        *
        * @return Lorentz vector squared module
        */
        inline double getSquaredModule(void) const { return (m_t*m_t - m_s.getSquaredModule()); }// TODO: performance issues
    
        /**
        * @brief Lorentz vector boost transformation method
        *
        * @param i_boostVector Boost transformation vector
        * @throw HadronizationException in case of squared beta factor larger or equal to one
        */
        void boost(const Vector3D& i_boostVector);

        /**
        * @brief Lorentz vector boost transformation vector get method
        *
        * This method returns the 3D vector corresponding to a boost 
        * transformation into a reference frame moving as given by the
        * current Lorentz vector 
        *
        * @return Boost transformation vector
        * @throw HadronizationException in case of null time component Lorentz vector
        */
        Vector3D getBoostVector(void) const;
    
        /**
        * @brief Overloaded (component by component) += operator
        *
        * @param i_lorentzVector Vector being summed to the present one
        * @return current vector + input vector
        */
        LorentzVector& operator+=(const LorentzVector& i_lorentzVector);

        /**
        * @brief Overloaded (component by component) -= operator
        *
        * @param i_lorentzVector Vector being subtracted to the present one
        * @return current vector - input vector
        */
        LorentzVector& operator-=(const LorentzVector& i_lorentzVector);
    
        /**
        * @brief Overloaded (component by component) *= operator
        *
        * @param i_factor Vector multiplicative factor
        * @return current vector times input factor
        */
        LorentzVector& operator*=(double i_factor);
    
    private:
    
        /**
        * Time component
        */
        double m_t;
    
        /**
        * Space components
        */
        Vector3D m_s;
};

/**
* @brief Overloaded (component by component) sum operator
*
* @param i_firstLorentzVector Vector sum first addend
* @param i_secondLorentzVector Vector sum second addend
* @return Sum vector
*/
LorentzVector operator+(const LorentzVector& i_firstLorentzVector,
                        const LorentzVector& i_secondLorentzVector);

/**
* @brief Overloaded (component by component) subtraction operator
*
* @param i_firstLorentzVector Vector subtraction first addend
* @param i_secondLorentzVector Vector subtraction second addend
* @return Difference vector
*/
LorentzVector operator-(const LorentzVector& i_firstLorentzVector,
                        const LorentzVector& i_secondLorentzVector);

#endif