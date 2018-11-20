#ifndef VECTOR3D_H
#define VECTOR3D_H

/**
* @brief 3D euclidean vector
*
* @author Christopher Bignamini
*/
class Vector3D
{
    public:
    
        /**
        * @brief Constructor
        */
        Vector3D(void);
        
        /**
        * @brief Constructor
        *
        * @param i_x X component value
        * @param i_y Y component value
        * @param i_z Z component value
        */
        Vector3D(double i_x,
                 double i_y,
                 double i_z);
    
        /**
        * @brief Destructor
        */
        ~Vector3D(void);
        
        // TODO: default copy constructor and assignement operator are being used
    
        /**
        * @brief X component value get method
        *
        * @return X component value
        */
        inline double getX(void) const { return m_x; }
        
        /**
        * @brief Y component value get method
        *
        * @return Y component value
        */
        inline double getY(void) const { return m_y; }
        
        /**
        * @brief Z component value get method
        *
        * @return Z component value
        */
        inline double getZ(void) const { return m_z; }
                
        /**
        * @brief X component value set method
        *
        * @param i_x X component value
        */
        inline void setX(const double i_x) { m_x = i_x; }
        
        /**
        * @brief Y component value set method
        *
        * @param i_y Y component value
        */
        inline void setY(const double i_y) { m_y = i_y; }
        
        /**
        * @brief Z component value set method
        *
        * @param i_z Z component value
        */
        inline void setZ(const double i_z) { m_z = i_z; }
        
        /**
        * @brief Vector squared module get method
        *
        * @return Vector squared module
        */
        inline double getSquaredModule(void) const { return (m_x*m_x + m_y*m_y + m_z*m_z); }// TODO: performance issues
            
        /**
        * @brief Overloaded (component by component) += operator
        *
        * @param i_vector3D Vector being summed to the present one 
        * @return current vector + input vector
        */
        Vector3D& operator+=(const Vector3D& i_vector3D);
        
        /**
        * @brief Overloaded (component by component) -= operator
        *
        * @param i_vector3D Vector being subtracted to the present one
        * @return current vector - input vector
        */
        Vector3D& operator-=(const Vector3D& i_vector3D);
    
        /**
        * @brief Overloaded (component by component) *= operator
        *
        * @param i_factor Vector multiplicative factor
        * @return current vector times input factor
        */
        Vector3D& operator*=(double i_factor);
    
    private:
    
        /**
        * X component
        */
        double m_x;

        /**
        * Y component
        */
        double m_y;

        /**
        * Z component
        */
        double m_z;

};

/**
* @brief Overloaded (component by component) sum operator
*
* @param i_firstVector3D Vector sum first addend
* @param i_secondVector3D Vector sum second addend
* @return Sum vector
*/
Vector3D operator+(const Vector3D& i_firstVector3D,
                   const Vector3D& i_secondVector3D);

/**
* @brief Overloaded (component by component) subtraction operator
*
* @param i_firstVector3D Vector subtraction first addend
* @param i_secondVector3D Vector subtraction second addend
* @return Difference vector
*/
Vector3D operator-(const Vector3D& i_firstVector3D,
                   const Vector3D& i_secondVector3D);

/**
* @brief Overloaded scalar product operator
*
* @param i_firstVector3D Scalar product first factor
* @param i_secondVector3D Scalar product second factor
* @return Difference vector
*/
double operator*(const Vector3D& i_firstVector3D,
                 const Vector3D& i_secondVector3D);

/**
* @brief Overloaded (component by component) product by factor operator
*
* @param i_factor Vector multiplicative factor
* @param i_vector3D Vector to be multiplied by input factor
* @return Input vector times input factor
*/
Vector3D operator*(double i_factor,
                   const Vector3D& i_vector3D);

#endif