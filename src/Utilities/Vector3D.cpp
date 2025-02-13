#include "../Include/Vector3D.h"

Vector3D::Vector3D(void)
                  :m_x(0.)
                  ,m_y(0.)
                  ,m_z(0.)
{
}

Vector3D::Vector3D(const double i_x,
                   const double i_y,
                   const double i_z)
                  :m_x(i_x)
                  ,m_y(i_y)
                  ,m_z(i_z)
{
}

Vector3D::~Vector3D(void)
{
}

Vector3D& Vector3D::operator+=(const Vector3D& i_vector3D)
{
    m_x += i_vector3D.m_x;
    m_y += i_vector3D.m_y;
    m_z += i_vector3D.m_z;
    
    return *this;
}

Vector3D& Vector3D::operator-=(const Vector3D& i_vector3D)
{
    m_x -= i_vector3D.m_x;
    m_y -= i_vector3D.m_y;
    m_z -= i_vector3D.m_z;
    
    return *this;
}

Vector3D& Vector3D::operator*=(double i_factor)
{
    m_x *= i_factor;
    m_y *= i_factor;
    m_z *= i_factor;
    
    return *this;
}

Vector3D operator+(const Vector3D& i_firstVector3D,
                   const Vector3D& i_secondVector3D)
{
    Vector3D o_vector3D(i_firstVector3D);
    o_vector3D += i_secondVector3D;
    
    return o_vector3D;
}

Vector3D operator-(const Vector3D& i_firstVector3D,
                   const Vector3D& i_secondVector3D)
{
    Vector3D o_vector3D(i_firstVector3D);
    o_vector3D -= i_secondVector3D;
    
    return o_vector3D;
}

double operator*(const Vector3D& i_firstVector3D,
                 const Vector3D& i_secondVector3D)
{
    return (i_firstVector3D.getX()*i_secondVector3D.getX() +
            i_firstVector3D.getY()*i_secondVector3D.getY() +
            i_firstVector3D.getZ()*i_secondVector3D.getZ());
}

Vector3D operator*(double i_factor,
                   const Vector3D& i_vector3D)
{
    return Vector3D(i_vector3D.getX()*i_factor,
                    i_vector3D.getY()*i_factor,
                    i_vector3D.getZ()*i_factor);
}
