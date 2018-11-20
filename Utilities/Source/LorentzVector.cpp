#include "../Include/LorentzVector.h"
#include "../Include/HadronizationException.h"
#include <cmath>

LorentzVector::LorentzVector(void)
                            :m_t(0.)
{
}
    
LorentzVector::LorentzVector(const double i_t,
                             const Vector3D& i_s)
                            :m_t(i_t)
                            ,m_s(i_s)
{
}

LorentzVector::LorentzVector(const double i_t,
                             const double i_x,
                             const double i_y,
                             const double i_z)
                            :m_t(i_t)
                            ,m_s(i_x,
                                 i_y,
                                 i_z)
{
}

LorentzVector::~LorentzVector(void)
{
}

void LorentzVector::boost(const Vector3D& i_boostVector)
{
    // TODO: add equations
    const double squaredBeta(i_boostVector.getSquaredModule());
    const double oneLessSquaredBeta(1. - squaredBeta);
    if(oneLessSquaredBeta>0.)
    {
        const double gamma(1./sqrt(oneLessSquaredBeta));
        const double boostTimesSpace(i_boostVector*m_s);
        const double gammaTimesTime(gamma*m_t);
        
        // Compute space components transformation
        m_s = m_s + ((gamma - 1.)*boostTimesSpace/squaredBeta + gammaTimesTime)*i_boostVector;

        // Compute time components transformation
        m_t = gammaTimesTime + gamma*boostTimesSpace;
    }
    else
    {
        throw HadronizationException("Error during boost calculation, squared beta factor larger or equal to one",
                                     __FUNCTION__,742);
    }
}

Vector3D LorentzVector::getBoostVector(void) const
{
    if(m_t!=0.)
    {
        return Vector3D((1./m_t)*m_s);
    }
    else
    {
        throw HadronizationException("Error during boost vector calculation from Lorentz vector, zero time component",
                                     __FUNCTION__,741);
    }
}

LorentzVector& LorentzVector::operator+=(const LorentzVector& i_lorentzVector)
{
    m_t += i_lorentzVector.m_t;
    m_s += i_lorentzVector.m_s;
    
    return *this;
}

LorentzVector& LorentzVector::operator-=(const LorentzVector& i_lorentzVector)
{
    m_t -= i_lorentzVector.m_t;
    m_s -= i_lorentzVector.m_s;
    
    return *this;
}

LorentzVector& LorentzVector::operator*=(const double i_factor)
{
    m_t *= i_factor;
    m_s *= i_factor;
    
    return *this;
}

LorentzVector operator+(const LorentzVector& i_firstLorentzVector,
                        const LorentzVector& i_secondLorentzVector)
{
    LorentzVector o_lorentzVector(i_firstLorentzVector);
    o_lorentzVector += i_secondLorentzVector;
    
    return o_lorentzVector;
}

LorentzVector operator-(const LorentzVector& i_firstLorentzVector,
                        const LorentzVector& i_secondLorentzVector)
{
    LorentzVector o_lorentzVector(i_firstLorentzVector);
    o_lorentzVector -= i_secondLorentzVector;
    
    return o_lorentzVector;
}
