#include "../Include/HadronData.h"
#include "../../../Utilities/Include/HadronizationException.h"

HadronData::HadronData(void)
                      :m_idCode(0)
                      ,m_strangeCharge(0)
                      ,m_charmCharge(0)
                      ,m_bottomCharge(0)
                      ,m_electricCharge(0.)
                      ,m_baryonicCharge(0.)
                      ,m_mass(1.)
                      ,m_spinMultiplicity(0)
                      ,m_strangeQuarkNumber(0)
                      ,m_ssbarComponentValue(0.)
                      ,m_oneMinusSSBARComponentValue(1.)
{
}

HadronData::HadronData(const int i_idCode,
                       const int i_strangeCharge,
                       const int i_charmCharge,
                       const int i_bottomCharge,
                       const double i_electricCharge,
                       const double i_baryonicCharge,
                       const double i_mass,
                       const unsigned int i_spinMultiplicity,
                       const unsigned int i_strangeQuarkNumber,
                       const double i_ssbarComponentValue)
                      :m_idCode(i_idCode)
                      ,m_strangeCharge(i_strangeCharge)
                      ,m_charmCharge(i_charmCharge)
                      ,m_bottomCharge(i_bottomCharge)
                      ,m_electricCharge(i_electricCharge)
                      ,m_baryonicCharge(i_baryonicCharge)
                      ,m_mass(i_mass)
                      ,m_spinMultiplicity(i_spinMultiplicity)
                      ,m_strangeQuarkNumber(i_strangeQuarkNumber)
                      ,m_ssbarComponentValue(i_ssbarComponentValue)
                      ,m_oneMinusSSBARComponentValue(1.-i_ssbarComponentValue)
{
    if(i_ssbarComponentValue<0. || i_ssbarComponentValue>1.)
    {
        throw HadronizationException("Error during hadron data setting, invalid ssbar wave function component provided",
                                     __FUNCTION__,221);
    }
    
    if(i_mass<=0.)
    {
        throw HadronizationException("Error during hadron data setting, invalid mass value provided",
                                     __FUNCTION__,222);
    }
}
    
HadronData::~HadronData(void)
{
}


void HadronData::setMass(const double i_mass)
{
    if(i_mass<=0.)
    {
        throw HadronizationException("Error during hadron data setting, invalid mass value provided",
                                     __FUNCTION__,222);
    }
    
    m_mass = i_mass;
}

void HadronData::setSSBARComponentValue(const double i_ssbarComponentValue)
{
    if(i_ssbarComponentValue<0. || i_ssbarComponentValue>1.)
    {
        throw HadronizationException("Error during hadron data setting, invalid ssbar wave function component provided",
                                     __FUNCTION__,221);
    }

    m_ssbarComponentValue = i_ssbarComponentValue;
    m_oneMinusSSBARComponentValue = 1. - i_ssbarComponentValue;
}
