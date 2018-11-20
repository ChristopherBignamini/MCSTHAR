#include "../Include/Particle.h"

Particle::Particle(void)
                  :m_spinMultiplicity(0)
                  ,m_ssbarComponentValue(0.)
{
}


Particle::Particle(const HadronData& i_hadronData,const unsigned int i_parentIndex,const unsigned int i_particleIndex)
                  :HadronizationObject(i_hadronData.getIdCode(),
                                       i_hadronData.getStrangeCharge(),
                                       i_hadronData.getCharmCharge(),
                                       i_hadronData.getBottomCharge(),
                                       i_hadronData.getElectricCharge(),
                                       i_hadronData.getBaryonicCharge(),
                                       i_hadronData.getMass(),
                                       i_particleIndex)
                  ,m_spinMultiplicity(i_hadronData.getSpinMultiplicity())
                  ,m_ssbarComponentValue(i_hadronData.getSSBARComponentValue())
                  ,m_parentIndex(i_parentIndex)
{
}

void Particle::setData(const HadronData& i_hadronData)
{
    m_idCode = i_hadronData.getIdCode();
    m_strangeCharge = i_hadronData.getStrangeCharge();
    m_charmCharge = i_hadronData.getCharmCharge();
    m_bottomCharge = i_hadronData.getBottomCharge();
    m_electricCharge = i_hadronData.getElectricCharge();
    m_baryonicCharge = i_hadronData.getBaryonicCharge();
    m_mass = i_hadronData.getMass();
    m_spinMultiplicity = i_hadronData.getSpinMultiplicity();
    m_ssbarComponentValue = i_hadronData.getSSBARComponentValue();
    m_p.SetPxPyPzE(0.,0.,0.,m_mass);
}

void Particle::setDataAndMomentum(const HadronData& i_hadronData,const TLorentzVector& i_momentum)
{
    setData(i_hadronData);
    m_p = i_momentum;
}

Particle::~Particle(void)
{
}




