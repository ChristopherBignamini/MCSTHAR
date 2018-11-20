#include "../Include/HadronizationObject.h"


HadronizationObject::HadronizationObject(void)
                                        :m_idCode(0)
                                        ,m_objectIndex(0)
                                        ,m_strangeCharge(0)
                                        ,m_charmCharge(0)
                                        ,m_bottomCharge(0)
                                        ,m_electricCharge(0.)
                                        ,m_baryonicCharge(0.)
                                        ,m_mass(0.)
                                        ,m_p(0.,0.,0.,0.)
{
}

HadronizationObject::HadronizationObject(const int i_idCode,
                                         const int i_strangeCharge,
                                         const int i_charmCharge,
                                         const int i_bottomCharge,
                                         const double i_electricCharge,
                                         const double i_baryonicCharge,
                                         const TLorentzVector& i_p,
                                         const unsigned int i_index)
                                        :m_idCode(i_idCode)
                                        ,m_objectIndex(i_index)
                                        ,m_strangeCharge(i_strangeCharge)
                                        ,m_charmCharge(i_charmCharge)
                                        ,m_bottomCharge(i_bottomCharge)
                                        ,m_electricCharge(i_electricCharge)
                                        ,m_baryonicCharge(i_baryonicCharge)
                                        ,m_mass(i_p.Mag())
                                        ,m_p(i_p)
{
}

HadronizationObject::HadronizationObject(const int i_idCode,
                                         const int i_strangeCharge,
                                         const int i_charmCharge,
                                         const int i_bottomCharge,
                                         const double i_electricCharge,
                                         const double i_baryonicCharge,
                                         const double i_mass,
                                         const unsigned int i_index)
                                        :m_idCode(i_idCode)
                                        ,m_objectIndex(i_index)
                                        ,m_strangeCharge(i_strangeCharge)
                                        ,m_charmCharge(i_charmCharge)
                                        ,m_bottomCharge(i_bottomCharge)
                                        ,m_electricCharge(i_electricCharge)
                                        ,m_baryonicCharge(i_baryonicCharge)
                                        ,m_mass(i_mass)
                                        ,m_p(TLorentzVector(0.,0.,0.,i_mass))
{
}


HadronizationObject::~HadronizationObject(void)
{
}

