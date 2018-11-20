#include "../Include/HadronSamplingGroups.h"
#include "../../../Utilities/Include/HadronizationException.h"

HadronSamplingGroups::HadronSamplingGroups(void)
{
}

HadronSamplingGroups::HadronSamplingGroups(const string& i_hadronListFile,const double i_lightHadronMaxMass)
                                          :HadronSet(i_hadronListFile,i_lightHadronMaxMass)
{
    try
    {
        buildHadronSamplingGroups();
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

HadronSamplingGroups::~HadronSamplingGroups(void)
{
}
    
void HadronSamplingGroups::buildHadronSamplingGroups(const string& i_hadronListFile,const double i_lightHadronMaxMass)
{    
    try
    {
        buildHadronSet(i_hadronListFile,i_lightHadronMaxMass);
        buildHadronSamplingGroups();
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}
    
const vector<const HadronData*>& HadronSamplingGroups::getHadronGroup(const HadronGroup i_hadronSamplingGroup) const
{
    map<HadronGroup, vector<const HadronData*> >::const_iterator samplingGroup =
        m_hadronSamplingGroups.find(i_hadronSamplingGroup);
    
    if(samplingGroup != m_hadronSamplingGroups.end())
    {
        return samplingGroup->second;
    }
    else
    {
        // Not identified hadron group, throw exception
        throw HadronizationException("Error during hadron sampling groups retrieval, not identified hadron sampling group",
                                     __FUNCTION__,215);
    }
}

vector<HadronGroup> HadronSamplingGroups::getHadronGroupList(void) const
{
    if(m_hadronSamplingGroups.size()>0)
    {
        vector<HadronGroup> o_hadronGroupList(m_hadronSamplingGroups.size());
        unsigned int hadronGroupIndex = 0;
        for(map<HadronGroup, vector<const HadronData*> >::const_iterator hadronGroup = m_hadronSamplingGroups.begin();
            hadronGroup != m_hadronSamplingGroups.end();++hadronGroup)
        {
            o_hadronGroupList[hadronGroupIndex] = hadronGroup->first;
            ++hadronGroupIndex;
        }
        return o_hadronGroupList;
    }
    else
    {
        // No available hadron group, throw exception
        throw HadronizationException("Error during hadron sampling group list retrieval, no available groups",
                                     __FUNCTION__,216);
    }
}

unsigned int HadronSamplingGroups::getNumberOfHadronsInGroup(const HadronGroup i_hadronSamplingGroup) const
{
    try
    {
        return (getHadronGroup(i_hadronSamplingGroup)).size();
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void HadronSamplingGroups::buildHadronSamplingGroups(void)
{
    m_hadronSamplingGroups.clear();
    
    try
    {        
        HadronGroup samplingGroup;
        
        // Build baryonic groups
        for(vector<HadronData>::iterator m_barIt = m_baryonSet.begin(); m_barIt != m_baryonSet.end();
            ++m_barIt)
        {
            
            if(m_barIt->getBaryonicCharge()==0)
            {
                // Meson found in baryon set, throw exception
                throw HadronizationException("Error during hadron sampling groups building, meson wrongly included in baryon set",
                                             __FUNCTION__,213);
            }
            
            if(m_barIt->getBaryonicCharge()==1. && m_barIt->getCharmCharge()==0 && m_barIt->getBottomCharge()==0)
            {
                samplingGroup = LightBaryons;
            }
            else if(m_barIt->getBaryonicCharge()==-1. && m_barIt->getCharmCharge()==0 && m_barIt->getBottomCharge()==0)
            {
                samplingGroup = LightAntibaryons;
            }
            else if(m_barIt->getCharmCharge()>0)
            {
                samplingGroup = CharmedHadrons;
            }
            else if(m_barIt->getCharmCharge()<0)
            {
                samplingGroup = AntiCharmedHadrons;
            }
            else if(m_barIt->getBottomCharge()<0)
            {
                samplingGroup = BottomedHadrons;
            }
            else if(m_barIt->getBottomCharge()>0)
            {
                samplingGroup = AntiBottomedHadrons;
            }
            else
            {
                // Not identified hadron type, throw exception
                throw HadronizationException("Error during hadron sampling groups building, not identified hadron group",
                                             __FUNCTION__,211);
            }
                        
            m_hadronSamplingGroups[samplingGroup].push_back(&(*m_barIt));
        }
        
        
        // Build mesonic groups
        for(vector<HadronData>::iterator m_mesIt = m_mesonSet.begin(); m_mesIt != m_mesonSet.end();
            ++m_mesIt)
        {
            if(m_mesIt->getBaryonicCharge())
            {
                // Baryons found in meson set, throw exception
                throw HadronizationException("Error during hadron sampling groups building, baryon wrongly included in meson set",
                                             __FUNCTION__,214);
            }
            
            if(m_mesIt->getStrangeCharge()==1 && m_mesIt->getCharmCharge()==0 && m_mesIt->getBottomCharge()==0)
            {
                samplingGroup = LightStrangeMesons;// TODO: these two groups must be inverted (and sampling procedure accordingly corrected)
            }
            else if(m_mesIt->getStrangeCharge()==-1 && m_mesIt->getCharmCharge()==0 && m_mesIt->getBottomCharge()==0)
            {
                samplingGroup = LightStrangeAntimesons;// TODO: these two groups must be inverted (and sampling procedure accordingly corrected)
            }
            else if(m_mesIt->getStrangeCharge()==0 && m_mesIt->getElectricCharge()>0. && m_mesIt->getCharmCharge()==0 &&
                    m_mesIt->getBottomCharge()==0)
            {
                samplingGroup = LightChargedMesons;
            }
            else if(m_mesIt->getStrangeCharge()==0 && m_mesIt->getElectricCharge()<0. && m_mesIt->getCharmCharge()==0 &&
                    m_mesIt->getBottomCharge()==0)
            {
                samplingGroup = LightChargedAntimesons;
            }
            else if(m_mesIt->getStrangeCharge()==0 && m_mesIt->getElectricCharge()==0. && m_mesIt->getCharmCharge()==0 &&
                    m_mesIt->getBottomCharge()==0)
            {
                samplingGroup = LightNeutralMesons;
            }
            else if(m_mesIt->getCharmCharge()>0)
            {
                samplingGroup = CharmedHadrons;
            }
            else if(m_mesIt->getCharmCharge()<0)
            {
                samplingGroup = AntiCharmedHadrons;
            }
            else if(m_mesIt->getBottomCharge()<0)
            {
                samplingGroup = BottomedHadrons;
            }
            else if(m_mesIt->getBottomCharge()>0)
            {
                samplingGroup = AntiBottomedHadrons;
            }
            else
            {
                // Not identified hadron type, throw exception
                throw HadronizationException("Error during hadron sampling groups building, not identified hadron group",
                                             __FUNCTION__,212);
            }
            
            m_hadronSamplingGroups[samplingGroup].push_back(&(*m_mesIt));
        }        
                        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

