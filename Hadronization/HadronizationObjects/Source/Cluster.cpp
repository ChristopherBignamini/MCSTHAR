#include "../Include/Cluster.h"
#include "../../../Utilities/Include/Constants.h"

Cluster::Cluster(void)
                :m_isHeavyFlavored(false)
                ,m_isHadronizable(true)
                ,m_energyDensity(0.)
                ,m_parentIndexes(0,0)
{
    m_idCode = clusterIdCode;
}


Cluster::Cluster(const bool i_isHadronizable,
                 const pair<unsigned int,unsigned int>& i_parentIndexes,
                 const int i_strangeCharge,
                 const int i_charmCharge,
                 const int i_bottomCharge,
                 const double i_electricCharge,
                 const double i_baryonicCharge,
                 const double i_energyDensity,
                 const TLorentzVector& i_p,
                 const unsigned int i_index)
                :HadronizationObject(clusterIdCode,
                                     i_strangeCharge,
                                     i_charmCharge,
                                     i_bottomCharge,
                                     i_electricCharge,
                                     i_baryonicCharge,
                                     i_p,
                                     i_index)
                ,m_isHadronizable(i_isHadronizable)
                ,m_energyDensity(i_energyDensity)
                ,m_parentIndexes(i_parentIndexes.first,i_parentIndexes.second)
{
    if(i_charmCharge || i_bottomCharge)
    {
        m_isHeavyFlavored = true;
    }
    else
    {
        m_isHeavyFlavored = false;
    }
}


Cluster::Cluster(const bool i_isHadronizable,
                 const pair<unsigned int,unsigned int>& i_parentIndexes,
                 const int i_strangeCharge,
                 const int i_charmCharge,
                 const int i_bottomCharge,
                 const double i_electricCharge,
                 const double i_baryonicCharge,
                 const double i_energyDensity,
                 const double i_mass,
                 const unsigned int i_index)
                :HadronizationObject(clusterIdCode,
                                     i_strangeCharge,
                                     i_charmCharge,
                                     i_bottomCharge,
                                     i_electricCharge,
                                     i_baryonicCharge,
                                     i_mass,
                                     i_index)
                ,m_isHadronizable(i_isHadronizable)
                ,m_energyDensity(i_energyDensity)
                ,m_parentIndexes(i_parentIndexes.first,i_parentIndexes.second)
{
    if(i_charmCharge || i_bottomCharge)
    {
        m_isHeavyFlavored = true;
    }
    else
    {
        m_isHeavyFlavored = false;
    }
}


Cluster::~Cluster(void)
{
}

void Cluster::setClusterData(const bool i_isHadronizable,
                             const pair<unsigned int,unsigned int>& i_parentIndexes,
                             const int i_strangeCharge,
                             const int i_charmCharge,
                             const int i_bottomCharge,
                             const double i_electricCharge,
                             const double i_baryonicCharge,
                             const double i_energyDensity,
                             const TLorentzVector& i_p)
{
	m_isHadronizable = i_isHadronizable;
    m_parentIndexes = i_parentIndexes;
	m_strangeCharge = i_strangeCharge;
	m_charmCharge = i_charmCharge;
	m_bottomCharge = i_bottomCharge;
	m_electricCharge = i_electricCharge;
	m_baryonicCharge = i_baryonicCharge;
	m_energyDensity = i_energyDensity;
    m_p = i_p;
    
    if(i_charmCharge || i_bottomCharge)
    {
        m_isHeavyFlavored = true;
    }
    else
    {
        m_isHeavyFlavored = false;
    }
}

Cluster operator+(const Cluster& i_firstCluster,const Cluster& i_secondCluster)
{
    
    Cluster o_sumCluster(true,
                         pair<unsigned int,unsigned int>(i_firstCluster.getIndex(),i_secondCluster.getIndex()),
                         i_firstCluster.getStrangeCharge()+i_secondCluster.getStrangeCharge(),
                         i_firstCluster.getCharmCharge()+i_secondCluster.getCharmCharge(),
                         i_firstCluster.getBottomCharge()+i_secondCluster.getBottomCharge(),
                         i_firstCluster.getElectricCharge()+i_secondCluster.getElectricCharge(),
                         i_firstCluster.getBaryonicCharge()+i_secondCluster.getBaryonicCharge(),
                         i_firstCluster.getEnergyDensity(),
                         i_firstCluster.getP()+i_secondCluster.getP());
    
    return o_sumCluster;
}