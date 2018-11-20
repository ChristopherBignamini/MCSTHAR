#include "../Include/HadronSet.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../HadronizationSetup/Include/HadronDataSetLoader.h"

#include <fstream>
using namespace std;

HadronSet::HadronSet(void)
{
}

HadronSet::HadronSet(const string& i_hadronListFile,const double i_lightHadronMaxMass)
{
    try
    {
        buildHadronSet(i_hadronListFile,i_lightHadronMaxMass);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

HadronSet::~HadronSet(void)
{
}

void HadronSet::buildHadronSet(const string& i_hadronListFile, const double i_lightHadronMaxMass)
{
    try
    {
        HadronDataSetLoader hadronDataSetLoader(i_hadronListFile,
                                                i_lightHadronMaxMass);

        ofstream outFile;
        outFile.open("blaBla.txt");

        m_mesonSet.clear();
        m_baryonSet.clear();
        for(unsigned int hadronIndex=0;hadronIndex<hadronDataSetLoader.getHadronDataSet().size();++hadronIndex)
        {
            if(hadronDataSetLoader.getHadronDataSet()[hadronIndex].getBaryonicCharge()!=0)
            {
                m_baryonSet.push_back(hadronDataSetLoader.getHadronDataSet()[hadronIndex]);
            }
            else
            {
                m_mesonSet.push_back(hadronDataSetLoader.getHadronDataSet()[hadronIndex]);        
            }
            outFile<<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getIdCode()<<" "
                <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getMass()<<" "
                <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getSpinMultiplicity()<<" "
                <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getElectricCharge()<<" "
                <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getStrangeCharge()<<" "
                <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getBaryonicCharge()<<" "
                <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getCharmCharge()<<" "
            <<hadronDataSetLoader.getHadronDataSet()[hadronIndex].getBottomCharge()<<endl;
        }
        
        outFile.close();
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    return;
}

const HadronData& HadronSet::getHadronFromIdCode(const int i_idCode) const
{
    
    // Search the hadron in the meson set
    // TODO: switch to map for the final version (vector are needed for comparison with old code)
//    map<int,HadronData>::const_iterator mapIt = m_mesonSet.find(i_idCode);
//    if(mapIt!=m_mesonSet.end())
//    {
//        return Particle(mapIt->second);
//    }
//        
//    // Search the hadron in the baryon set
//    mapIt = m_baryonSet.find(i_idCode);
//    if(mapIt!=m_baryonSet.end())
//    {
//        return Particle(mapIt->second);
//    }
    for(unsigned int hadronIndex=0;hadronIndex<m_mesonSet.size();++hadronIndex)
    {
        if(m_mesonSet[hadronIndex].getIdCode() == i_idCode)
        {
            return m_mesonSet[hadronIndex];
        }
    }

    for(unsigned int hadronIndex=0;hadronIndex<m_baryonSet.size();++hadronIndex)
    {
        if(m_baryonSet[hadronIndex].getIdCode() == i_idCode)
        {
            return m_baryonSet[hadronIndex];
        }
    }

    
    // Hadron not found, throw exception
    throw HadronizationException("Error during hadron search, ID not identified",
                                 __FUNCTION__,201);
        
}

const vector<HadronData>& HadronSet::getHadronGroup(const HadronType i_hadronType) const
{
//    vector<Particle> o_hadronGroup;
    // TODO: switch to map for the final version (vector are needed for comparison with old code)
//    map<int,HadronData>::const_iterator mapIt;
//	if(i_hadronType=='m')
//    {
//        for(mapIt = m_mesonSet.begin();mapIt!=m_mesonSet.end();mapIt++)
//        {
//            o_hadronGroup.push_back(Particle(mapIt->second));
//        }
//		return o_hadronGroup;			
//	}
//	else if(i_hadronType=='b')
//    {
//        for(mapIt = m_baryonSet.begin();mapIt!=m_baryonSet.end();mapIt++)
//        {
//            o_hadronGroup.push_back(Particle(mapIt->second));
//        }
//		return o_hadronGroup;
//	}

	if(i_hadronType==Meson)
    {
        return m_mesonSet;
	}
	else if(i_hadronType==Baryon)
    {
        return m_baryonSet;
	}
    
    // Hadron group not found, throw exception
    throw HadronizationException("Error during hadron group search, group not identified",
                                 __FUNCTION__,202);
}


