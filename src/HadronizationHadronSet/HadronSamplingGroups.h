#ifndef HADRONSAMPLINGGROUPS_H
#define HADRONSAMPLINGGROUPS_H

#include "HadronSet.h"
#include <map>

using namespace std;

/**
* Hadron sampling group list
*/
enum HadronGroup
{
    LightBaryons,
    LightAntibaryons,
    LightStrangeMesons,
    LightStrangeAntimesons,
    LightChargedMesons,
    LightChargedAntimesons,
    LightNeutralMesons,
    CharmedHadrons,
    AntiCharmedHadrons,
    BottomedHadrons,
    AntiBottomedHadrons
};

/**
* @brief Flavor ordered/grouped hadron set storing class
*
* Class inheriting from base HadronSet class and used to introduce a 
* flavor dependent hadron grouping and ordering as required by the 
* currently adopted hadronization channel sampling algorithm.
* Hadrons stored in the hadron set are divided into the groups listed
* in the HadronGroup enumerator
*
* @author Christopher Bignamini
*/
// TODO: remove inheritance???
class HadronSamplingGroups : public HadronSet
{
    
    public:
    
        /**
        * @brief Constructor
        *
        */
        HadronSamplingGroups(void);
        
        /**
        * @brief Constructor
        *
        * @param i_hadronListFile Hadron list file
        * @param i_lightHadronMaxMass Light flavored hadron maximum mass value (hadron with larger mass are excluded)
        * @throw HadronizationException if an error occurs during hadron set loading (not existing file/empty hadron set) or ordering/grouping
        */
        HadronSamplingGroups(const string& i_hadronListFile,double i_lightHadronMaxMass);
        
        /**
        * @brief Destructor
        *
        */
        ~HadronSamplingGroups(void);
        
        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Hadron sampling groups builder
        *
        * Method building the sampling meson and baryon groups starting
        * from the provided input hadron list file. A call to this
        * method causes the deletion of previously built hadron groups.
        *
        * @param i_hadronListFile Hadron list file
        * @param i_lightHadronMaxMass Light flavored hadron maximum mass value (hadron with larger mass are excluded)
        * @throw HadronizationException if an error occurs during hadron set loading (not existing file/empty hadron set) or ordering/grouping
        */
        void buildHadronSamplingGroups(const string& i_hadronListFile,double i_lightHadronMaxMass);
        
        /**
        * @brief Hadron sampling group getter
        *
        * @param i_hadronSamplingGroup Hadron sampling group 
        * @return Vector of HadronData(s) corresponding to the hadrons of the specified sampling group
        * @throw HadronizationException in case of unidentified hadron group
        */
        // TODO: add constness to this method
        const vector<const HadronData*>& getHadronGroup(HadronGroup i_hadronSamplingGroup) const;

        /**
        * @brief Hadron sampling group list getter
        *
        * @return List of available (not empty) hadron groups
        * @throw HadronizationException if no hadron group is available
        */
        vector<HadronGroup> getHadronGroupList(void) const;

        /**
        * @brief Hadron sampling group multiplicity getter
        *
        * @return Number of hadrons belonging to the provided sampling group
        * @throw HadronizationException in case of unidentified hadron group
        */
        unsigned int getNumberOfHadronsInGroup(HadronGroup i_hadronSamplingGroup) const;
    
    private:

        /**
        * @brief Hadron sampling groups building method
        *
        * Method building the sampling meson and baryon groups 
        * (private method called by constructors and buildHadronSamplingGroups public method)
        *
        * @throw HadronizationException if an error occurs during hadron set loading (not existing file/empty hadron set) or ordering/grouping
        */
        void buildHadronSamplingGroups(void);
    
        /**
        * Hadron sampling group map
        */
        map<HadronGroup, vector<const HadronData*> > m_hadronSamplingGroups;
    
};

#endif
