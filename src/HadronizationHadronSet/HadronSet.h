#ifndef HADRONSET_H
#define HADRONSET_H

#include <vector>
#include <string>
#include "HadronData.h"

using namespace std;

/**
* Hadron type list
*/
enum HadronType
{
    Meson,
    Baryon
};

/**
* @brief Hadron set storing class
*
* Class used to store the set of hadrons used during the
* hadronization process. The list of hadrons is loaded
* from a user specified text file.
*
* NOTE: - K_short and K_long are excluded from the set of
*         hadrons during the hadron list loading process.
*       - Light flavored hadrons with mass value larger than 
*         i_lightHadronMaxMass are excluded from the hadron set.
*       - Hadrons containing more than one heavy flavored 
*         constituent are excluded from the hadron set.
*       - Hadrons with a top quark constituent are excluded 
*         from the hadron set.
*       - Only two and three quark states are identified
*         and therefore loaded.
*
* @author Christopher Bignamini
*/
class HadronSet
{
    
    public:
        
        /**
        * @brief Constructor
        *
        */
        HadronSet(void);
    
        // TODO: use new constructor without data file, and with hadron data vector
        // TODO: switch to xml for hadron data
        /**
        * @brief Constructor
        *
        * @param i_hadronListFile Hadron list file
        * @param i_lightHadronMaxMass Light flavored hadron maximum mass value (hadron with larger mass are excluded)
        * @throw HadronizationException if an error occurs during hadron set loading (not existing file/empty hadron set)
        */
        HadronSet(const string& i_hadronListFile,double i_lightHadronMaxMass);

        /**
        * @brief Destructor
        *
        */
        virtual ~HadronSet(void);

        // Default copy constructor and overloaded assignement operator are being used
    
        /**
        * @brief Hadron set building method
        *
        * This method builds the meson and baryon set starting
        * from the provided input hadron list file. A call to this 
        * method causes the deletion of a previously built hadron set.
        *
        * @param i_hadronListFile Hadron list file
        * @param i_lightHadronMaxMass Light flavored hadron maximum mass value (hadron with larger mass are excluded)
        * @throw HadronizationException if an error occurs during hadron set loading (not existing file/empty hadron set)
        */
        void buildHadronSet(const string& i_hadronListFile,double i_lightHadronMaxMass);

        /**
        * @brief Hadron data getter from id code
        *
        * @param i_idCode Hadron id code
        * @return Hadron data corresponding to the required hadron
        * @throw HadronizationException if no hadron with the provided id code exists in the built hadron set
        */
        const HadronData& getHadronFromIdCode(int i_idCode) const;
            
        /**
        * @brief Hadron group (mesons or baryons) getter
        *
        * This method returns the loaded set of mesons or baryons
        *
        * @param i_hadronType Hadron type (Meson or Baryon)
        * @return Hadron data set corresponding to the provided hadron type
        * @throw HadronizationException in case of not identified hadron type
        */
        const vector<HadronData>& getHadronGroup(HadronType i_hadronType) const;
        
        /**
        * @brief Number of mesons getter
        *
        * @return Number of loaded mesons
        */
        inline unsigned int getNumberOfMesons(void) const { return m_mesonSet.size(); }

        /**
        * @brief Number of baryons getter
        *
        * @return Number of loaded baryons
        */
        inline unsigned int getNumberOfBaryons(void) const { return m_baryonSet.size(); }

    protected:
        
        /**
        * Loaded meson set
        */
        vector<HadronData> m_mesonSet;
        // TODO: switch to map for the final version (vector are needed for comparison with old code)
        // map<int,HadronData> m_mesonSet;
    
        /**
        * Loaded baryon set
        */
        vector<HadronData> m_baryonSet;
        // TODO: switch to map for the final version (vector are needed for comparison with old code)
        // map<int,HadronData> m_baryonSet;        
};

#endif
