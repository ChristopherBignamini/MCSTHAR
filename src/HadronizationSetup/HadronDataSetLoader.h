#ifndef HADRONDATASETLOADER_H
#define HADRONDATASETLOADER_H

#include "PartitionFunctionCalculationSetup.h"
#include "../../HadronizationHadronSet/Include/HadronData.h"// TODO: move class to new location
#include <vector>
#include <string>

using namespace std;

/**
* @brief Hadron data set file loading class
*
* @author Christopher Bignamini
*/
class HadronDataSetLoader
{
    public:
        
        /**
        * @brief Constructor
        *
        * @param i_hadronDataSetFileName Hadron data set file name
        * @param i_lightHadronMaxMass Light flavored hadron maximum mass value (hadron with larger mass are excluded)
        * @throw HadronizationException in case of error during data set file loading or empty data set
        */
        HadronDataSetLoader(const string& i_hadronDataSetFileName,
                            double i_lightHadronMaxMass);
        
        /**
        * @brief Destructor
        */
        ~HadronDataSetLoader(void);
        
        // Default copy constructor and assignement operator are being used
        
        /**
        * @brief Hadron data setp data get method
        *
        * @return Loaded hadron data set
        * @throw HadronizationException in case of hadron data set unavailability
        */
        const vector<HadronData>& getHadronDataSet(void) const;
        
    private:
    
    
        /**
        * @brief Hadron abelian charge calculation method
        *
        * This method computes the abelian charges of a hadron
        * starting from its flavor composition and id code
        *
        * @param i_flavor Hadron flavor composition
        * @param i_pid Hadron id code
        * @param io_strangeCharge Hadron electric charge (I/O parameter)
        * @param io_strangeCharge Hadron strange charge (I/O parameter)
        * @param io_charmCharge Hadron charm charge (I/O parameter)
        * @param io_bottomCharge Hadron bottom charge (I/O parameter)
        * @param io_topCharge Hadron top charge (I/O parameter)
        * @param io_heavyFlavoredConstituentNumber Number of heavy flavored hadron constituents (I/O parameter)
        * @throw HadronizationException in case of not identified flavor composition
        */
        void computeHadronCharges(unsigned int i_flavor,
                                  int i_pid,
                                  int& io_electricCharge,
                                  int& io_strangeCharge,
                                  int& io_charmCharge,
                                  int& io_bottomCharge,
                                  int& io_topCharge,
                                  unsigned int& io_heavyFlavoredConstituentNumber) const;

        /**
        * Hadron data set availability status
        */
        bool m_hadronDataSetAvailability;
        
        /**
        * Hadron data set
        */
        vector<HadronData> m_hadronDataSet;
};

#endif
