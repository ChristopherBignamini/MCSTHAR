#ifndef CHARGEDATA_H
#define CHARGEDATA_H

/**
* @brief Herwig cluster charge configuration structure
*
* @author Christopher Bignamini
*/
struct ChargeData
{
    /**
    * @brief Constructor
    *
    * @param i_strangeCharge Strange charge value
    * @param i_charmCharge Charm charge value
    * @param i_bottomCharge Bottom charge value
    * @param i_electricCharge Electric charge value
    * @param i_baryonicCharge Baryonic charge value
    */
    ChargeData(const int i_strangeCharge,
               const int i_charmCharge,
               const int i_bottomCharge,
               const double i_electricCharge,
               const double i_baryonicCharge)
              :strangeCharge(i_strangeCharge)
              ,charmCharge(i_charmCharge)
              ,bottomCharge(i_bottomCharge)
              ,electricCharge(i_electricCharge)
              ,baryonicCharge(i_baryonicCharge)
    {
    }
    
    // Default copy constructor and assignement operator are being used
    
    /**
    * Strange charge value
    */
    int strangeCharge;

    /**
    * Charm charge value
    */
    int charmCharge;
    
    /**
    * Bottom charge value
    */
    int bottomCharge;
    
    /**
    * Electric charge value
    */
    double electricCharge;
    
    /**
    * Baryonic charge value
    */
    double baryonicCharge;
};

#endif
