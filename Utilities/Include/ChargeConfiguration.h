#ifndef CHARGECONFIGURATION_H
#define CHARGECONFIGURATION_H

namespace MCSTHAR
{
    namespace Utilities
    {
        
        /**
        * @brief Charge configuration storage structure
        *
        * @author Christopher Bignamini
        */
        struct ChargeConfiguration
        {
            /**
            * @brief Constructor
            */
            ChargeConfiguration(void)
                               :strangeCharge(0)
                               ,charmCharge(0)
                               ,bottomCharge(0)
                               ,electricCharge(0.)
                               ,baryonicCharge(0.)
            {
            }
            
            /**
            * @brief Constructor
            *
            * @param i_strangeCharge Strange charge
            * @param i_charmCharge Charm charge
            * @param i_bottomCharge Bottom charge
            * @param i_electricCharge Electric charge
            * @param i_baryonicCharge Baryonic charge
            */
            ChargeConfiguration(const int i_strangeCharge,
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
            
            /**
            * @brief Anti-charge configuration get method
            */
            inline ChargeConfiguration getAntiConfiguration(void) const
            {
                return ChargeConfiguration(-strangeCharge,-charmCharge,-bottomCharge,-electricCharge,-baryonicCharge);
            }
            
            /**
            * Strange charge
            */
            int strangeCharge;
            
            /**
            * Charm charge
            */
            int charmCharge;
            
            /**
            * Bottom charge
            */
            int bottomCharge;
            
            /**
            * Electric charge
            */
            double electricCharge;
            
            /**
            * Baryonic charge
            */
            double baryonicCharge;
            
            /**
            * @brief ChargeConfiguration comparison operator
            *
            * @param i_chargeConfiguration Charge configuration used for comparison
            * @return True if *this is less than i_chargeConfiguration (see ordering in method implementation), false otherwise
            */
            bool operator<(const ChargeConfiguration& i_chargeConfiguration) const;
            
        };

    }
}
    
#endif
