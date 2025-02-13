#ifndef HWCFORMCSTHAR_H
#define HWCFORMCSTHAR_H

/**
* @brief Fortran HWCFORMCSTHAR subroutine wrapper header
*
* @author Christopher Bignamini (wrapper only)
*/
#define hwcforMCSTHAR hwcformcsthar_
extern "C"
{
    /**
    * @brief Herwig6510 cluster formation subroutine
    *
    * Herwig6510 cluster formation subroutine modified
    * to avoid single cluster decay. The modified
    * parts are identified in .f file by the CCC MCSTHAR++ CCC 
    * tag.
    *
    * @author Christopher Bignamini (wrapper only)
    */
	void hwcforMCSTHAR(void);
}

#endif