#ifndef CONVERTPARTICLEID_H
#define CONVERTPARTICLEID_H

/**
* @brief Fortran CONVERTPARTICLEID routine wrapper header
*
* @author Christopher Bignamini
*/
#define convertParticleId convertparticleid_
extern "C"
{
    /**
    * @brief Particle ID conversion subroutine
    *
    * Subroutine used to convert hadron IDs used by
    * MCSTHAR++ into those adopted by Herwig6510
    *
    * @author Christopher Bignamini
    */
	void convertParticleId(void);
}

#endif