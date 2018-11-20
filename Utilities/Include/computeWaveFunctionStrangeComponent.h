#ifndef COMPUTEWAVEFUNCTIONSTRANGECOMPONENT_H
#define COMPUTEWAVEFUNCTIONSTRANGECOMPONENT_H

/**
* @brief Meson series enumerator
* 
* @author Christopher Bignamini
*/
enum MesonSeries
{
    Octect,
    Singlet
};

/**
* @brief Wave function strange component square module calculation function
*
*  This function performs the calculation of the square module of the strange component
*  of the light unflavoured mesons: given the wave function of the form
*
*          Cu*|uubar> + Cd*|ddbar> + Cs*|ssbar>
*
*  this function returns the |Cs|^2 value.
*
* @param i_mesonType Meson series (octect or singlet)
* @param i_mixingAngle Mixing angle (degree)
* @return Wave function strange component square modulus
* @throw HadronizationException in case of not identified meson series 
*
*
* @author Christopher Bignamini
*/
double computeWaveFunctionStrangeComponent(MesonSeries i_mesonSeries,double i_mixingAngle);

#endif