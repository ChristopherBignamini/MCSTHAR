#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

/**
* @brief Constant parameter header
*
* @author Christopher Bignamini
*/

// Physical/Mathematical constants
// TODO: use MATH_PI???
/**
* pi
*/
static const double pi(3.1415926536);

/**
* pi^2
*/
static const double piSquared(9.86960438);

/**
* pi/180
*/
static const double piOn180(pi/180.);

/**
* 2*pi
*/
static const double twoPi(2.*pi);

/**
* 2*pi^2
*/
static const double twoPiSquared(2.*piSquared);

/**
* (2*pi)^3
*/
static const double twoPiCube(pow(2.*pi,3));

/**
* 4*pi
*/
static const double fourPi(4.*pi);

/**
* Two times neutral pion mass
*/
static const double twoPiZeroMass(2.*0.13498e0);

/**
* GeV^3 to fm^3 volume conversion factor
*/
static const double gevCubeToFermiCube((1./0.197326968)*(1./0.197326968)*(1./0.197326968));

/**
* GeV^3 to fm^3 volume conversion factor/(2*pi)^3
*/
static const double gevCubeToFermiCubeOnTwoPiCube(gevCubeToFermiCube/twoPiCube);

/**
* Flavour basis constant mixing angle
*/
static const double constantMixingAngle(std::atan(std::sqrt(2.)));

/**
* Quark/diquark electric charges
*/
static const double oneThird(1./3.);
static const double twoThird(2./3.);
static const double fourThird(4./3.);


// TODO: separate from physical consts???
// MCSTHAR++ constants

/**
* Cluster ID code
*/
static const int clusterIdCode(91);

/**
* Tolerance parameter used in phase space generation for squared
* mass comparison with zero
*/
static const double phaseSpaceSquaredMassTolerance(1.e-13);

/**
* Tolerance parameter used in phase space generation for 3-momentum
* module comparison with zero
*/
static const double phaseSpace3MomentumModTolerance(1.e-13);

/**
* Tolerance parameter used in phase space generation/integration for
* energy availability: this parameter is used to accept a decay channel
* in terms of phase space availability. If the cluster mass is larger than
* the particle mass sum of a quantity smaller or equal to this parameter
* the decay channel itself is not accepted (The phase space integral in the partition
* function calculation and the phase space sampling weigth during event generation
* would be close to zero anyway)
*/
// TODO: as for other parameters, find final recipe
static const double phaseSpaceEnergyAvailabilityTolerance(1.e-13);

// TODO: is there a rigorous recipe for this? This, or even larger values,
// can be in principle be acceptable due to the very large partition function
// values observed...
/**
* A partition function squared error with abs value less or equal
* to this parameter is considered (and set) equal to zero
*/
static const double partitionFunctionSquaredErrorZeroRange(1.e-5);


#endif