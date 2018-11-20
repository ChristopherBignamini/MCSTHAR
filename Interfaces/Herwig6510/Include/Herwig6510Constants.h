#ifndef HERWIG6510CONSTANTS_H
#define HERWIG6510CONSTANTS_H

#include "ChargeData.h"
#include "../../../Utilities/Include/Constants.h"
#include <map>

using namespace std;

/**
* @brief Herwig6510 constants (see http://www.hep.phy.cam.ac.uk/theory/webber/Herwig/hw65_manual_rev.pdf )
*
* @author Christopher Bignamini 
*/

/**
* Herwig6501 cluster id code (IDPDG)
*/
static const int herwigClusterIdCode(91);

/**
* Herwig6501 cluster internal code (IDHW)
*/
static const int herwigInternalClusterIdCode(19);

/**
* Herwig6501 beam cluster before soft process status 
*/
static const int herwigBeamClusterBeforeSoftProcessStatus(164);

/**
* Herwig6501 target cluster before soft process status
*/
static const int herwigTargetClusterBeforeSoftProcessStatus(165);

/**
* Herwig6501 hard process cluster (hadronized) status
*/
static const int herwigHadronizedClusterStatus(183);

/**
* Herwig6501 direct unstable hadron (2-body clus.) before decay
*/
static const int herwigUnstableHadronFromClusterStatus(192);// TODO: check this (maybe indirect should be used for errors with herwig6510.f)


// TODO: check correctness of this strategy
/**
* Quark/Diquark to charge configurations map initialization (see static map below)
*/
static pair<int,ChargeData> quarkChargesPairs[] =
{
    make_pair(1,ChargeData(0,0,0,-oneThird,oneThird)),
    make_pair(-1,ChargeData(0,0,0,oneThird,-oneThird)),
    make_pair(2,ChargeData(0,0,0,twoThird,oneThird)),
    make_pair(-2,ChargeData(0,0,0,-twoThird,-oneThird)),
    make_pair(3,ChargeData(-1,0,0,-oneThird,oneThird)),
    make_pair(-3,ChargeData(1,0,0,oneThird,-oneThird)),
    make_pair(4,ChargeData(0,1,0,twoThird,oneThird)),
    make_pair(-4,ChargeData(0,-1,0,-twoThird,-oneThird)),
    make_pair(5,ChargeData(0,0,-1,-oneThird,oneThird)),
    make_pair(-5,ChargeData(0,0,1,oneThird,-oneThird)),
    make_pair(2203,ChargeData(0,0,0,fourThird,twoThird)),
    make_pair(-2203,ChargeData(0,0,0,-fourThird,-twoThird)),
    make_pair(2101,ChargeData(0,0,0,oneThird,twoThird)),
    make_pair(-2101,ChargeData(0,0,0,-oneThird,-twoThird)),
    make_pair(1103,ChargeData(0,0,0,-twoThird,twoThird)),
    make_pair(-1103,ChargeData(0,0,0,twoThird,-twoThird)),
    make_pair(3101,ChargeData(-1,0,0,-twoThird,twoThird)),
    make_pair(-3101,ChargeData(1,0,0,twoThird,-twoThird)),
    make_pair(3201,ChargeData(-1,0,0,oneThird,twoThird)),
    make_pair(-3201,ChargeData(1,0,0,-oneThird,-twoThird)),
    make_pair(3303,ChargeData(-1,0,0,-twoThird,twoThird)),
    make_pair(-3303,ChargeData(1,0,0,twoThird,-twoThird))
};

/**
* Quark/Diquark to charge configurations map (map key being given by quark/diquark id)
*/
static const map<int,ChargeData> quarkChargesMap(quarkChargesPairs,quarkChargesPairs+sizeof(quarkChargesPairs));

#endif
