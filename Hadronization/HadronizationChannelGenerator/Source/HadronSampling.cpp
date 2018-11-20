#include "../Include/HadronSampling.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/Constants.h"
#include <cassert>
#include "TMath.h"

// TODO: Test
#include <iostream>
#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_sf_bessel.h"
struct rparams
{
    vector<double> lightHadronSpinMultiplicities;
    vector<double> lightHadronMasses;
    vector<double> lightHadronElectricCharges;
    vector<double> lightHadronBaryonicCharges;
    vector<double> lightHadronStrangeCharges;
    double clusterSamplingVolume;
    double clusterMass;
    double electricCharge;
    double baryonicCharge;
    double strangeCharge;
};

void print_state4(unsigned int iter, gsl_multiroot_fsolver * s)
{
    printf ("iter = %3u x = % .4f % .4f % .4f % .4f "
            "f(x) = % .3e % .3e % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1),
            gsl_vector_get (s->f, 2),
            gsl_vector_get (s->f, 3));
}

void print_state4(unsigned int iter, gsl_multiroot_fdfsolver * s)
{
    printf ("iter = %3u x = % .4f % .4f % .4f % .4f "
            "f(x) = % .3e % .3e % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1),
            gsl_vector_get (s->f, 2),
            gsl_vector_get (s->f, 3));
}

void print_state4Ex(unsigned int iter, gsl_multiroot_fsolver * s)
{
    printf ("iter = %3u x = % .4f % .4f % .4f % .4f "
            "f(x) = % .3e % .3e % .3e % .3e\n",
            iter,
            std::exp(gsl_vector_get (s->x, 0)),
            std::exp(gsl_vector_get (s->x, 1)),
            std::exp(gsl_vector_get (s->x, 2)),
            gsl_vector_get (s->x, 3),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1),
            gsl_vector_get (s->f, 2),
            gsl_vector_get (s->f, 3));
}

void print_state4Ex(unsigned int iter, gsl_multiroot_fdfsolver * s)
{
    printf ("iter = %3u x = % .4f % .4f % .4f % .4f "
            "f(x) = % .3e % .3e % .3e % .3e\n",
            iter,
            std::exp(gsl_vector_get (s->x, 0)),
            std::exp(gsl_vector_get (s->x, 1)),
            std::exp(gsl_vector_get (s->x, 2)),
            gsl_vector_get (s->x, 3),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1),
            gsl_vector_get (s->f, 2),
            gsl_vector_get (s->f, 3));
}

void print_state3(unsigned int iter, gsl_multiroot_fsolver * s)
{
    printf ("iter = %3u x = % .3f % .3f % .3f "
            "f(x) = % .3e % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1),
            gsl_vector_get (s->f, 2));
}


int computeSaddlePointEquations4(const gsl_vector* x,
                                void* params,
                                gsl_vector* f)
{
    const vector<double> lightHadronSpinMultiplicities(((struct rparams *) params)->lightHadronSpinMultiplicities);
    const vector<double> lightHadronMasses(((struct rparams *) params)->lightHadronMasses);
    const vector<double> lightHadronElectricCharges(((struct rparams *) params)->lightHadronElectricCharges);
    const vector<double> lightHadronBaryonicCharges(((struct rparams *) params)->lightHadronBaryonicCharges);
    const vector<double> lightHadronStrangeCharges(((struct rparams *) params)->lightHadronStrangeCharges);
    
    const double lambdaQ(gsl_vector_get(x,0));
    const double lambdaB(gsl_vector_get(x,1));
    const double lambdaS(gsl_vector_get(x,2));
    const double T(gsl_vector_get(x,3));
    const double samplingVolumeTimesTOn2piSquared(((struct rparams *) params)->clusterSamplingVolume*T/twoPiSquared);
    const double samplingVolumeTimesSquaredTOn2piSquared(samplingVolumeTimesTOn2piSquared*T);
    
    double yElectricCharge(0.);
    double yBaryonicCharge(0.);
    double yStrangeCharge(0.);
    double mass(0.);
    
    const unsigned int numberOfHadrons(lightHadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(lightHadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(lightHadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(lightHadronStrangeCharges[hadronIndex]);
        const double& currentSpinMultiplicity(lightHadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(lightHadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMass(currentMass*currentMass);
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = pow(lambdaQ,currentElectricCharge)*
                       pow(lambdaB,currentBaryonicCharge)*
                       pow(lambdaS,currentStrangeCharge);
        
        yElectricCharge += currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;// TODO: factorize
        yBaryonicCharge += currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        yStrangeCharge += currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
//        mass += squaredMass*currentSpinMultiplicity*(massOnT*gsl_sf_bessel_Kn(3,massOnT)-bessel2)*lambdaFactor;
        mass += squaredMass*currentSpinMultiplicity*(massOnT*gsl_sf_bessel_Kn(1,massOnT)+3.*bessel2)*lambdaFactor;
    }
    
    yElectricCharge = ((struct rparams *) params)->electricCharge -
                      yElectricCharge*samplingVolumeTimesTOn2piSquared;
    yBaryonicCharge = ((struct rparams *) params)->baryonicCharge -
                      yBaryonicCharge*samplingVolumeTimesTOn2piSquared;
    yStrangeCharge = ((struct rparams *) params)->strangeCharge -
                      yStrangeCharge*samplingVolumeTimesTOn2piSquared;
    mass = ((struct rparams *) params)->clusterMass - mass*samplingVolumeTimesSquaredTOn2piSquared;
    
    gsl_vector_set(f,0,yElectricCharge);
    gsl_vector_set(f,1,yBaryonicCharge);
    gsl_vector_set(f,2,yStrangeCharge);
    gsl_vector_set(f,3,mass);
    
    return GSL_SUCCESS;
}

int computeSaddlePointEquations4_df(const gsl_vector* x,
                                 void* params,
                                 gsl_matrix* J)
{
    const vector<double> lightHadronSpinMultiplicities(((struct rparams *) params)->lightHadronSpinMultiplicities);
    const vector<double> lightHadronMasses(((struct rparams *) params)->lightHadronMasses);
    const vector<double> lightHadronElectricCharges(((struct rparams *) params)->lightHadronElectricCharges);
    const vector<double> lightHadronBaryonicCharges(((struct rparams *) params)->lightHadronBaryonicCharges);
    const vector<double> lightHadronStrangeCharges(((struct rparams *) params)->lightHadronStrangeCharges);
    
    const double lambdaQ(gsl_vector_get(x,0));
    const double lambdaB(gsl_vector_get(x,1));
    const double lambdaS(gsl_vector_get(x,2));
    const double T(gsl_vector_get(x,3));
    const double samplingVolumeOn2piSquared(((struct rparams *) params)->clusterSamplingVolume/twoPiSquared);
    const double samplingVolumeTimesTOn2piSquared(samplingVolumeOn2piSquared*T);
    const double samplingVolumeTimesSquaredTOn2piSquared(samplingVolumeTimesTOn2piSquared*T);
    
    
    double dElectricChargeDLambdaQ(0.);
    double dElectricChargeDLambdaB(0.);
    double dElectricChargeDLambdaS(0.);
    double dElectricChargeDT(0.);
    double dBaryonicChargeDLambdaQ(0.);
    double dBaryonicChargeDLambdaB(0.);
    double dBaryonicChargeDLambdaS(0.);
    double dBaryonicChargeDT(0.);
    double dStrangeChargeDLambdaQ(0.);
    double dStrangeChargeDLambdaB(0.);
    double dStrangeChargeDLambdaS(0.);
    double dStrangeChargeDT(0.);
    double dMassDLambdaQ(0.);
    double dMassDLambdaB(0.);
    double dMassDLambdaS(0.);
    double dMassDT(0.);
    
    const unsigned int numberOfHadrons(lightHadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(lightHadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(lightHadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(lightHadronStrangeCharges[hadronIndex]);
        const double& currentSpinMultiplicity(lightHadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(lightHadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMassOnT(massOnT*massOnT);
        const double squaredMass(currentMass*currentMass);
        const double cubeMass(squaredMass*currentMass);
        const double cubeMassOnT(cubeMass/T);
        const double bessel1(gsl_sf_bessel_Kn(1,massOnT));
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = pow(lambdaQ,currentElectricCharge)*
                       pow(lambdaB,currentBaryonicCharge)*
                       pow(lambdaS,currentStrangeCharge);
        
        
        
        dElectricChargeDLambdaQ +=
            currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge/lambdaQ;// TODO: factorize
        dElectricChargeDLambdaB +=
            currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge/lambdaB;
        dElectricChargeDLambdaS +=
            currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge/lambdaS;
        dElectricChargeDT +=
            currentElectricCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        dBaryonicChargeDLambdaQ +=
            currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge/lambdaQ;
        dBaryonicChargeDLambdaB +=
            currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge/lambdaB;
        dBaryonicChargeDLambdaS +=
            currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge/lambdaS;
        dBaryonicChargeDT +=
            currentBaryonicCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        dStrangeChargeDLambdaQ +=
            currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge/lambdaQ;
        dStrangeChargeDLambdaB +=
            currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge/lambdaB;
        dStrangeChargeDLambdaS +=
            currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge/lambdaS;
        dStrangeChargeDT +=
            currentStrangeCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        dMassDLambdaQ +=
            squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentElectricCharge/lambdaQ;
        dMassDLambdaB +=
            squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentBaryonicCharge/lambdaB;
        dMassDLambdaS +=
            squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentStrangeCharge/lambdaS;
        dMassDT +=
            squaredMass*currentSpinMultiplicity*(2.*massOnT*bessel1+(12.+squaredMassOnT)*bessel2)*lambdaFactor;
    }

    
    dElectricChargeDLambdaQ = -dElectricChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaB = -dElectricChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaS = -dElectricChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDT = -dElectricChargeDT*samplingVolumeOn2piSquared;
    dBaryonicChargeDLambdaQ = -dBaryonicChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaB = -dBaryonicChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaS = -dBaryonicChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDT = -dBaryonicChargeDT*samplingVolumeOn2piSquared;
    dStrangeChargeDLambdaQ = -dStrangeChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaB = -dStrangeChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaS = -dStrangeChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDT = -dStrangeChargeDT*samplingVolumeOn2piSquared;
    dMassDLambdaQ = -dMassDLambdaQ*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaB = -dMassDLambdaB*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaS = -dMassDLambdaS*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDT = -dMassDT*samplingVolumeTimesTOn2piSquared;
    
    
    gsl_matrix_set (J, 0, 0, dElectricChargeDLambdaQ);
    gsl_matrix_set (J, 0, 1, dElectricChargeDLambdaB);
    gsl_matrix_set (J, 0, 2, dElectricChargeDLambdaS);
    gsl_matrix_set (J, 0, 3, dElectricChargeDT);
    gsl_matrix_set (J, 1, 0, dBaryonicChargeDLambdaQ);
    gsl_matrix_set (J, 1, 1, dBaryonicChargeDLambdaB);
    gsl_matrix_set (J, 1, 2, dBaryonicChargeDLambdaS);
    gsl_matrix_set (J, 1, 3, dBaryonicChargeDT);
    gsl_matrix_set (J, 2, 0, dStrangeChargeDLambdaQ);
    gsl_matrix_set (J, 2, 1, dStrangeChargeDLambdaB);
    gsl_matrix_set (J, 2, 2, dStrangeChargeDLambdaS);
    gsl_matrix_set (J, 2, 3, dStrangeChargeDT);
    gsl_matrix_set (J, 3, 0, dMassDLambdaQ);
    gsl_matrix_set (J, 3, 1, dMassDLambdaB);
    gsl_matrix_set (J, 3, 2, dMassDLambdaS);
    gsl_matrix_set (J, 3, 3, dMassDT);
    
    return GSL_SUCCESS;
}

int computeSaddlePointEquations4_fdf(const gsl_vector* x,
                                    void* params,
                                     gsl_vector* f,
                                    gsl_matrix* J)
{
    computeSaddlePointEquations4(x,params,f);
    computeSaddlePointEquations4_df(x,params,J);
    return GSL_SUCCESS;
}

int computeSaddlePointEquations4Ex(const gsl_vector* x,
                                 void* params,
                                 gsl_vector* f)
{
    const vector<double> lightHadronSpinMultiplicities(((struct rparams *) params)->lightHadronSpinMultiplicities);
    const vector<double> lightHadronMasses(((struct rparams *) params)->lightHadronMasses);
    const vector<double> lightHadronElectricCharges(((struct rparams *) params)->lightHadronElectricCharges);
    const vector<double> lightHadronBaryonicCharges(((struct rparams *) params)->lightHadronBaryonicCharges);
    const vector<double> lightHadronStrangeCharges(((struct rparams *) params)->lightHadronStrangeCharges);
    
    const double muQ(gsl_vector_get(x,0));
    const double muB(gsl_vector_get(x,1));
    const double muS(gsl_vector_get(x,2));
    const double T(gsl_vector_get(x,3));
    const double samplingVolumeTimesTOn2piSquared(((struct rparams *) params)->clusterSamplingVolume*T/twoPiSquared);
    const double samplingVolumeTimesSquaredTOn2piSquared(samplingVolumeTimesTOn2piSquared*T);
    
    double yElectricCharge(0.);
    double yBaryonicCharge(0.);
    double yStrangeCharge(0.);
    double mass(0.);
    
    const unsigned int numberOfHadrons(lightHadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(lightHadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(lightHadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(lightHadronStrangeCharges[hadronIndex]);
        const double& currentSpinMultiplicity(lightHadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(lightHadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMass(currentMass*currentMass);
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = std::exp(muQ*currentElectricCharge)*
        std::exp(muB*currentBaryonicCharge)*
        std::exp(muS*currentStrangeCharge);
        
        yElectricCharge += currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;// TODO: factorize
        yBaryonicCharge += currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        yStrangeCharge += currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        //        mass += squaredMass*currentSpinMultiplicity*(massOnT*gsl_sf_bessel_Kn(3,massOnT)-bessel2)*lambdaFactor;
        mass += squaredMass*currentSpinMultiplicity*(massOnT*gsl_sf_bessel_Kn(1,massOnT)+3.*bessel2)*lambdaFactor;
    }
    
    yElectricCharge = ((struct rparams *) params)->electricCharge -
    yElectricCharge*samplingVolumeTimesTOn2piSquared;
    yBaryonicCharge = ((struct rparams *) params)->baryonicCharge -
    yBaryonicCharge*samplingVolumeTimesTOn2piSquared;
    yStrangeCharge = ((struct rparams *) params)->strangeCharge -
    yStrangeCharge*samplingVolumeTimesTOn2piSquared;
    mass = ((struct rparams *) params)->clusterMass - mass*samplingVolumeTimesSquaredTOn2piSquared;
    
    gsl_vector_set(f,0,yElectricCharge);
    gsl_vector_set(f,1,yBaryonicCharge);
    gsl_vector_set(f,2,yStrangeCharge);
    gsl_vector_set(f,3,mass);
    
    return GSL_SUCCESS;
}

int computeSaddlePointEquations4Ex_df(const gsl_vector* x,
                                    void* params,
                                    gsl_matrix* J)
{
    const vector<double> lightHadronSpinMultiplicities(((struct rparams *) params)->lightHadronSpinMultiplicities);
    const vector<double> lightHadronMasses(((struct rparams *) params)->lightHadronMasses);
    const vector<double> lightHadronElectricCharges(((struct rparams *) params)->lightHadronElectricCharges);
    const vector<double> lightHadronBaryonicCharges(((struct rparams *) params)->lightHadronBaryonicCharges);
    const vector<double> lightHadronStrangeCharges(((struct rparams *) params)->lightHadronStrangeCharges);
    
    const double muQ(gsl_vector_get(x,0));
    const double muB(gsl_vector_get(x,1));
    const double muS(gsl_vector_get(x,2));
    const double T(gsl_vector_get(x,3));
    const double samplingVolumeOn2piSquared(((struct rparams *) params)->clusterSamplingVolume/twoPiSquared);
    const double samplingVolumeTimesTOn2piSquared(samplingVolumeOn2piSquared*T);
    const double samplingVolumeTimesSquaredTOn2piSquared(samplingVolumeTimesTOn2piSquared*T);
    
    
    double dElectricChargeDLambdaQ(0.);
    double dElectricChargeDLambdaB(0.);
    double dElectricChargeDLambdaS(0.);
    double dElectricChargeDT(0.);
    double dBaryonicChargeDLambdaQ(0.);
    double dBaryonicChargeDLambdaB(0.);
    double dBaryonicChargeDLambdaS(0.);
    double dBaryonicChargeDT(0.);
    double dStrangeChargeDLambdaQ(0.);
    double dStrangeChargeDLambdaB(0.);
    double dStrangeChargeDLambdaS(0.);
    double dStrangeChargeDT(0.);
    double dMassDLambdaQ(0.);
    double dMassDLambdaB(0.);
    double dMassDLambdaS(0.);
    double dMassDT(0.);
    
    const unsigned int numberOfHadrons(lightHadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(lightHadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(lightHadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(lightHadronStrangeCharges[hadronIndex]);
        const double& currentSpinMultiplicity(lightHadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(lightHadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMassOnT(massOnT*massOnT);
        const double squaredMass(currentMass*currentMass);
        const double cubeMass(squaredMass*currentMass);
        const double cubeMassOnT(cubeMass/T);
        const double bessel1(gsl_sf_bessel_Kn(1,massOnT));
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = std::exp(muQ*currentElectricCharge)*
        std::exp(muB*currentBaryonicCharge)*
        std::exp(muS*currentStrangeCharge);
        
        
        
        dElectricChargeDLambdaQ +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;// TODO: factorize
        dElectricChargeDLambdaB +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dElectricChargeDLambdaS +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dElectricChargeDT +=
        currentElectricCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        dBaryonicChargeDLambdaQ +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;
        dBaryonicChargeDLambdaB +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dBaryonicChargeDLambdaS +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dBaryonicChargeDT +=
        currentBaryonicCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        dStrangeChargeDLambdaQ +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;
        dStrangeChargeDLambdaB +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dStrangeChargeDLambdaS +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dStrangeChargeDT +=
        currentStrangeCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        dMassDLambdaQ +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentElectricCharge;
        dMassDLambdaB +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentBaryonicCharge;
        dMassDLambdaS +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentStrangeCharge;
        dMassDT +=
        squaredMass*currentSpinMultiplicity*(2.*massOnT*bessel1+(12.+squaredMassOnT)*bessel2)*lambdaFactor;
    }
    
    
    dElectricChargeDLambdaQ = -dElectricChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaB = -dElectricChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaS = -dElectricChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDT = -dElectricChargeDT*samplingVolumeOn2piSquared;
    dBaryonicChargeDLambdaQ = -dBaryonicChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaB = -dBaryonicChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaS = -dBaryonicChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDT = -dBaryonicChargeDT*samplingVolumeOn2piSquared;
    dStrangeChargeDLambdaQ = -dStrangeChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaB = -dStrangeChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaS = -dStrangeChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDT = -dStrangeChargeDT*samplingVolumeOn2piSquared;
    dMassDLambdaQ = -dMassDLambdaQ*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaB = -dMassDLambdaB*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaS = -dMassDLambdaS*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDT = -dMassDT*samplingVolumeTimesTOn2piSquared;
    
    
    gsl_matrix_set (J, 0, 0, dElectricChargeDLambdaQ);
    gsl_matrix_set (J, 0, 1, dElectricChargeDLambdaB);
    gsl_matrix_set (J, 0, 2, dElectricChargeDLambdaS);
    gsl_matrix_set (J, 0, 3, dElectricChargeDT);
    gsl_matrix_set (J, 1, 0, dBaryonicChargeDLambdaQ);
    gsl_matrix_set (J, 1, 1, dBaryonicChargeDLambdaB);
    gsl_matrix_set (J, 1, 2, dBaryonicChargeDLambdaS);
    gsl_matrix_set (J, 1, 3, dBaryonicChargeDT);
    gsl_matrix_set (J, 2, 0, dStrangeChargeDLambdaQ);
    gsl_matrix_set (J, 2, 1, dStrangeChargeDLambdaB);
    gsl_matrix_set (J, 2, 2, dStrangeChargeDLambdaS);
    gsl_matrix_set (J, 2, 3, dStrangeChargeDT);
    gsl_matrix_set (J, 3, 0, dMassDLambdaQ);
    gsl_matrix_set (J, 3, 1, dMassDLambdaB);
    gsl_matrix_set (J, 3, 2, dMassDLambdaS);
    gsl_matrix_set (J, 3, 3, dMassDT);
    
    return GSL_SUCCESS;
}

int computeSaddlePointEquations4Ex_fdf(const gsl_vector* x,
                                     void* params,
                                     gsl_vector* f,
                                     gsl_matrix* J)
{
    computeSaddlePointEquations4Ex(x,params,f);
    computeSaddlePointEquations4Ex_df(x,params,J);
    return GSL_SUCCESS;
}


// TODO: Test

int computeSaddlePointEquations3(const gsl_vector* x,
                                 void* params,
                                 gsl_vector* f)
{
    const vector<double> lightHadronSpinMultiplicities(((struct rparams *) params)->lightHadronSpinMultiplicities);
    const vector<double> lightHadronMasses(((struct rparams *) params)->lightHadronMasses);
    const vector<double> lightHadronElectricCharges(((struct rparams *) params)->lightHadronElectricCharges);
    const vector<double> lightHadronBaryonicCharges(((struct rparams *) params)->lightHadronBaryonicCharges);
    const vector<double> lightHadronStrangeCharges(((struct rparams *) params)->lightHadronStrangeCharges);
    
    const double lambdaQ(gsl_vector_get(x,0));
    const double lambdaB(gsl_vector_get(x,1));
    const double lambdaS(gsl_vector_get(x,2));
    //    const double T(gsl_vector_get(x,3));
    const double T(0.160);
    const double squaredT(T*T);
    const double samplingVolumeTimesSquaredTOn2piSquared(((struct rparams *) params)->clusterSamplingVolume*squaredT/twoPiSquared);
    
    double yElectricCharge(0.);
    double yBaryonicCharge(0.);
    double yStrangeCharge(0.);
    double mass(0.);
    
    const unsigned int numberOfHadrons(lightHadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(lightHadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(lightHadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(lightHadronStrangeCharges[hadronIndex]);
        const double& currentSpinMultiplicity(lightHadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(lightHadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMass(currentMass*currentMass);
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = pow(lambdaQ,currentElectricCharge)*
        pow(lambdaB,currentBaryonicCharge)*
        pow(lambdaS,currentStrangeCharge);
        
        yElectricCharge += currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;// TODO: factorize
        yBaryonicCharge += currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        yStrangeCharge += currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        mass += squaredMass*currentSpinMultiplicity*(massOnT*gsl_sf_bessel_Kn(3,massOnT)-bessel2)*lambdaFactor;
    }
    
    yElectricCharge = ((struct rparams *) params)->electricCharge -
    yElectricCharge*samplingVolumeTimesSquaredTOn2piSquared;
    yBaryonicCharge = ((struct rparams *) params)->baryonicCharge -
    yBaryonicCharge*samplingVolumeTimesSquaredTOn2piSquared;
    yStrangeCharge = ((struct rparams *) params)->strangeCharge -
    yStrangeCharge*samplingVolumeTimesSquaredTOn2piSquared;
    mass = ((struct rparams *) params)->clusterMass - mass*samplingVolumeTimesSquaredTOn2piSquared;
    
    gsl_vector_set(f,0,yElectricCharge);
    gsl_vector_set(f,1,yBaryonicCharge);
    gsl_vector_set(f,2,yStrangeCharge);
    
    return GSL_SUCCESS;
}
// TODO: Test



HadronSampling::HadronSampling(const HadronSamplingGroups& i_hadronSamplingGroups,
                               const double i_samplingTemperature,
                               const double i_samplingEnergyDensity,
                               RandomNumberGenerator& io_randomGenerator)
                              :m_isHadronSamplingReady(false)
                              ,m_isHadronSetAvailable(false)
                              ,m_weight(0.)
                              ,m_hadronSamplingGroups(&i_hadronSamplingGroups)
                              ,m_randomGenerator(&io_randomGenerator)
                              ,m_samplingTemperature(i_samplingTemperature)
                              ,m_samplingEnergyDensity(i_samplingEnergyDensity/gevCubeToFermiCube)
                              ,m_clusterSamplingVolume(0.)
                              ,m_clusterMass(0.)
                              ,m_residualStrangeCharge(0)
                              ,m_residualElectricCharge(0.)
                              ,m_residualBaryonicCharge(0.)
                              ,m_lightHadronExponentialWeightFactor(0.)
                              ,m_acceptedHeavyHadronNumber(0)
                              ,m_hadronMassSum(0.)
                              ,m_runNewCalculation(false)// TEST
                              ,m_lambdaElectricCharge(1.)// TEST
                              ,m_lambdaStrangeCharge(1.)// TEST
                              ,m_lambdaBaryonicCharge(1.)
                              ,m_samplingT(0.160)// TEST
{
    try
    {
        if(m_samplingTemperature<=0.)
        {
            throw HadronizationException("Error during hadron sampling object creation, non positive sampling temperature provided",
                                         __FUNCTION__,111);
        }
        if(m_samplingEnergyDensity<=0.)
        {
            throw HadronizationException("Error during hadron sampling object creation, non positive sampling energy density provided",
                                         __FUNCTION__,112);            
        }
        if(m_randomGenerator==NULL)
        {
            throw HadronizationException("Error during phase space sampling object creation, NULL pointer to random number generator provided",
                                         __FUNCTION__,113);
        }
        computeSamplingData();
        m_isHadronSamplingReady = true;
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

HadronSampling::~HadronSampling(void)
{
}

unsigned int HadronSampling::run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
                                 const double i_mass,
                                 const double i_energyDensity)
{
    if(!m_isHadronSamplingReady)
    {
        throw HadronizationException("Error during hadron sampling, sampling class not correctly initialized",
                                     __FUNCTION__,114);
    }
    
    m_isHadronSetAvailable = false;
    
    try
    {
        
        // Store cluster data
        m_clusterMass = i_mass;
        const int clusterStrangeCharge(i_chargeConfiguration.strangeCharge);
        const double clusterBaryonicCharge(i_chargeConfiguration.baryonicCharge);
        const double clusterElectricCharge(i_chargeConfiguration.electricCharge);

        // Compute cluster mass less phase space availability tolerance value (for channel acceptance test only)
        // Checking the decay channel total mass wrt this reduced cluster mass value, instead of the
        // real one, allows to include in the simulation and partition function calculation only channels with
        // a reasonable phase space integral value.
        const double clusterMassLessTolerance(m_clusterMass-phaseSpaceEnergyAvailabilityTolerance);
        if(clusterMassLessTolerance<=0.)
        {
            throw HadronizationException("Error during hadron sampling, non positive cluster mass less phase space tolerance",
                                         __FUNCTION__,117);
        }
        
        // Set cluster sampling volume parameter
        m_clusterSamplingVolume = m_clusterMass/m_samplingEnergyDensity;
        
        
        // If cluster is heavy flavored set the corresponding heavy flavored hadron sampling probabilities
        const bool isClusterHeavyFlavored(((i_chargeConfiguration.charmCharge!=0) ||
                                          (i_chargeConfiguration.bottomCharge!=0)));
        if(isClusterHeavyFlavored)
        {            
            const int clusterCharmCharge(i_chargeConfiguration.charmCharge);
            const int clusterBottomCharge(i_chargeConfiguration.bottomCharge);

            if(clusterCharmCharge!=0 && clusterBottomCharge!=0)
            {
                // Double heavy flavored cluster, throw exception
                return 3; // TODO: switch to exit code static table/enumerator (also for exception)
            }
            
            // Set heavy flavored sampling probability (and hadron group) for the provided cluster composition
            const bool samplingProbabilitySetStatus = setHeavyHadronSamplingProbabilities(clusterCharmCharge,clusterBottomCharge);
            if(samplingProbabilitySetStatus)
            {
                return 2;
            }
        }
        else
        {
            // TODO: TEST
            if(m_runNewCalculation)
            {
                m_residualMass = m_clusterMass;
                computeFugacities(i_chargeConfiguration.electricCharge,
                                  i_chargeConfiguration.strangeCharge,
                                  i_chargeConfiguration.baryonicCharge,
                                  m_lambdaElectricCharge,
                                  m_lambdaStrangeCharge,
                                  m_lambdaBaryonicCharge,
                                  m_samplingT);
                m_runNewCalculation = false;
            }
            // TODO: TEST
        }
        
        
        // TODO: check channel composition (e.g., does neutral channel contain ssbar particles??)
        // TODO: set parameter for sampling attempt count limit
        unsigned long long int samplingAttempt=0;// TODO: long long or long?
        unsigned long long int maxSamplingAttempt;
        if(isClusterHeavyFlavored)
        {
            // TODO: fix this
//            maxSamplingAttempt = 1000000000;
            maxSamplingAttempt = 100000;
        }
        else
        {
            // TODO: fix this
//            maxSamplingAttempt = 1000000000;
            maxSamplingAttempt = 100000;
        }
        
        do
        {
            m_weight = 1.;
            m_hadronMassSum = 0.;
            ++samplingAttempt;// TODO: see comment above concerning sampling attempt handling (generation)
                              // TODO: no sampling limit is foreseen for part. function calculation, check class usage case! 
            
            // TODO: there is no reason to remove the map keys I think.. Close the weight calculation
            //       open question and then optimize this (maybe old hadron create problem in new event generation)
            // TODO: make similar check for all other containers used
            m_hadronSet.clear();
            m_hadronSetPositionMap.clear();
                        
            // Store cluster charges
            m_residualStrangeCharge = clusterStrangeCharge;
            m_residualElectricCharge = clusterElectricCharge;
            m_residualBaryonicCharge = clusterBaryonicCharge;
                        
            if(isClusterHeavyFlavored)
            {
                // Heavy flavored cluster: run heavy flavored hadron sampling first
                runHadronSampling(m_heavyHadronSamplingGroup,true);
                
                computeFugacities(m_residualElectricCharge,
                                  m_residualStrangeCharge,
                                  m_residualBaryonicCharge,
                                  m_lambdaElectricCharge,
                                  m_lambdaStrangeCharge,
                                  m_lambdaBaryonicCharge,
                                  m_samplingT);
            }
            
            // Run light flavored hadron sampling
            runLightHadronSampling();
            
            // Check validity of the sampled hadronization channel
            if(m_isHadronSetAvailable)
            {
                // Check number of sampled hadrons and their mass sum
                // TODO: add double comparison tolerance
                if((m_hadronSet.size()<2) || (clusterMassLessTolerance<=m_hadronMassSum))
                {
                    m_isHadronSetAvailable = false;
                }
            }
            else// TODO: debug, for comparison with old code only, to be removed!
            {
//                if(!isClusterHeavyFlavored)
//                {
//                    samplingAttempt -= 1;// TODO: no sampling limit is foreseen for part. function calculation, check class usage case!
//                }
            }
            
        }while(m_isHadronSetAvailable==false && samplingAttempt<maxSamplingAttempt);
        
        // TODO: see comment above concerning sampling attempt handling
        if(m_isHadronSetAvailable==false)
        {
            // Maximum number of sampling attempt exceeded
            return 1;// TODO: use static table??
        }
        
        // Compute sampling weight
        computeSampledChannelWeight(i_energyDensity);
        
        // Sampling correctly performed
        return 0;// TODO: use static table??
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
}

// TODO: debug, after reordering removal this method must be const
const vector<HadronData>& HadronSampling::getHadronSet(void)
{
    if(m_isHadronSetAvailable)
    {
        // TODO: debug, remove after comparison with old code
        if(m_hadronSet.begin()->getCharmCharge() || m_hadronSet.begin()->getBottomCharge())
        {
            vector<HadronData> orderedHadronSet(m_hadronSet);
            orderedHadronSet.erase(orderedHadronSet.begin());
            orderedHadronSet.push_back(*m_hadronSet.begin());
            m_hadronSet = orderedHadronSet;
        }
        
        return m_hadronSet;
    }
    else
    {
        throw HadronizationException("Error during sampled hadron set retrieval, hadrons not available",
                                     __FUNCTION__,115);
    }
}

unsigned int HadronSampling::getNumberOfHadrons(void) const
{
    if(m_isHadronSetAvailable)
    {
        return m_hadronSet.size();
    }
    else
    {
        throw HadronizationException("Error during hadron sampling weight retrieval, hadrons not available",
                                     __FUNCTION__,115);
    }
}

double HadronSampling::getSamplingWeight(void) const
{
    if(m_isHadronSetAvailable)
    {
        return m_weight;
    }
    else
    {
        throw HadronizationException("Error during hadron sampling weight retrieval, hadrons not available",
                                     __FUNCTION__,115);
    }
}

void HadronSampling::computeSamplingData(void)
{
    
    try
    {
        // Retrieve hadron sampling group list
        const vector<HadronGroup> groupList(m_hadronSamplingGroups->getHadronGroupList());
       
        // Reset m_lightHadronExponentialWeightFactor (updated in computeHadronSamplingWeights)
        m_lightHadronExponentialWeightFactor = 0.;
        
        // Loop over sampling group
        for(vector<HadronGroup>::const_iterator hadronGroup = groupList.begin();
            hadronGroup != groupList.end();++hadronGroup)
        {
            
            // Compute hadron sampling weights
            switch(*hadronGroup)
            {                    
                // Light hadron groups
                case LightBaryons:
                case LightAntibaryons:
                case LightStrangeMesons:
                case LightStrangeAntimesons:
                case LightChargedMesons:
                case LightChargedAntimesons:
                case LightNeutralMesons:
                    computeHadronSamplingWeights(*hadronGroup,true);
                    break;
                    
                // Heavy flavored hadron groups
                case CharmedHadrons:
                case AntiCharmedHadrons:
                case BottomedHadrons:
                case AntiBottomedHadrons:
                    computeHadronSamplingWeights(*hadronGroup,false);
                    break;

                default:
                    throw HadronizationException("Error during hadron sampling probability calculation, hadron group not identified",
                                                 __FUNCTION__,116);
            }
            
        }
        
        // Close light hadron exponential sampling factor calculation
        m_lightHadronExponentialWeightFactor = std::exp(m_lightHadronExponentialWeightFactor);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void HadronSampling::computeHadronSamplingWeights(const HadronGroup i_hadronGroup,const bool i_isLightGroup)
{
    try
    {
        
        const vector<const HadronData*> hadronData(m_hadronSamplingGroups->getHadronGroup(i_hadronGroup));
        const unsigned int numberOfHadrons(hadronData.size());
        double meanMultiplicity = 0.;
        
        // TODO: remove for loop duplication (Functor or pointer to function)
        if(i_isLightGroup)
        {
            // TEST
            const double lambda_q(1.0);
            const double lambda_s(1.0);
            const double lambda_b(1.0);
//            const double lambda_q(2.96032390071886464);
//            const double lambda_s(0.32777547018731618);
//            const double lambda_b(0.02845072155233162);
            // TEST
            
            double meanMultiplicitySum = 0.;
            unsigned int hadronSpinMultiplicity = 0;
            
            // Loop over hadrons and compute mean multiplicities
            vector<double> meanMultiplicities(numberOfHadrons,0.);
            vector<double> samplingWeightFactors(numberOfHadrons,0.);
            for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons;++hadronIndex)
            {
                hadronSpinMultiplicity = hadronData[hadronIndex]->getSpinMultiplicity();
                // This weight is the grandcanonical mean multiplicity value except for the multiplicative cluster volume factor
                meanMultiplicity = computeLightHadronMeanMultiplicity(hadronData[hadronIndex]->getMass(),hadronSpinMultiplicity);
                
                
                meanMultiplicity *= pow(lambda_q,hadronData[hadronIndex]->getElectricCharge())*
                                    pow(lambda_s,hadronData[hadronIndex]->getStrangeCharge())*
                                    pow(lambda_b,hadronData[hadronIndex]->getBaryonicCharge());
                
                
                if(meanMultiplicity>0.)
                {
                    meanMultiplicities[hadronIndex] = meanMultiplicity;
                    meanMultiplicitySum += meanMultiplicity;
                    m_lightHadronExponentialWeightFactor += meanMultiplicity;
                    samplingWeightFactors[hadronIndex] = hadronSpinMultiplicity/meanMultiplicity;
                    
                    // TODO: test
//                    if(hadronData[hadronIndex]->getMass()<m_clusterMass)
//                    {
                        m_lightHadronSpinMultiplicities.push_back(hadronData[hadronIndex]->getSpinMultiplicity());
                        m_lightHadronMasses.push_back(hadronData[hadronIndex]->getMass());
                        m_lightHadronElectricCharges.push_back(hadronData[hadronIndex]->getElectricCharge());
                        m_lightHadronBaryonicCharges.push_back(hadronData[hadronIndex]->getBaryonicCharge());
                        m_lightHadronStrangeCharges.push_back(static_cast<double>(hadronData[hadronIndex]->getStrangeCharge()));
//                    }
                    // TODO: test
                    
                }
                else
                {
                    throw HadronizationException("Error during hadron sampling probability calculation, light hadron with null mean multiplicity found",
                                                 __FUNCTION__,116);
                }
            }

            // TODO: all hadron sampling weights are very small, is there the possibility to
            // identify a common factor causing the weight range and to factorize it??
        
            // This is the sampling group grandcanonical mean multiplicity value except for
            // the multiplicative cluster volume factor
            m_hadronGroupMeanMultiplicities[i_hadronGroup] = meanMultiplicitySum;

            // Set light flavored hadron sampling weight factor (see .h for details)
            m_lightHadronMultiplicityWeightFactors[i_hadronGroup] = samplingWeightFactors;
            
            // Loop over hadrons and compute cumulative sampling probabilities
            double cumulativeMean = 0.;
            vector<double> cumulativeProbabilities(numberOfHadrons,0.);
            for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons-1;++hadronIndex)
            {
                // Compute cumulative probability
                cumulativeMean += meanMultiplicities[hadronIndex]/meanMultiplicitySum;
                
                // Set cumulative probability
                cumulativeProbabilities[hadronIndex] = cumulativeMean;
            }
            cumulativeProbabilities[numberOfHadrons-1] = 1.;
            
            // Set light flavored hadron cumulative probabilities
            m_samplingCumulativeProbabilities[i_hadronGroup] = cumulativeProbabilities;

        }
        else
        {
                        
            double hadronMass;
            vector<double> hadronMasses(numberOfHadrons,0.);
            vector<double> samplingWeights(numberOfHadrons,0.);
            // Loop over hadrons and compute sampling weights (normalization is included setHeavyHadronSamplingProbabilities)
            // Heavy hadron masses are also stored for faster sampling probability calculation (in setHeavyHadronSamplingProbabilities)
            for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons;++hadronIndex)
            {
                hadronMass = hadronData[hadronIndex]->getMass();
                hadronMasses[hadronIndex] = hadronMass;
                meanMultiplicity = computeHeavyHadronSamplingWeight(hadronMass);
                if(meanMultiplicity>0.)
                {
                    samplingWeights[hadronIndex] = meanMultiplicity;
                }
                else
                {
                    throw HadronizationException("Error during hadron sampling probability calculation, heavy hadron with null mean multiplicity found",
                                                 __FUNCTION__,116);                
                }
            }
            m_heavyHadronMasses[i_hadronGroup] = hadronMasses;

            // Set heavy flavored hadron sampling weights
            m_heavyHadronSamplingWeights[i_hadronGroup] = samplingWeights;
        }
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
}

void HadronSampling::runLightHadronSampling(void)
{
    // TODO: check the possibility to perform a single random number generation call
    
    // TODO: mean multiplicities, even including sampling volume, are very small...why? Can we refine this???
    //       a lot of sampling attempts lead to zero hadrons at present!!!
    
    // TODO: can be the mutiplication per m_clusterSamplingVolume optimized
    // (in case of light flavored cluster m_clusterSamplingVolume is fixed during sampling attempts)
    
    try
    {
        // Generate the number of baryons and antibaryons
        int numberOfHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightBaryons]);
        int numberOfAntiHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightAntibaryons]);
        // Check baryonic charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualBaryonicCharge)
        {
            // Baryonic charge violation, re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(LightBaryons,false,numberOfHadron);
        runHadronSampling(LightAntibaryons,false,numberOfAntiHadron);
        
        // Generate the number of strange mesons and antimesons
        numberOfHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightStrangeMesons]);
        numberOfAntiHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightStrangeAntimesons]);
        // Check strange charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualStrangeCharge)
        {
            // Strange charge violation, re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(LightStrangeMesons,false,numberOfHadron);
        runHadronSampling(LightStrangeAntimesons,false,numberOfAntiHadron);
        
        // Generate the number of electrically charged non strange mesons and antimesons
        numberOfHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightChargedMesons]);
        numberOfAntiHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightChargedAntimesons]);
        // Check electric charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualElectricCharge)
        {
            // Strange charge violation re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(LightChargedMesons,false,numberOfHadron);
        runHadronSampling(LightChargedAntimesons,false,numberOfAntiHadron);
        
        // Generate the number of neutral mesons
        numberOfHadron = m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightNeutralMesons]);// TODO: avoid local
        // Run hadron/antihadron sampling
        runHadronSampling(LightNeutralMesons,false,numberOfHadron);
                
        // Hadron sampling completed
        m_isHadronSetAvailable = true;
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void HadronSampling::runHadronSampling(const HadronGroup i_hadronGroup,const bool i_isHeavyGroup,const unsigned int i_numberOfHadrons)
{
    try
    {
        // Check hadron availability
        assert(m_samplingCumulativeProbabilities[i_hadronGroup].size()>0);
        m_hadronSetPositionMap[i_hadronGroup].resize(i_numberOfHadrons);
        // Retrieve number of hadron minus one of the provided group
        unsigned int numberOfHadronsMinusOne;
        if(i_isHeavyGroup)
        {
            numberOfHadronsMinusOne = m_acceptedHeavyHadronNumber - 1;
        }
        else
        {
            numberOfHadronsMinusOne = m_samplingCumulativeProbabilities[i_hadronGroup].size() - 1;
        }
        
        // Perform i_numberOfHadrons hadron samplings
        for(unsigned int hadronCount = 0;hadronCount<i_numberOfHadrons;++hadronCount)
        {
            // Generate random number from uniform distribution
            const double randomUniform(m_randomGenerator->getUniformRandom());
            
            // Select the hadron using the stored cumulative probabilities
            unsigned int hadronIndex = 0;
            vector<double>::const_iterator hadronIt = m_samplingCumulativeProbabilities[i_hadronGroup].begin();
            // Cumulative function value comparison with the sampled random number is not required for
            // the last hadron. For this reason numberOfHadronsMinusOne is used.
            for(;hadronIndex<numberOfHadronsMinusOne;++hadronIndex)
            {
                if(randomUniform <= *hadronIt)
                {
                    break;
                }
                ++hadronIt;
            }
        
            // TODO: is this check correct? What about closed heavy flavored hadrons??? Same question for other similar checks of this class!
            if(i_isHeavyGroup)
            {
                // Heavy hadron sampling weight update
                m_weight /= m_heavyHadronSamplingProbabilities[hadronIndex];
                // Retrieve sampled hadron data
                hadronIndex = m_availableHeavyHadronPositions[hadronIndex];
            }
            
            // Retrieve sampled hadron data
            const HadronData& sampledHadronData(*((m_hadronSamplingGroups->getHadronGroup(i_hadronGroup))[hadronIndex]));

            // TODO: is this check correct? What about closed heavy flavored hadrons??? Same question for other similar checks of this class!
            if(i_isHeavyGroup)
            {
                // With this sampling weight update the heavy flavor constibution is complete 
                m_weight *= sampledHadronData.getSpinMultiplicity();
                
                // Update sampling volume fol light hadron sampling
                m_clusterSamplingVolume = (m_clusterMass-sampledHadronData.getMass())/m_samplingEnergyDensity;
                
                // TODO: TEST
                m_residualMass = m_clusterMass-sampledHadronData.getMass();
                // TODO: TEST
            }
            
            // Update cluster residual data
            m_residualStrangeCharge -= sampledHadronData.getStrangeCharge();
            m_residualElectricCharge -= sampledHadronData.getElectricCharge();
            m_residualBaryonicCharge -= sampledHadronData.getBaryonicCharge();
            
            // Update sampled hadron set (sampling probabilities are ordered according to hadron sampling group vector)
            // TODO: can we make a guess about vector lenght???? (use resize with i_numberOfHadrons)
            m_hadronSet.push_back(sampledHadronData);
            
            // TODO: this is needed for light hadrons only
            m_hadronSetPositionMap[i_hadronGroup][hadronCount] = hadronIndex;
                        
            // Update hadron mass sum
            m_hadronMassSum += sampledHadronData.getMass();
        }
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
    
}

bool HadronSampling::setHeavyHadronSamplingProbabilities(const int i_clusterCharmCharge,const int i_clusterBottomCharge)
{
    // TODO: after all checks it would be nice to have mass ordered heavy hadrons
    // for faster selection
    
    // Identify which heavy flavored hadron group sampling weights must be normalized
    if(i_clusterCharmCharge>0)
    {
        m_heavyHadronSamplingGroup = CharmedHadrons;
    }
    else if(i_clusterCharmCharge<0)
    {
        m_heavyHadronSamplingGroup = AntiCharmedHadrons;
    }
    else if(i_clusterBottomCharge>0)
    {
        m_heavyHadronSamplingGroup = AntiBottomedHadrons;
    }
    else
    {
        m_heavyHadronSamplingGroup = BottomedHadrons;
    }

    // Access pointers to heavy hadron sampling weights and masses
    const vector<double>* hadronSamplingWeights(&m_heavyHadronSamplingWeights[m_heavyHadronSamplingGroup]);
    const vector<double>* hadronMasses(&m_heavyHadronMasses[m_heavyHadronSamplingGroup]);
    
    double weight = 0.;
    double weightSum = 0.;
    const unsigned int numberOfHadrons(hadronMasses->size());
    m_heavyHadronSamplingProbabilities.resize(numberOfHadrons,0.);
    m_availableHeavyHadronPositions.resize(numberOfHadrons,0.);
    // Compute normalization for hadron accepted for sampling
    m_acceptedHeavyHadronNumber = 0;
    for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        // TODO: after all the checks remove "="
        if((*hadronMasses)[hadronIndex]<=m_clusterMass)
        {
            weight = (*hadronSamplingWeights)[hadronIndex];
            m_heavyHadronSamplingProbabilities[m_acceptedHeavyHadronNumber] = weight;
            m_availableHeavyHadronPositions[m_acceptedHeavyHadronNumber] = hadronIndex;
            weightSum += weight;
            ++m_acceptedHeavyHadronNumber;
        }
    }

    if(m_acceptedHeavyHadronNumber==0)
    {
        // TODO: find final strategy for anomalous flown, error code or exception?
        return true;
    }
    
    // Normalize sampling probabilities
    double cumulativeProbability = 0.;
    vector<double>* samplingCumulativeProbabilities(&m_samplingCumulativeProbabilities[m_heavyHadronSamplingGroup]);
    samplingCumulativeProbabilities->resize(m_acceptedHeavyHadronNumber,0.);
    for(unsigned int hadronIndex = 0;hadronIndex<m_acceptedHeavyHadronNumber;++hadronIndex)
    {
        weight = m_heavyHadronSamplingProbabilities[hadronIndex]/weightSum;
        m_heavyHadronSamplingProbabilities[hadronIndex] = weight;
        cumulativeProbability += weight;
        (*samplingCumulativeProbabilities)[hadronIndex] = cumulativeProbability;
    }
    
    return false;
}

void HadronSampling::computeSampledChannelWeight(const double i_clusterEnergyDensity)
{
    try
    {
        // TODO: add detailed description of weight calculation
        
        // Light flavored hadron sampling weight
        int numberOfLightHadrons = 0;
        for(map<HadronGroup,vector<unsigned int> >::iterator hadronGroupIt = m_hadronSetPositionMap.begin();
            hadronGroupIt->first<=LightNeutralMesons;++hadronGroupIt)
        {
            const unsigned int numberOfHadrons(hadronGroupIt->second.size());
            if(numberOfHadrons>0)
            {
                // Retrieve stored weight factors
                const vector<double>* lightHadronMultiplicityWeightFactors(&m_lightHadronMultiplicityWeightFactors[hadronGroupIt->first]);
                
                // Reorder hadron positions
                sort(hadronGroupIt->second.begin(),hadronGroupIt->second.end());
                
                // Count multiplicity and update weight
                map<unsigned int, int> hadronMultiplicityMap;
                unsigned int referencePosition = hadronGroupIt->second[0];
                hadronMultiplicityMap[referencePosition] = 1;
                for(unsigned int hadronIndex = 1;hadronIndex<numberOfHadrons;++hadronIndex)
                {
                    if(hadronGroupIt->second[hadronIndex] == referencePosition)
                    {
                        // Update number of identical hadron
                        ++hadronMultiplicityMap[referencePosition];
                    }
                    else
                    {                        
                        // Set new reference hadron position
                        referencePosition = hadronGroupIt->second[hadronIndex];
                        
                        // Reset hadron multiplicity
                        hadronMultiplicityMap[referencePosition] = 1;
                    }
                }
                
                // Update sampling weight
                const map<unsigned int,int>::const_iterator mapEnd(hadronMultiplicityMap.end());
                for(map<unsigned int,int>::const_iterator mapIt = hadronMultiplicityMap.begin(); mapIt!=mapEnd; ++mapIt )
                {
                    m_weight *= pow((*lightHadronMultiplicityWeightFactors)[mapIt->first],mapIt->second);
                }

                
                // Update number of light hadrons
                numberOfLightHadrons += numberOfHadrons;
            }
        }
        
        // Update sampling
        m_weight *= pow(m_lightHadronExponentialWeightFactor,m_clusterSamplingVolume);
        // TODO: check if division with positive exponent is faster
        m_weight /= pow(m_clusterSamplingVolume,numberOfLightHadrons);
        m_weight *= pow(m_clusterMass*gevCubeToFermiCubeOnTwoPiCube/i_clusterEnergyDensity,static_cast<int>(m_hadronSet.size()));
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

double HadronSampling::computeLightHadronMeanMultiplicity(const double i_hadronMass,const unsigned int i_hadronSpinMultiplicity) const
{
    // Light hadron mean multiplicities, used for the sampling procedure, are defined using the
    // corresponding grandcanonical equation (except for cluster volume, which is included at the
    // end of the sampling procedure for performance related reasons)
	using TMath::BesselK;
	return (i_hadronSpinMultiplicity/(2.*piSquared))*i_hadronMass*
    i_hadronMass*m_samplingTemperature*
    BesselK(2,i_hadronMass/m_samplingTemperature);
}

double HadronSampling::computeHeavyHadronSamplingWeight(const double i_hadronMass) const
{
    // TODO: why do we need this sampling? Heavy hadrons are independently sampled
	return std::exp(-i_hadronMass/m_samplingTemperature);
}

void HadronSampling::setHadronSamplingTemperature(const double i_hadronSamplingTemperature)
{
    m_isHadronSamplingReady = false;
    if(i_hadronSamplingTemperature<=0.)
    {
        throw HadronizationException("Error in hadron sampling setup, non positive sampling temperature pareameter provided",
                                     __FUNCTION__,111);
    }
    else
    {
        m_samplingTemperature = i_hadronSamplingTemperature;
        computeSamplingData();
        m_isHadronSamplingReady = true;
    }
}

void HadronSampling::setHadronSamplingEnergyDensity(const double i_hadronSamplingEnergyDensity)
{
    m_isHadronSamplingReady = false;
    if(i_hadronSamplingEnergyDensity<=0.)
    {
        throw HadronizationException("Error in hadron sampling setup, non positive sampling energy density provided",
                                     __FUNCTION__,112);
    }
    else
    {
        m_samplingEnergyDensity = i_hadronSamplingEnergyDensity/gevCubeToFermiCube;
        computeSamplingData();
        m_isHadronSamplingReady = true;
    }
}

// TODO: Test
void HadronSampling::computeFugacities(const double i_electricCharge,
                                       const int i_strangeCharge,
                                       const double i_baryonicCharge,
                                       double& io_lambdaElectricCharge,
                                       double& io_lambdaStrangeCharge,
                                       double& io_lambdaBaryonicCharge,
                                       double& io_samplingT)
{
#if 1
    
    #if 0
    
        #if 0

            m_electricCharge = i_electricCharge;
            m_strangeCharge = static_cast<double>(i_strangeCharge);
            m_baryonicCharge = i_baryonicCharge;
            
            const gsl_multiroot_fsolver_type *T;
            gsl_multiroot_fsolver *s;
            
            int status;
            size_t i, iter = 0;
            
            const size_t n = 4;
            struct rparams p = {m_lightHadronSpinMultiplicities,
                                m_lightHadronMasses,
                                m_lightHadronElectricCharges,
                                m_lightHadronBaryonicCharges,
                                m_lightHadronStrangeCharges,
                                m_clusterSamplingVolume,
                                m_residualMass,
                                m_electricCharge,
                                m_baryonicCharge,
                                m_strangeCharge};
            
            cout<<"m_residualMass "<<m_residualMass<<endl;
            cout<<"m_clusterSamplingVolume "<<m_clusterSamplingVolume<<endl;
            cout<<"m_electricCharge "<<m_electricCharge<<endl;
            cout<<"m_baryonicCharge "<<m_baryonicCharge<<endl;
            cout<<"m_strangeCharge "<<m_strangeCharge<<endl;
    
            gsl_multiroot_function f = {&computeSaddlePointEquations4,
                                        n,
                                        &p};
            
            double x_init[4] = {1.0,1.0,1.0,0.160};
    //        double x_init[4] = {1.4,18.0,0.7,0.15};
        
        //    x[1]=18.0; B
        //    x[2]=0.7; S
        //    x[3]=1.4 Q
            
            gsl_vector *x = gsl_vector_alloc(n);
            
            gsl_vector_set(x, 0, x_init[0]);
            gsl_vector_set(x, 1, x_init[1]);
            gsl_vector_set(x, 2, x_init[2]);
            gsl_vector_set(x, 3, x_init[3]);
            
            T = gsl_multiroot_fsolver_hybrids;
            s = gsl_multiroot_fsolver_alloc(T,
                                            n);
            gsl_multiroot_fsolver_set(s,
                                      &f,
                                      x);
            
            print_state4(iter, s);
            
            do
            {
                iter++;
                status = gsl_multiroot_fsolver_iterate(s);
                
                print_state4(iter,
                            s);
                
                if (status)   /* check if solver is stuck */
                    break;
                
                status =
                gsl_multiroot_test_residual (s->f, 1e-9);
            }
            while (status == GSL_CONTINUE && iter < 1000);
            
            printf ("status = %s\n", gsl_strerror (status));
            
            gsl_multiroot_fsolver_free(s);
            gsl_vector_free(x);
            return;
    
        #else

            m_electricCharge = i_electricCharge;
            m_strangeCharge = static_cast<double>(i_strangeCharge);
            m_baryonicCharge = i_baryonicCharge;
            
            const gsl_multiroot_fsolver_type *T;
            gsl_multiroot_fsolver *s;
            
            int status;
            size_t i, iter = 0;
            
            const size_t n = 4;
            struct rparams p = {m_lightHadronSpinMultiplicities,
                m_lightHadronMasses,
                m_lightHadronElectricCharges,
                m_lightHadronBaryonicCharges,
                m_lightHadronStrangeCharges,
                m_clusterSamplingVolume,
                m_residualMass,
                m_electricCharge,
                m_baryonicCharge,
                m_strangeCharge};
            
            cout<<"m_residualMass "<<m_residualMass<<endl;
            cout<<"m_clusterSamplingVolume "<<m_clusterSamplingVolume<<endl;
            cout<<"m_electricCharge "<<m_electricCharge<<endl;
            cout<<"m_baryonicCharge "<<m_baryonicCharge<<endl;
            cout<<"m_strangeCharge "<<m_strangeCharge<<endl;
    
            gsl_multiroot_function f = {&computeSaddlePointEquations4Ex,
                n,
                &p};
            
            double x_init[4] = {0.,0.,0.,0.160};
            //        double x_init[4] = {1.4,18.0,0.7,0.15};
            
            //    x[1]=18.0; B
            //    x[2]=0.7; S
            //    x[3]=1.4 Q
            
            gsl_vector *x = gsl_vector_alloc(n);
            
            gsl_vector_set(x, 0, x_init[0]);
            gsl_vector_set(x, 1, x_init[1]);
            gsl_vector_set(x, 2, x_init[2]);
            gsl_vector_set(x, 3, x_init[3]);
            
            T = gsl_multiroot_fsolver_hybrids;
            s = gsl_multiroot_fsolver_alloc(T,
                                            n);
            gsl_multiroot_fsolver_set(s,
                                      &f,
                                      x);
            
            print_state4Ex(iter, s);
            
            do
            {
                iter++;
                status = gsl_multiroot_fsolver_iterate(s);
                
                print_state4Ex(iter,
                             s);
                
                if (status)   /* check if solver is stuck */
                    break;
                
                status =
                gsl_multiroot_test_residual (s->f, 1e-9);
            }
            while (status == GSL_CONTINUE && iter < 1000);
            
            printf ("status = %s\n", gsl_strerror (status));
            
            gsl_multiroot_fsolver_free(s);
            gsl_vector_free(x);
            return;

        #endif
    
    #else
    
        #if 0
    
            m_electricCharge = i_electricCharge;
            m_strangeCharge = static_cast<double>(i_strangeCharge);
            m_baryonicCharge = i_baryonicCharge;
        
            const gsl_multiroot_fdfsolver_type *T;
            gsl_multiroot_fdfsolver *s;
        
            int status;
            size_t i, iter = 0;
        
            const size_t n = 4;
            struct rparams p = {m_lightHadronSpinMultiplicities,
                m_lightHadronMasses,
                m_lightHadronElectricCharges,
                m_lightHadronBaryonicCharges,
                m_lightHadronStrangeCharges,
                m_clusterSamplingVolume,
                m_residualMass,
                m_electricCharge,
                m_baryonicCharge,
                m_strangeCharge};
            
            cout<<"m_residualMass "<<m_residualMass<<endl;
            cout<<"m_clusterSamplingVolume "<<m_clusterSamplingVolume<<endl;
            cout<<"m_electricCharge "<<m_electricCharge<<endl;
            cout<<"m_baryonicCharge "<<m_baryonicCharge<<endl;
            cout<<"m_strangeCharge "<<m_strangeCharge<<endl;
    
            gsl_multiroot_function_fdf f = {&computeSaddlePointEquations4,
                &computeSaddlePointEquations4_df,
                &computeSaddlePointEquations4_fdf,
                n, &p};
        
//            double x_init[4] = {1.0,1.0,1.0,0.160};
            double x_init[4] = {17.424,0.238,0.238,0.122};
    
            //    x[1]=18.0; B
            //    x[2]=0.7; S
            //    x[3]=1.4 Q
            
            gsl_vector *x = gsl_vector_alloc(n);
            
            gsl_vector_set(x, 0, x_init[0]);
            gsl_vector_set(x, 1, x_init[1]);
            gsl_vector_set(x, 2, x_init[2]);
            gsl_vector_set(x, 3, x_init[3]);
            
            T =  gsl_multiroot_fdfsolver_hybridsj;
            s = gsl_multiroot_fdfsolver_alloc (T,
                                            n);
            gsl_multiroot_fdfsolver_set(s,
                                      &f,
                                      x);
            
            print_state4(iter, s);
            
            do
            {
                iter++;
                status = gsl_multiroot_fdfsolver_iterate(s);
                
                print_state4(iter,
                             s);
                
                if (status)   /* check if solver is stuck */
                    break;
                
                status =
                gsl_multiroot_test_residual (s->f, 1e-9);
            }
            while (status == GSL_CONTINUE && iter < 1000);
            
            printf ("status = %s\n", gsl_strerror (status));
            
            gsl_multiroot_fdfsolver_free(s);
            gsl_vector_free(x);
            return;
    
        #else

            m_electricCharge = i_electricCharge;
            m_strangeCharge = static_cast<double>(i_strangeCharge);
            m_baryonicCharge = i_baryonicCharge;
            
            const gsl_multiroot_fdfsolver_type *T;
            gsl_multiroot_fdfsolver *s;
            
            int status;
            size_t i, iter = 0;
            
            const size_t n = 4;
            struct rparams p = {m_lightHadronSpinMultiplicities,
                m_lightHadronMasses,
                m_lightHadronElectricCharges,
                m_lightHadronBaryonicCharges,
                m_lightHadronStrangeCharges,
                m_clusterSamplingVolume,
                m_residualMass,
                m_electricCharge,
                m_baryonicCharge,
                m_strangeCharge};
    
            cout<<"m_residualMass "<<m_residualMass<<endl;
            cout<<"m_clusterSamplingVolume "<<m_clusterSamplingVolume<<endl;
            cout<<"m_electricCharge "<<m_electricCharge<<endl;
            cout<<"m_baryonicCharge "<<m_baryonicCharge<<endl;
            cout<<"m_strangeCharge "<<m_strangeCharge<<endl;
    
            gsl_multiroot_function_fdf f = {&computeSaddlePointEquations4Ex,
                &computeSaddlePointEquations4Ex_df,
                &computeSaddlePointEquations4Ex_fdf,
                n, &p};
            
    double x_init[4] = {std::log(1.0),std::log(1.0),std::log(1.0),0.160};
//            double x_init[4] = {17.424,0.238,0.238,0.122};
    
            //    x[1]=18.0; B
            //    x[2]=0.7; S
            //    x[3]=1.4 Q
            
            gsl_vector *x = gsl_vector_alloc(n);
            
            gsl_vector_set(x, 0, x_init[0]);
            gsl_vector_set(x, 1, x_init[1]);
            gsl_vector_set(x, 2, x_init[2]);
            gsl_vector_set(x, 3, x_init[3]);
            
            T =  gsl_multiroot_fdfsolver_hybridsj;
            s = gsl_multiroot_fdfsolver_alloc (T,
                                               n);
            gsl_multiroot_fdfsolver_set(s,
                                        &f,
                                        x);
            
            print_state4Ex(iter, s);
            
            do
            {
                iter++;
                status = gsl_multiroot_fdfsolver_iterate(s);
                
                print_state4Ex(iter,
                             s);
                
                if (status)   /* check if solver is stuck */
                    break;
                
                status =
                gsl_multiroot_test_residual (s->f, 1e-9);
            }
            while (status == GSL_CONTINUE && iter < 1000);
            
            printf ("status = %s\n", gsl_strerror (status));
            
            gsl_multiroot_fdfsolver_free(s);
            gsl_vector_free(x);
            return;

        #endif
    
    #endif
    
#else
    
    m_electricCharge = i_electricCharge;
    m_strangeCharge = static_cast<double>(i_strangeCharge);
    m_baryonicCharge = i_baryonicCharge;
    
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t i, iter = 0;
    
    //    const size_t n = 4;
    const size_t n = 3;
    struct rparams p = {m_lightHadronSpinMultiplicities,
        m_lightHadronMasses,
        m_lightHadronElectricCharges,
        m_lightHadronBaryonicCharges,
        m_lightHadronStrangeCharges,
        m_clusterSamplingVolume,
        m_residualMass,
        m_electricCharge,
        m_baryonicCharge,
        m_strangeCharge};
    
    gsl_multiroot_function f = {&computeSaddlePointEquations3,
                                n,
                                &p};
    
    double x_init[3] = {1.0,1.0,1.0};
    gsl_vector *x = gsl_vector_alloc(n);
    
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);
    gsl_vector_set(x, 2, x_init[2]);
    
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,
                                    n);
    gsl_multiroot_fsolver_set(s,
                              &f,
                              x);
    
    print_state3(iter, s);
    
    do
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        
        print_state3(iter,
                    s);
        
        if (status)   /* check if solver is stuck */
            break;
        
        status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
    while (status == GSL_CONTINUE && iter < 1000);
    
    printf ("status = %s\n", gsl_strerror (status));
    
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    return;
    
#endif
}
// TODO: Test

