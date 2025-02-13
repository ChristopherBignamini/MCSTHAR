#include "../Include/HadronSamplingNew.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/Constants.h"
#include <cassert>
#include "TMath.h"

// TODO: Test
#include <iostream>
#include <iomanip>
#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_sf_bessel.h"
struct rparams
{
    vector<double> hadronSpinMultiplicities;
    vector<double> hadronMasses;
    vector<double> hadronElectricCharges;
    vector<double> hadronBaryonicCharges;
    vector<double> hadronStrangeCharges;
    vector<double> hadronCharmCharges;
    vector<double> hadronBottomCharges;
    double clusterSamplingVolume;
    double clusterMass;
    double electricCharge;
    double baryonicCharge;
    double strangeCharge;
    double charmCharge;
    double bottomCharge;
};

void print_state(unsigned int iter, gsl_multiroot_fdfsolver * s)
{
    cout<<setprecision(7);
//    cout<<"iter = "<<iter<<" x = "
//        <<gsl_vector_get(s->x, 0)<<" "
//        <<gsl_vector_get(s->x, 1)<<" "
//        <<gsl_vector_get(s->x, 2)<<" "
//        <<gsl_vector_get(s->x, 3)<<" "
//        <<gsl_vector_get(s->x, 4)<<" "
//        <<gsl_vector_get(s->x, 5)<<endl<<" f(x) = "
//        <<gsl_vector_get (s->f, 0)<<" "
//        <<gsl_vector_get (s->f, 1)<<" "
//        <<gsl_vector_get (s->f, 2)<<" "
//        <<gsl_vector_get (s->f, 3)<<" "
//        <<gsl_vector_get (s->f, 4)<<" "
//        <<gsl_vector_get (s->f, 5)<<endl;
    cout<<"iter = "<<iter<<" x = "
    <<std::exp(gsl_vector_get(s->x, 0))<<" "
    <<std::exp(gsl_vector_get(s->x, 1))<<" "
    <<std::exp(gsl_vector_get(s->x, 2))<<" "
    <<std::exp(gsl_vector_get(s->x, 3))<<" "
    <<std::exp(gsl_vector_get(s->x, 4))<<" "
    <<gsl_vector_get(s->x, 5)<<" f(x) = "
    <<gsl_vector_get (s->f, 0)<<" "
    <<gsl_vector_get (s->f, 1)<<" "
    <<gsl_vector_get (s->f, 2)<<" "
    <<gsl_vector_get (s->f, 3)<<" "
    <<gsl_vector_get (s->f, 4)<<" "
    <<gsl_vector_get (s->f, 5)<<endl;
    
    
//    printf ("iter = %3u x = % .6f % .6f % .6f % .6f % .6f % .6f "
//            "f(x) = % .6e % .6e % .6e % .6e % .6e % .6e\n",
//            iter,
//            std::exp(gsl_vector_get (s->x, 0)),
//            std::exp(gsl_vector_get (s->x, 1)),
//            std::exp(gsl_vector_get (s->x, 2)),
//            std::exp(gsl_vector_get (s->x, 3)),
//            std::exp(gsl_vector_get (s->x, 4)),
//            gsl_vector_get (s->x, 5),
//            gsl_vector_get (s->f, 0),
//            gsl_vector_get (s->f, 1),
//            gsl_vector_get (s->f, 2),
//            gsl_vector_get (s->f, 3),
//            gsl_vector_get (s->f, 4),
//            gsl_vector_get (s->f, 5));
}


int computeSaddlePointEquations(const gsl_vector* x,
                                void* params,
                                gsl_vector* f)
{
//    cout<<endl;
//    cout<<"computeSaddlePointEquations"<<endl;
    const vector<double> hadronSpinMultiplicities(((struct rparams *) params)->hadronSpinMultiplicities);
    const vector<double> hadronMasses(((struct rparams *) params)->hadronMasses);
    const vector<double> hadronElectricCharges(((struct rparams *) params)->hadronElectricCharges);
    const vector<double> hadronBaryonicCharges(((struct rparams *) params)->hadronBaryonicCharges);
    const vector<double> hadronStrangeCharges(((struct rparams *) params)->hadronStrangeCharges);
    const vector<double> hadronCharmCharges(((struct rparams *) params)->hadronCharmCharges);
    const vector<double> hadronBottomCharges(((struct rparams *) params)->hadronBottomCharges);
    
    const double muQ(gsl_vector_get(x,0));
    const double muB(gsl_vector_get(x,1));
    const double muS(gsl_vector_get(x,2));
    const double muC(gsl_vector_get(x,3));
    const double muBt(gsl_vector_get(x,4));
    const double T(gsl_vector_get(x,5));
    
//    cout<<"muQ "<<muQ<<endl;
//    cout<<"muB "<<muB<<endl;
//    cout<<"muS "<<muS<<endl;
//    cout<<"muC "<<muC<<endl;
//    cout<<"muBt "<<muBt<<endl;
//    cout<<"T "<<T<<endl;
    
    const double samplingVolumeTimesTOn2piSquared(((struct rparams *) params)->clusterSamplingVolume*T/twoPiSquared);
    const double samplingVolumeTimesSquaredTOn2piSquared(samplingVolumeTimesTOn2piSquared*T);
    
    double yElectricCharge(0.);
    double yBaryonicCharge(0.);
    double yStrangeCharge(0.);
    double yCharmCharge(0.);
    double yBottomCharge(0.);
    double mass(0.);
    
    const unsigned int numberOfHadrons(hadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(hadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(hadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(hadronStrangeCharges[hadronIndex]);
        const double& currentCharmCharge(hadronCharmCharges[hadronIndex]);
        const double& currentBottomCharge(hadronBottomCharges[hadronIndex]);
        const double& currentSpinMultiplicity(hadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(hadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMass(currentMass*currentMass);
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = std::exp(muQ*currentElectricCharge)*
                       std::exp(muB*currentBaryonicCharge)*
                       std::exp(muS*currentStrangeCharge)*
                       std::exp(muC*currentCharmCharge)*
                       std::exp(muBt*currentBottomCharge);
        
//        cout<<"currentBottomCharge "<<currentBottomCharge
//            <<" lambdaFactor "<<lambdaFactor
//            <<" massOnT*gsl_sf_bessel_Kn(1,massOnT) "<<massOnT*gsl_sf_bessel_Kn(1,massOnT)
//            <<" bessel2 "<<bessel2<<endl;
        
        yElectricCharge += currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;// TODO: factorize
        yBaryonicCharge += currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        yStrangeCharge += currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        yCharmCharge += currentCharmCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        yBottomCharge += currentBottomCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor;
        mass += squaredMass*currentSpinMultiplicity*(massOnT*gsl_sf_bessel_Kn(1,massOnT)+3.*bessel2)*lambdaFactor;
    }
    
    yElectricCharge = ((struct rparams *) params)->electricCharge -
                      yElectricCharge*samplingVolumeTimesTOn2piSquared;
    yBaryonicCharge = ((struct rparams *) params)->baryonicCharge -
                      yBaryonicCharge*samplingVolumeTimesTOn2piSquared;
    yStrangeCharge = ((struct rparams *) params)->strangeCharge -
                     yStrangeCharge*samplingVolumeTimesTOn2piSquared;
    yCharmCharge = ((struct rparams *) params)->charmCharge -
                   yCharmCharge*samplingVolumeTimesTOn2piSquared;
    yBottomCharge = ((struct rparams *) params)->bottomCharge -
                    yBottomCharge*samplingVolumeTimesTOn2piSquared;
    mass = ((struct rparams *) params)->clusterMass -
           mass*samplingVolumeTimesSquaredTOn2piSquared;
    
    
//    cout<<"yElectricCharge "<<yElectricCharge<<endl;
//    cout<<"yBaryonicCharge "<<yBaryonicCharge<<endl;
//    cout<<"yStrangeCharge "<<yStrangeCharge<<endl;
//    cout<<" yCharmCharge "<<yCharmCharge<<endl;
//    cout<<" yBottomCharge "<<yBottomCharge<<endl;
//    cout<<" mass "<<mass<<endl;
    
    gsl_vector_set(f,0,yElectricCharge);
    f->data[0*f->stride] = yElectricCharge;
    gsl_vector_set(f,1,yBaryonicCharge);
    f->data[1*f->stride] = yBaryonicCharge;
    gsl_vector_set(f,2,yStrangeCharge);
    f->data[2*f->stride] = yStrangeCharge;
    gsl_vector_set(f,3,yCharmCharge);
    f->data[3*f->stride] = yCharmCharge;
    gsl_vector_set(f,4,yBottomCharge);
    f->data[4*f->stride] = yBottomCharge;
    gsl_vector_set(f,5,mass);
    f->data[5*f->stride] = mass;
    
    
    return GSL_SUCCESS;
}

int computeSaddlePointEquations_df(const gsl_vector* x,
                                   void* params,
                                   gsl_matrix* J)
{
//    cout<<"computeSaddlePointEquations_df"<<endl;

    const vector<double> hadronSpinMultiplicities(((struct rparams *) params)->hadronSpinMultiplicities);
    const vector<double> hadronMasses(((struct rparams *) params)->hadronMasses);
    const vector<double> hadronElectricCharges(((struct rparams *) params)->hadronElectricCharges);
    const vector<double> hadronBaryonicCharges(((struct rparams *) params)->hadronBaryonicCharges);
    const vector<double> hadronStrangeCharges(((struct rparams *) params)->hadronStrangeCharges);
    const vector<double> hadronCharmCharges(((struct rparams *) params)->hadronCharmCharges);
    const vector<double> hadronBottomCharges(((struct rparams *) params)->hadronBottomCharges);

    
    const double muQ(gsl_vector_get(x,0));
    const double muB(gsl_vector_get(x,1));
    const double muS(gsl_vector_get(x,2));
    const double muC(gsl_vector_get(x,3));
    const double muBt(gsl_vector_get(x,4));
    const double T(gsl_vector_get(x,5));
    
//    cout<<"muQ "<<muQ<<endl;
//    cout<<"muB "<<muB<<endl;
//    cout<<"muS "<<muS<<endl;
//    cout<<"muC "<<muC<<endl;
//    cout<<"muBt "<<muBt<<endl;
//    cout<<"T "<<T<<endl;
    
    const double samplingVolumeOn2piSquared(((struct rparams *) params)->clusterSamplingVolume/twoPiSquared);
    const double samplingVolumeTimesTOn2piSquared(samplingVolumeOn2piSquared*T);
    const double samplingVolumeTimesSquaredTOn2piSquared(samplingVolumeTimesTOn2piSquared*T);
    
    
    double dElectricChargeDLambdaQ(0.);
    double dElectricChargeDLambdaB(0.);
    double dElectricChargeDLambdaS(0.);
    double dElectricChargeDLambdaC(0.);
    double dElectricChargeDLambdaBt(0.);
    double dElectricChargeDT(0.);
    
    double dBaryonicChargeDLambdaQ(0.);
    double dBaryonicChargeDLambdaB(0.);
    double dBaryonicChargeDLambdaS(0.);
    double dBaryonicChargeDLambdaC(0.);
    double dBaryonicChargeDLambdaBt(0.);
    double dBaryonicChargeDT(0.);
    
    double dStrangeChargeDLambdaQ(0.);
    double dStrangeChargeDLambdaB(0.);
    double dStrangeChargeDLambdaS(0.);
    double dStrangeChargeDLambdaC(0.);
    double dStrangeChargeDLambdaBt(0.);
    double dStrangeChargeDT(0.);
    
    double dCharmChargeDLambdaQ(0.);
    double dCharmChargeDLambdaB(0.);
    double dCharmChargeDLambdaS(0.);
    double dCharmChargeDLambdaC(0.);
    double dCharmChargeDLambdaBt(0.);
    double dCharmChargeDT(0.);
    
    double dBottomChargeDLambdaQ(0.);
    double dBottomChargeDLambdaB(0.);
    double dBottomChargeDLambdaS(0.);
    double dBottomChargeDLambdaC(0.);
    double dBottomChargeDLambdaBt(0.);
    double dBottomChargeDT(0.);
    
    double dMassDLambdaQ(0.);
    double dMassDLambdaB(0.);
    double dMassDLambdaS(0.);
    double dMassDLambdaC(0.);
    double dMassDLambdaBt(0.);
    double dMassDT(0.);
    
    const unsigned int numberOfHadrons(hadronSpinMultiplicities.size());
    double lambdaFactor(1.);
    
    for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        const double& currentElectricCharge(hadronElectricCharges[hadronIndex]);
        const double& currentBaryonicCharge(hadronBaryonicCharges[hadronIndex]);
        const double& currentStrangeCharge(hadronStrangeCharges[hadronIndex]);
        const double& currentCharmCharge(hadronCharmCharges[hadronIndex]);
        const double& currentBottomCharge(hadronBottomCharges[hadronIndex]);
        const double& currentSpinMultiplicity(hadronSpinMultiplicities[hadronIndex]);
        const double& currentMass(hadronMasses[hadronIndex]);
        
        const double massOnT(currentMass/T);
        const double squaredMassOnT(massOnT*massOnT);
        const double squaredMass(currentMass*currentMass);
        const double cubeMass(squaredMass*currentMass);
        const double cubeMassOnT(cubeMass/T);
        const double bessel1(gsl_sf_bessel_Kn(1,massOnT));
        const double bessel2(gsl_sf_bessel_Kn(2,massOnT));
        
        lambdaFactor = std::exp(muQ*currentElectricCharge)*
                       std::exp(muB*currentBaryonicCharge)*
                       std::exp(muS*currentStrangeCharge)*
                       std::exp(muC*currentCharmCharge)*
                       std::exp(muBt*currentBottomCharge);
                
        dElectricChargeDLambdaQ +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;// TODO: factorize
        dElectricChargeDLambdaB +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dElectricChargeDLambdaS +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dElectricChargeDLambdaC +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentCharmCharge;
        dElectricChargeDLambdaBt +=
        currentElectricCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBottomCharge;
        dElectricChargeDT +=
        currentElectricCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        
        dBaryonicChargeDLambdaQ +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;
        dBaryonicChargeDLambdaB +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dBaryonicChargeDLambdaS +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dBaryonicChargeDLambdaC +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentCharmCharge;
        dBaryonicChargeDLambdaBt +=
        currentBaryonicCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBottomCharge;
        dBaryonicChargeDT +=
        currentBaryonicCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        
        dStrangeChargeDLambdaQ +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;
        dStrangeChargeDLambdaB +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dStrangeChargeDLambdaS +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dStrangeChargeDLambdaC +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentCharmCharge;
        dStrangeChargeDLambdaBt +=
        currentStrangeCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBottomCharge;
        dStrangeChargeDT +=
        currentStrangeCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;

        dCharmChargeDLambdaQ +=
        currentCharmCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;
        dCharmChargeDLambdaB +=
        currentCharmCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dCharmChargeDLambdaS +=
        currentCharmCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dCharmChargeDLambdaC +=
        currentCharmCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentCharmCharge;
        dCharmChargeDLambdaBt +=
        currentCharmCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBottomCharge;
        dCharmChargeDT +=
        currentCharmCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;
        
        dBottomChargeDLambdaQ +=
        currentBottomCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentElectricCharge;
        dBottomChargeDLambdaB +=
        currentBottomCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBaryonicCharge;
        dBottomChargeDLambdaS +=
        currentBottomCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentStrangeCharge;
        dBottomChargeDLambdaC +=
        currentBottomCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentCharmCharge;
        dBottomChargeDLambdaBt +=
        currentBottomCharge*currentSpinMultiplicity*squaredMass*bessel2*lambdaFactor*currentBottomCharge;
        dBottomChargeDT +=
        currentBottomCharge*currentSpinMultiplicity*(cubeMassOnT*bessel1+3.*squaredMass*bessel2)*lambdaFactor;        
        
        dMassDLambdaQ +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentElectricCharge;
        dMassDLambdaB +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentBaryonicCharge;
        dMassDLambdaS +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentStrangeCharge;
        dMassDLambdaC +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentCharmCharge;
        dMassDLambdaBt +=
        squaredMass*currentSpinMultiplicity*(massOnT*bessel1+3.*bessel2)*lambdaFactor*currentBottomCharge;
        dMassDT +=
        squaredMass*currentSpinMultiplicity*(2.*massOnT*bessel1+(12.+squaredMassOnT)*bessel2)*lambdaFactor;
    }
        
    dElectricChargeDLambdaQ = -dElectricChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaB = -dElectricChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaS = -dElectricChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaC = -dElectricChargeDLambdaC*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDLambdaBt = -dElectricChargeDLambdaBt*samplingVolumeTimesTOn2piSquared;
    dElectricChargeDT = -dElectricChargeDT*samplingVolumeOn2piSquared;
    
    dBaryonicChargeDLambdaQ = -dBaryonicChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaB = -dBaryonicChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaS = -dBaryonicChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaC = -dBaryonicChargeDLambdaC*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDLambdaBt = -dBaryonicChargeDLambdaBt*samplingVolumeTimesTOn2piSquared;
    dBaryonicChargeDT = -dBaryonicChargeDT*samplingVolumeOn2piSquared;
    
    dStrangeChargeDLambdaQ = -dStrangeChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaB = -dStrangeChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaS = -dStrangeChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaC = -dStrangeChargeDLambdaC*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDLambdaBt = -dStrangeChargeDLambdaBt*samplingVolumeTimesTOn2piSquared;
    dStrangeChargeDT = -dStrangeChargeDT*samplingVolumeOn2piSquared;
    
    dCharmChargeDLambdaQ = -dCharmChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dCharmChargeDLambdaB = -dCharmChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dCharmChargeDLambdaS = -dCharmChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dCharmChargeDLambdaC = -dCharmChargeDLambdaC*samplingVolumeTimesTOn2piSquared;
    dCharmChargeDLambdaBt = -dCharmChargeDLambdaBt*samplingVolumeTimesTOn2piSquared;
    dCharmChargeDT = -dCharmChargeDT*samplingVolumeOn2piSquared;

    dBottomChargeDLambdaQ = -dBottomChargeDLambdaQ*samplingVolumeTimesTOn2piSquared;
    dBottomChargeDLambdaB = -dBottomChargeDLambdaB*samplingVolumeTimesTOn2piSquared;
    dBottomChargeDLambdaS = -dBottomChargeDLambdaS*samplingVolumeTimesTOn2piSquared;
    dBottomChargeDLambdaC = -dBottomChargeDLambdaC*samplingVolumeTimesTOn2piSquared;
    dBottomChargeDLambdaBt = -dBottomChargeDLambdaBt*samplingVolumeTimesTOn2piSquared;
    dBottomChargeDT = -dBottomChargeDT*samplingVolumeOn2piSquared;

    dMassDLambdaQ = -dMassDLambdaQ*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaB = -dMassDLambdaB*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaS = -dMassDLambdaS*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaC = -dMassDLambdaC*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDLambdaBt = -dMassDLambdaBt*samplingVolumeTimesSquaredTOn2piSquared;
    dMassDT = -dMassDT*samplingVolumeTimesTOn2piSquared;
    
//    cout<<"dElectricChargeDLambdaQ "<<dElectricChargeDLambdaQ<<endl;
//    cout<<"dElectricChargeDLambdaB "<<dElectricChargeDLambdaB<<endl;
//    cout<<"dElectricChargeDLambdaS "<<dElectricChargeDLambdaS<<endl;
//    cout<<"dElectricChargeDLambdaC "<<dElectricChargeDLambdaC<<endl;
//    cout<<"dElectricChargeDLambdaBt "<<dElectricChargeDLambdaBt<<endl;
//    cout<<"dElectricChargeDT "<<dElectricChargeDT<<endl;
//    cout<<"dBaryonicChargeDLambdaQ "<<dBaryonicChargeDLambdaQ<<endl;
//    cout<<"dBaryonicChargeDLambdaB "<<dBaryonicChargeDLambdaB<<endl;
//    cout<<"dBaryonicChargeDLambdaS "<<dBaryonicChargeDLambdaS<<endl;
//    cout<<"dBaryonicChargeDLambdaC "<<dBaryonicChargeDLambdaC<<endl;
//    cout<<"dBaryonicChargeDLambdaBt "<<dBaryonicChargeDLambdaBt<<endl;
//    cout<<"dBaryonicChargeDT "<<dBaryonicChargeDT<<endl;
//    cout<<"dStrangeChargeDLambdaQ "<<dStrangeChargeDLambdaQ<<endl;
//    cout<<"dStrangeChargeDLambdaB "<<dStrangeChargeDLambdaB<<endl;
//    cout<<"dStrangeChargeDLambdaS "<<dStrangeChargeDLambdaS<<endl;
//    cout<<"dStrangeChargeDLambdaC "<<dStrangeChargeDLambdaC<<endl;
//    cout<<"dStrangeChargeDLambdaBt "<<dStrangeChargeDLambdaBt<<endl;
//    cout<<"dStrangeChargeDT "<<dStrangeChargeDT<<endl;
//    cout<<"dCharmChargeDLambdaQ "<<dCharmChargeDLambdaQ<<endl;
//    cout<<"dCharmChargeDLambdaB "<<dCharmChargeDLambdaB<<endl;
//    cout<<"dCharmChargeDLambdaS "<<dCharmChargeDLambdaS<<endl;
//    cout<<"dCharmChargeDLambdaC "<<dCharmChargeDLambdaC<<endl;
//    cout<<"dCharmChargeDLambdaBt "<<dCharmChargeDLambdaBt<<endl;
//    cout<<"dCharmChargeDT "<<dCharmChargeDT<<endl;
//    cout<<"dBottomChargeDLambdaQ "<<dBottomChargeDLambdaQ<<endl;
//    cout<<"dBottomChargeDLambdaB "<<dBottomChargeDLambdaB<<endl;
//    cout<<"dBottomChargeDLambdaS "<<dBottomChargeDLambdaS<<endl;
//    cout<<"dBottomChargeDLambdaC "<<dBottomChargeDLambdaC<<endl;
//    cout<<"dBottomChargeDLambdaBt "<<dBottomChargeDLambdaBt<<endl;
//    cout<<"dBottomChargeDT "<<dBottomChargeDT<<endl;
//    cout<<"dMassDLambdaQ "<<dMassDLambdaQ<<endl;
//    cout<<"dMassDLambdaB "<<dMassDLambdaB<<endl;
//    cout<<"dMassDLambdaS "<<dMassDLambdaS<<endl;
//    cout<<"dMassDLambdaC "<<dMassDLambdaC<<endl;
//    cout<<"dMassDLambdaBt "<<dMassDLambdaBt<<endl;
//    cout<<"dMassDT "<<dMassDT<<endl;
    
    gsl_matrix_set (J, 0, 0, dElectricChargeDLambdaQ);
    gsl_matrix_set (J, 0, 1, dElectricChargeDLambdaB);
    gsl_matrix_set (J, 0, 2, dElectricChargeDLambdaS);
    gsl_matrix_set (J, 0, 3, dElectricChargeDLambdaC);
    gsl_matrix_set (J, 0, 4, dElectricChargeDLambdaBt);
    gsl_matrix_set (J, 0, 5, dElectricChargeDT);
    
    gsl_matrix_set (J, 1, 0, dBaryonicChargeDLambdaQ);
    gsl_matrix_set (J, 1, 1, dBaryonicChargeDLambdaB);
    gsl_matrix_set (J, 1, 2, dBaryonicChargeDLambdaS);
    gsl_matrix_set (J, 1, 3, dBaryonicChargeDLambdaC);
    gsl_matrix_set (J, 1, 4, dBaryonicChargeDLambdaBt);
    gsl_matrix_set (J, 1, 5, dBaryonicChargeDT);
    
    gsl_matrix_set (J, 2, 0, dStrangeChargeDLambdaQ);
    gsl_matrix_set (J, 2, 1, dStrangeChargeDLambdaB);
    gsl_matrix_set (J, 2, 2, dStrangeChargeDLambdaS);
    gsl_matrix_set (J, 2, 3, dStrangeChargeDLambdaC);
    gsl_matrix_set (J, 2, 4, dStrangeChargeDLambdaBt);
    gsl_matrix_set (J, 2, 5, dStrangeChargeDT);

    gsl_matrix_set (J, 3, 0, dCharmChargeDLambdaQ);
    gsl_matrix_set (J, 3, 1, dCharmChargeDLambdaB);
    gsl_matrix_set (J, 3, 2, dCharmChargeDLambdaS);
    gsl_matrix_set (J, 3, 3, dCharmChargeDLambdaC);
    gsl_matrix_set (J, 3, 4, dCharmChargeDLambdaBt);
    gsl_matrix_set (J, 3, 5, dCharmChargeDT);

    gsl_matrix_set (J, 4, 0, dBottomChargeDLambdaQ);
    gsl_matrix_set (J, 4, 1, dBottomChargeDLambdaB);
    gsl_matrix_set (J, 4, 2, dBottomChargeDLambdaS);
    gsl_matrix_set (J, 4, 3, dBottomChargeDLambdaC);
    gsl_matrix_set (J, 4, 4, dBottomChargeDLambdaBt);
    gsl_matrix_set (J, 4, 5, dBottomChargeDT);

    gsl_matrix_set (J, 5, 0, dMassDLambdaQ);
    gsl_matrix_set (J, 5, 1, dMassDLambdaB);
    gsl_matrix_set (J, 5, 2, dMassDLambdaS);
    gsl_matrix_set (J, 5, 3, dMassDLambdaC);
    gsl_matrix_set (J, 5, 4, dMassDLambdaBt);
    gsl_matrix_set (J, 5, 5, dMassDT);
    
    return GSL_SUCCESS;
}

int computeSaddlePointEquations_fdf(const gsl_vector* x,
                                       void* params,
                                       gsl_vector* f,
                                       gsl_matrix* J)
{
    computeSaddlePointEquations(x,params,f);
    computeSaddlePointEquations_df(x,params,J);
    return GSL_SUCCESS;
}


HadronSamplingNew::HadronSamplingNew(const HadronSamplingGroups& i_hadronSamplingGroups,
                                     RandomNumberGenerator& io_randomGenerator)
                                    :m_isHadronSamplingReady(false)
                                    ,m_isHadronSetAvailable(false)
                                    ,m_weight(0.)
                                    ,m_hadronSamplingGroups(&i_hadronSamplingGroups)
                                    ,m_randomGenerator(&io_randomGenerator)
                                    ,m_samplingTemperature(0.160)
                                    ,m_clusterSamplingVolume(0.)
                                    ,m_energyDensity(0.)
                                    ,m_clusterMass(0.)
                                    ,m_residualStrangeCharge(0)
                                    ,m_residualElectricCharge(0.)
                                    ,m_residualBaryonicCharge(0.)
                                    ,m_residualCharmCharge(0.)
                                    ,m_residualBottomCharge(0.)
                                    ,m_hadronExponentialWeightFactor(0.)
                                    ,m_hadronMassSum(0.)
                                    ,m_lambdaElectricCharge(1.)// TEST
                                    ,m_lambdaStrangeCharge(1.)// TEST
                                    ,m_lambdaBaryonicCharge(1.)
                                    ,m_lambdaCharmCharge(1.)
                                    ,m_lambdaBottomCharge(1.)
                                    ,m_electricCharge(0.)
                                    ,m_strangeCharge(0.)
                                    ,m_baryonicCharge(0.)
                                    ,m_charmCharge(0.)
                                    ,m_bottomCharge(0.)
                                    ,m_verboseLog(false)


{
    try
    {
        if(m_randomGenerator==NULL)
        {
            throw HadronizationException("Error during phase space sampling object creation, NULL pointer to random number generator provided",
                                         __FUNCTION__,113);
        }

        retrieveHadronData();
        m_isHadronSamplingReady = true;

    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

HadronSamplingNew::~HadronSamplingNew(void)
{
}

void HadronSamplingNew::retrieveHadronData(void)
{
    // Retrieve hadron sampling group list
    const vector<HadronGroup> groupList(m_hadronSamplingGroups->getHadronGroupList());
    
    // Loop over sampling group
    for(vector<HadronGroup>::const_iterator hadronGroup = groupList.begin();
        hadronGroup != groupList.end();++hadronGroup)
    {
        const vector<const HadronData*> hadronData(m_hadronSamplingGroups->getHadronGroup(*hadronGroup));
        const unsigned int numberOfHadrons(hadronData.size());
        for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons;++hadronIndex)
        {
            m_hadronSpinMultiplicities.push_back(hadronData[hadronIndex]->getSpinMultiplicity());
            
            if(hadronData[hadronIndex]->getBottomCharge()!=0)
            {
                m_hadronMasses.push_back(hadronData[hadronIndex]->getMass());
            }
            else
            {
                m_hadronMasses.push_back(hadronData[hadronIndex]->getMass());
            }
            m_hadronElectricCharges.push_back(hadronData[hadronIndex]->getElectricCharge());
            m_hadronBaryonicCharges.push_back(hadronData[hadronIndex]->getBaryonicCharge());
            m_hadronStrangeCharges.push_back(static_cast<double>(hadronData[hadronIndex]->getStrangeCharge()));
            m_hadronCharmCharges.push_back(static_cast<double>(hadronData[hadronIndex]->getCharmCharge()));
            m_hadronBottomCharges.push_back(static_cast<double>(hadronData[hadronIndex]->getBottomCharge()));
        }
    }
    
    return;
}

unsigned int HadronSamplingNew::run(const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration,
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
        
        // Compute cluster mass less phase space availability tolerance value (for channel acceptance test only)
        // Checking the decay channel total mass wrt this reduced cluster mass value, instead of the
        // real one, allows to include in the simulation and partition function calculation only channels with
        // a reasonable phase space integral value.
        const double clusterMassLessTolerance(i_mass-phaseSpaceEnergyAvailabilityTolerance);
        if(clusterMassLessTolerance<=0.)
        {
            throw HadronizationException("Error during hadron sampling, non positive cluster mass less phase space tolerance",
                                         __FUNCTION__,117);
        }
        
        if((m_clusterMass != i_mass) ||
           (m_strangeCharge != static_cast<double>(i_chargeConfiguration.strangeCharge)) ||
           (m_electricCharge != i_chargeConfiguration.electricCharge) ||
           (m_baryonicCharge != i_chargeConfiguration.baryonicCharge) ||
           (m_charmCharge != static_cast<double>(i_chargeConfiguration.charmCharge)) ||
           (m_bottomCharge != static_cast<double>(i_chargeConfiguration.bottomCharge)) ||
           (m_energyDensity != i_energyDensity))
        {
            m_clusterMass = i_mass;
            m_strangeCharge = static_cast<double>(i_chargeConfiguration.strangeCharge);
            m_electricCharge = i_chargeConfiguration.electricCharge;
            m_baryonicCharge = i_chargeConfiguration.baryonicCharge;
            m_charmCharge = static_cast<double>(i_chargeConfiguration.charmCharge);
            m_bottomCharge = static_cast<double>(i_chargeConfiguration.bottomCharge);
            m_energyDensity = i_energyDensity;
            m_clusterSamplingVolume = m_clusterMass/m_energyDensity*gevCubeToFermiCube;
            computeFugacities();
            computeSamplingData();
        }

        
        // TODO: check channel composition (e.g., does neutral channel contain ssbar particles??)
        // TODO: set parameter for sampling attempt count limit
        unsigned long long int samplingAttempt=0;// TODO: long long or long?
        unsigned long long int maxSamplingAttempt(100000);
        
        do
        {            
            m_residualStrangeCharge = m_strangeCharge;
            m_residualElectricCharge = m_electricCharge;
            m_residualBaryonicCharge = m_baryonicCharge;
            m_residualCharmCharge = m_charmCharge;
            m_residualBottomCharge = m_bottomCharge;
            
            m_weight = 1.;
            m_hadronMassSum = 0.;
            ++samplingAttempt;// TODO: see comment above concerning sampling attempt handling (generation)
                              // TODO: no sampling limit is foreseen for part. function calculation, check class usage case! 
            
            // TODO: there is no reason to remove the map keys I think.. Close the weight calculation
            //       open question and then optimize this (maybe old hadron create problem in new event generation)
            // TODO: make similar check for all other containers used
            m_hadronSet.clear();
            m_hadronSetPositionMap.clear();
                        
            // Run light flavored hadron sampling
            runHadronSampling();
            
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
        computeSampledChannelWeight();
        
        // Sampling correctly performed
        return 0;// TODO: use static table??
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void HadronSamplingNew::computeFugacities(void)
{
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;
    
    int status;
    size_t i, iter = 0;
    
    const size_t n = 6;
    struct rparams p = {m_hadronSpinMultiplicities,
                        m_hadronMasses,
                        m_hadronElectricCharges,
                        m_hadronBaryonicCharges,
                        m_hadronStrangeCharges,
                        m_hadronCharmCharges,
                        m_hadronBottomCharges,
                        m_clusterSamplingVolume,
                        m_clusterMass,
                        m_electricCharge,
                        m_baryonicCharge,
                        m_strangeCharge,
                        m_charmCharge,
                        m_bottomCharge};
        
    gsl_multiroot_function_fdf f = {&computeSaddlePointEquations,
                                    &computeSaddlePointEquations_df,
                                    &computeSaddlePointEquations_fdf,
                                    n,
                                    &p};
    
    double x_init[6] = {0.,0.,0.,0.,0., 1.0};
    
    gsl_vector *x = gsl_vector_alloc(n);
    
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);
    gsl_vector_set(x, 2, x_init[2]);
    gsl_vector_set(x, 3, x_init[3]);
    gsl_vector_set(x, 4, x_init[4]);
    gsl_vector_set(x, 5, x_init[5]);
    
    T =  gsl_multiroot_fdfsolver_gnewton;//gsl_multiroot_fdfsolver_gnewton;//gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc(T,
                                      n);
    gsl_multiroot_fdfsolver_set(s,
                                &f,
                                x);
    if(m_verboseLog)
    {
        print_state(iter, s);
    }
    
    do
    {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(s);
        
        if(m_verboseLog)
        {
            print_state(iter,s);
        }
        
        if (status)   /* check if solver is stuck */
            break;
        
        status = gsl_multiroot_test_residual (s->f, 1e-9);
    }
    while (status == GSL_CONTINUE && iter < 10000);
    
    if(m_verboseLog)
    {
        printf ("status = %s\n", gsl_strerror (status));
    }
    
    m_lambdaElectricCharge = std::exp(gsl_vector_get(s->x, 0));
    m_lambdaBaryonicCharge = std::exp(gsl_vector_get(s->x, 1));
    m_lambdaStrangeCharge = std::exp(gsl_vector_get(s->x, 2));
    m_lambdaCharmCharge = std::exp(gsl_vector_get(s->x, 3));
    m_lambdaBottomCharge = std::exp(gsl_vector_get(s->x, 4));
    m_samplingTemperature = gsl_vector_get(s->x, 5);
    
    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);
    return;
}
    
void HadronSamplingNew::computeSamplingData(void)
{
    try
    {
        m_hadronExponentialWeightFactor = 0.;
        
        // Retrieve hadron sampling group list
        const vector<HadronGroup> groupList(m_hadronSamplingGroups->getHadronGroupList());
                
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
                case CharmedHadrons:
                case AntiCharmedHadrons:
                case BottomedHadrons:
                case AntiBottomedHadrons:
                    computeHadronSamplingWeights(*hadronGroup);
                    break;
                default:
                    throw HadronizationException("Error during hadron sampling probability calculation, hadron group not identified",
                                                 __FUNCTION__,116);
            }
            
        }
        
        // Close light hadron exponential sampling factor calculation
        m_hadronExponentialWeightFactor = std::exp(m_hadronExponentialWeightFactor);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }

    return;
}
                                
void HadronSamplingNew::computeHadronSamplingWeights(const HadronGroup i_hadronGroup)
{
    try
    {
        const vector<const HadronData*> hadronData(m_hadronSamplingGroups->getHadronGroup(i_hadronGroup));
        const unsigned int numberOfHadrons(hadronData.size());
        double meanMultiplicity = 0.;
        double meanMultiplicitySum = 0.;
        unsigned int hadronSpinMultiplicity = 0;
        
        // Loop over hadrons and compute mean multiplicities
        vector<double> meanMultiplicities(numberOfHadrons,0.);
        vector<double> samplingWeightFactors(numberOfHadrons,0.);
        for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons;++hadronIndex)
        {
            hadronSpinMultiplicity = hadronData[hadronIndex]->getSpinMultiplicity();
            meanMultiplicity = computeHadronMeanMultiplicity(hadronData[hadronIndex]->getMass(),hadronSpinMultiplicity);
            
//            cout<<"el "<<hadronData[hadronIndex]->getElectricCharge()
//                <<" st "<<hadronData[hadronIndex]->getStrangeCharge()
//                <<" bar "<<hadronData[hadronIndex]->getBaryonicCharge()
//                <<" ch "<<hadronData[hadronIndex]->getCharmCharge()
//                <<" bt "<<hadronData[hadronIndex]->getBottomCharge()
//                <<" meanMultiplicity "<<meanMultiplicity
//            <<" lambda factor "<<pow(m_lambdaElectricCharge,hadronData[hadronIndex]->getElectricCharge())*
//            pow(m_lambdaStrangeCharge,hadronData[hadronIndex]->getStrangeCharge())*
//            pow(m_lambdaBaryonicCharge,hadronData[hadronIndex]->getBaryonicCharge())*
//            pow(m_lambdaCharmCharge,hadronData[hadronIndex]->getCharmCharge())*
//            pow(m_lambdaBottomCharge,hadronData[hadronIndex]->getBottomCharge())<<endl;
            
            meanMultiplicity *= pow(m_lambdaElectricCharge,hadronData[hadronIndex]->getElectricCharge())*
                                pow(m_lambdaStrangeCharge,hadronData[hadronIndex]->getStrangeCharge())*
                                pow(m_lambdaBaryonicCharge,hadronData[hadronIndex]->getBaryonicCharge())*
                                pow(m_lambdaCharmCharge,hadronData[hadronIndex]->getCharmCharge())*
                                pow(m_lambdaBottomCharge,hadronData[hadronIndex]->getBottomCharge());

//            cout<<" meanMultiplicity "<<meanMultiplicity<<endl;

            
            if(meanMultiplicity>0.)
            {
                meanMultiplicities[hadronIndex] = meanMultiplicity;
                meanMultiplicitySum += meanMultiplicity;
                m_hadronExponentialWeightFactor += meanMultiplicity;
                samplingWeightFactors[hadronIndex] = hadronSpinMultiplicity/meanMultiplicity;
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
        m_hadronMultiplicityWeightFactors[i_hadronGroup] = samplingWeightFactors;
        
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
    catch(HadronizationException& ex)
    {
        throw ex;
    }

}

double HadronSamplingNew::computeHadronMeanMultiplicity(const double i_hadronMass,const unsigned int i_hadronSpinMultiplicity) const
{
    // Light hadron mean multiplicities, used for the sampling procedure, are defined using the
    // corresponding grandcanonical equation (except for cluster volume, which is included at the
    // end of the sampling procedure for performance related reasons)
    using TMath::BesselK;
    return (i_hadronSpinMultiplicity/(2.*piSquared))*i_hadronMass*
            i_hadronMass*m_samplingTemperature*
            BesselK(2,i_hadronMass/m_samplingTemperature);
}
                               
void HadronSamplingNew::runHadronSampling(void)
{
    // TODO: check the possibility to perform a single random number generation call
    
    // TODO: mean multiplicities, even including sampling volume, are very small...why? Can we refine this???
    //       a lot of sampling attempts lead to zero hadrons at present!!!
    
    // TODO: can be the mutiplication per m_clusterSamplingVolume optimized
    // (in case of light flavored cluster m_clusterSamplingVolume is fixed during sampling attempts)
    
    try
    {
        // Generate the number of Bottomed and Antibottomed hadrons
        int numberOfHadron =
        m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[BottomedHadrons]);
        int numberOfAntiHadron =
        m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[AntiBottomedHadrons]);
        // Check bottom charge conservation
        if(numberOfAntiHadron-numberOfHadron != m_residualBottomCharge)
        {
            // Bottom charge violation, re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(BottomedHadrons,numberOfHadron);
        runHadronSampling(AntiBottomedHadrons,numberOfAntiHadron);
        
        
        
        // Generate the number of Charmed and Anticharmed hadrons
        numberOfHadron =
        m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[CharmedHadrons]);
        numberOfAntiHadron =
        m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[AntiCharmedHadrons]);
        // Check baryonic charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualCharmCharge)
        {
            // Bottom charge violation, re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(CharmedHadrons,numberOfHadron);
        runHadronSampling(AntiCharmedHadrons,numberOfAntiHadron);
        
        
        
        // Generate the number of baryons and antibaryons
        numberOfHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightBaryons]);
        numberOfAntiHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightAntibaryons]);
        // Check baryonic charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualBaryonicCharge)
        {
            // Baryonic charge violation, re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(LightBaryons,numberOfHadron);
        runHadronSampling(LightAntibaryons,numberOfAntiHadron);
        

        
        
        // Generate the number of strange mesons and antimesons
        numberOfHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightStrangeMesons]);
        numberOfAntiHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightStrangeAntimesons]);
        // Check strange charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualStrangeCharge)
        {
            // Strange charge violation, re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(LightStrangeMesons,numberOfHadron);
        runHadronSampling(LightStrangeAntimesons,numberOfAntiHadron);
                
        
        // Generate the number of electrically charged non strange mesons and antimesons
        numberOfHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightChargedMesons]);
        numberOfAntiHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightChargedAntimesons]);
        // Check electric charge conservation
        if(numberOfHadron-numberOfAntiHadron != m_residualElectricCharge)
        {
            // Strange charge violation re-start sampling
            return;
        }
        // Run hadron/antihadron sampling
        runHadronSampling(LightChargedMesons,numberOfHadron);
        runHadronSampling(LightChargedAntimesons,numberOfAntiHadron);
                
        
        // Generate the number of neutral mesons
        numberOfHadron =
            m_randomGenerator->getPoissonRandom(m_clusterSamplingVolume*m_hadronGroupMeanMultiplicities[LightNeutralMesons]);
        // TODO: avoid local
        // Run hadron/antihadron sampling
        runHadronSampling(LightNeutralMesons,numberOfHadron);
        
          // Hadron sampling completed
        m_isHadronSetAvailable = true;
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

void HadronSamplingNew::runHadronSampling(const HadronGroup i_hadronGroup,const unsigned int i_numberOfHadrons)
{
    try
    {
        // Check hadron availability
        assert(m_samplingCumulativeProbabilities[i_hadronGroup].size()>0);
        m_hadronSetPositionMap[i_hadronGroup].resize(i_numberOfHadrons);
        // Retrieve number of hadron minus one of the provided group
        unsigned int numberOfHadronsMinusOne(m_samplingCumulativeProbabilities[i_hadronGroup].size() - 1);
        
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
                        
            // Retrieve sampled hadron data
            const HadronData& sampledHadronData(*((m_hadronSamplingGroups->getHadronGroup(i_hadronGroup))[hadronIndex]));
                        
            // Update cluster residual data
            m_residualStrangeCharge -= sampledHadronData.getStrangeCharge();
            m_residualElectricCharge -= sampledHadronData.getElectricCharge();
            m_residualBaryonicCharge -= sampledHadronData.getBaryonicCharge();
            m_residualCharmCharge -= sampledHadronData.getCharmCharge();
            m_residualBottomCharge -= sampledHadronData.getBottomCharge();
            
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

void HadronSamplingNew::computeSampledChannelWeight(void)
{
    try
    {
        // TODO: add detailed description of weight calculation
        
        // Light flavored hadron sampling weight
        int totalNumberOfHadrons = 0;
        for(map<HadronGroup,vector<unsigned int> >::iterator hadronGroupIt = m_hadronSetPositionMap.begin();
            hadronGroupIt!=m_hadronSetPositionMap.end();++hadronGroupIt)
        {
            const unsigned int numberOfHadrons(hadronGroupIt->second.size());
            if(numberOfHadrons>0)
            {
                // Retrieve stored weight factors
                const vector<double>* hadronMultiplicityWeightFactors(&m_hadronMultiplicityWeightFactors[hadronGroupIt->first]);
                
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
                    m_weight *= pow((*hadronMultiplicityWeightFactors)[mapIt->first],mapIt->second);
                }
                
                // Update number of light hadrons
                totalNumberOfHadrons += numberOfHadrons;
            }
        }
        
        // Update sampling
        m_weight *= pow(m_hadronExponentialWeightFactor,m_clusterSamplingVolume);
        // TODO: check if division with positive exponent is faster
        m_weight /= pow(m_clusterSamplingVolume,totalNumberOfHadrons);
        m_weight *= pow(m_clusterMass*gevCubeToFermiCubeOnTwoPiCube/m_energyDensity,static_cast<int>(m_hadronSet.size()));
        
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

                               
// TODO: debug, after reordering removal this method must be const
const vector<HadronData>& HadronSamplingNew::getHadronSet(void)
{
    if(m_isHadronSetAvailable)
    {
        // TODO: debug, remove after comparison with old code
//        if(m_hadronSet.begin()->getCharmCharge() || m_hadronSet.begin()->getBottomCharge())
//        {
//            vector<HadronData> orderedHadronSet(m_hadronSet);
//            orderedHadronSet.erase(orderedHadronSet.begin());
//            orderedHadronSet.push_back(*m_hadronSet.begin());
//            m_hadronSet = orderedHadronSet;
//        }
        
        return m_hadronSet;
    }
    else
    {
        throw HadronizationException("Error during sampled hadron set retrieval, hadrons not available",
                                     __FUNCTION__,115);
    }
}

unsigned int HadronSamplingNew::getNumberOfHadrons(void) const
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

double HadronSamplingNew::getSamplingWeight(void) const
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
