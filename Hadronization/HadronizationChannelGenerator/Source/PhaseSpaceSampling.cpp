#include <cmath>
#include "../Include/PhaseSpaceSampling.h"
#include "../../../Utilities/Include/Constants.h"
#include "../../../Utilities/Include/HadronizationException.h"

PhaseSpaceSampling::PhaseSpaceSampling(void)
                                      :m_isPhaseSpaceAvailable(false)
                                      ,m_weight(0.)
                                      ,m_randomGenerator(NULL)
                                      ,m_jacobianFactor(0.)
{
}

PhaseSpaceSampling::PhaseSpaceSampling(RandomNumberGenerator& i_randomGenerator)
                                      :m_isPhaseSpaceAvailable(false)
                                      ,m_weight(0.)
                                      ,m_randomGenerator(&i_randomGenerator)
                                      ,m_jacobianFactor(0.)
{
    if(m_randomGenerator==NULL)
    {
        throw HadronizationException("Error during phase space sampling object creation, invalid random number generator provided",
                                     __FUNCTION__,128);
    }
}

PhaseSpaceSampling::~PhaseSpaceSampling(void)
{
}

void PhaseSpaceSampling::setRandomNumberGenerator(RandomNumberGenerator& io_randomGenerator)
{
    m_randomGenerator = &io_randomGenerator;
    if(m_randomGenerator==NULL)
    {
        throw HadronizationException("Error during phase space sampling object creation, invalid random number generator provided",
                                     __FUNCTION__,128);
    }
}


const vector<TLorentzVector>& PhaseSpaceSampling::getPhaseSpace(void) const
{
    if(m_isPhaseSpaceAvailable)
    {
        return m_hadronMomenta;
    }
    else
    {
        throw HadronizationException("Error during phase space configuration retrieval, configuration not available",
                                     __FUNCTION__,121);
    }
}

double PhaseSpaceSampling::getSamplingWeight(void) const
{
    if(m_isPhaseSpaceAvailable)
    {
        return m_weight;
    }
    else
    {
        throw HadronizationException("Error during phase space configuration sampling weight retrieval, configuration not available",
                                     __FUNCTION__,122);
    }
}

void PhaseSpaceSampling::run(const TLorentzVector& i_clusterMomentum, const vector<double>& i_hadronMasses)
{
    if(m_randomGenerator)
    {        
        m_isPhaseSpaceAvailable = false;
        
        // Check number of input hadron masses
        const unsigned int numberOfHadrons(i_hadronMasses.size());
        if(numberOfHadrons<2)
        {
            throw HadronizationException("Error during phase space configuration sampling, less than 2 hadron masses provided",
                                         __FUNCTION__,123);
        }

        // Check energy availability for the provided input
        double hadronMassSum = i_hadronMasses[0];
        for(unsigned int hadronIndex=1;hadronIndex<numberOfHadrons;++hadronIndex)
        {
            hadronMassSum += i_hadronMasses[hadronIndex];
        }
        
        // TODO: refine check, add zero kinetic energy case
        if(i_clusterMomentum.Mag()<=hadronMassSum)
        {
            throw HadronizationException("Error during phase space configuration sampling, hadron mas sum is larger than cluster one",
                                         __FUNCTION__,124);
        }
        
        // Random residual squared mass for 2-body split
        double randomSquaredMass;
        // Angular variables of the sampled hadron momentum for each two-body splitting
        double cosTheta;
        double phi;
        
        // Reset sampling framework
        m_hadronMomenta.resize(numberOfHadrons);
        m_weight = 1.;
        m_residualMomentum = i_clusterMomentum;
        m_jacobianFactor = 1.;
        
        // Phase space configuration generation splitting loop
        for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons-2;++hadronIndex)
        {
            
            // Generate random squared residual mass
            try
            {
                randomSquaredMass = generateSquaredResidualMass(hadronIndex,i_hadronMasses);
            }
            catch(HadronizationException& exception)
            {
                throw exception;
            }
            
            // Generate angular variables for the sampled hadron momentum
            // TODO: perform generation of the full set of needed random numbers in one call
            cosTheta = m_randomGenerator->getUniformRandom(-1.,1.);
            phi = m_randomGenerator->getUniformRandom(0.,twoPi);
                    
            try
            {
                // Perform two-body split
                m_hadronMomenta[hadronIndex] = performTwoBodySplit(i_hadronMasses[hadronIndex],randomSquaredMass,cosTheta,phi);
                
            }
            catch(HadronizationException& exception)
            {
                throw exception;
            }
                    
            // Update total available momentum
            m_residualMomentum -= m_hadronMomenta[hadronIndex];
                
        }

        // Perform last split
        // Generate angular variables for the sampled hadron momentum
        cosTheta = m_randomGenerator->getUniformRandom(-1.,1.);
        phi = m_randomGenerator->getUniformRandom(0.,twoPi);
        
        try
        {
            // Perform two-body split
            m_hadronMomenta[numberOfHadrons-2] =
                performTwoBodySplit(i_hadronMasses[numberOfHadrons-2],
                                    i_hadronMasses[numberOfHadrons-1]*i_hadronMasses[numberOfHadrons-1],cosTheta,phi);
        }
        catch(HadronizationException& exception)
        {
            throw exception;
        }
        m_hadronMomenta[numberOfHadrons-1] = m_residualMomentum - m_hadronMomenta[numberOfHadrons-2];
                
        // Update phase space sampling weight
        m_weight *= 2.*pow(twoPi,static_cast<int>(numberOfHadrons-1))*m_jacobianFactor;
        
        // Update phase space sampling weight with the hadron energy contribution:
        // the hadron energy must be computed in the hadronizing cluster rest frame
        // TODO: avoid this step, hadron energy in the right rf are already available
        const TVector3 boostVector(-i_clusterMomentum.BoostVector());
        TLorentzVector hadronMomentum;
        for(unsigned int hadronIndex=0;hadronIndex<numberOfHadrons;++hadronIndex)
        {
            hadronMomentum = m_hadronMomenta[hadronIndex];
            hadronMomentum.Boost(boostVector);
            m_weight *= hadronMomentum.E();
        }
        
        m_isPhaseSpaceAvailable = true;
        
    }
    else
    {
        throw HadronizationException("Error during cluster hadronization, random number generator not available",
                                     __FUNCTION__,127);
    }
}

void PhaseSpaceSampling::run(const Cluster& i_cluster, const vector<HadronData>& i_hadronData)
{
    // Retrieve hadron masses
    const unsigned int numberOfHadrons(i_hadronData.size());
    vector<double> hadronMasses(numberOfHadrons,0.);
    for(unsigned int hadronIndex = 0;hadronIndex<numberOfHadrons;++hadronIndex)
    {
        hadronMasses[hadronIndex] = i_hadronData[hadronIndex].getMass();
    }
    
    try
    {
        // Run phase space sampling
        run(i_cluster.getP(),hadronMasses);
    }
    catch(HadronizationException& exception)
    {
        throw exception;
    }
}


TLorentzVector PhaseSpaceSampling::performTwoBodySplit(const double i_hadronMass,const double i_residualSquaredMass,
                                                       const double i_cosTheta,const double i_phi)
{
    // TODO: optimize calculation
    const double startSquaredMass(m_residualMomentum.Mag2());
    
    // Compute the 3-momentum module and sin(theta) for the back-to-back two-body split
	double hadronMomentumModule =
        (1./(4.*startSquaredMass))*(startSquaredMass - pow(i_hadronMass + sqrt(i_residualSquaredMass),2))*
        (startSquaredMass - pow(i_hadronMass - sqrt(i_residualSquaredMass),2));
	   
    if(hadronMomentumModule>=0.)
    {
		hadronMomentumModule = sqrt(hadronMomentumModule);
	}
    else
    {
        // TODO: veirify and refine this check
        if(-hadronMomentumModule<phaseSpace3MomentumModTolerance)
        {
            hadronMomentumModule = 0.;
        }
        else
        {
            // The error is too large, throw exception
            throw HadronizationException("Error during phase space configuration sampling, negative 3-momentum squared modulus generated",
                                         __FUNCTION__,126);
        }
    }
	const double sinTheta(sqrt(1. - i_cosTheta*i_cosTheta));

    // Sampled momentum corresponding to the above 3-momentum
	TLorentzVector o_hadronMomentum(hadronMomentumModule*sinTheta*cos(i_phi),
                                    hadronMomentumModule*sinTheta*sin(i_phi),
                                    hadronMomentumModule*i_cosTheta,
                                    sqrt(i_hadronMass*i_hadronMass + hadronMomentumModule*hadronMomentumModule));
    
    // TODO: add check to avoid division by zero
    // Compute the phase space integral jacobian contribution corresponding to the sampled momentum
	m_jacobianFactor = m_jacobianFactor*hadronMomentumModule/
        (o_hadronMomentum.E()+sqrt(hadronMomentumModule*hadronMomentumModule+i_residualSquaredMass));
    
    // Boost the generated momentum from the two-body back-to-back reference frame to the original collision one
	const TVector3 boostVector(m_residualMomentum.BoostVector());
	o_hadronMomentum.Boost(boostVector);
    
    // Return generated momentum
    return o_hadronMomentum;
}

double PhaseSpaceSampling::generateSquaredResidualMass(const unsigned int i_hadronIndex,
                                                       const vector<double>& i_hadronMasses)
{    
    // Compute minimum squared residual invariant mass
    double minSquaredMass = i_hadronMasses[i_hadronIndex+1];
    for(unsigned int residualHadronIndex=i_hadronIndex+2;residualHadronIndex<i_hadronMasses.size();
        ++residualHadronIndex)
    {
        // TODO: avoid this sum at every loop iteration
        minSquaredMass += i_hadronMasses[residualHadronIndex];
    }
    minSquaredMass *= minSquaredMass;
    
    // Compute maximum squared residual invariant mass
    double maxSquaredMass = (m_residualMomentum.Mag() - i_hadronMasses[i_hadronIndex]);
    maxSquaredMass *= maxSquaredMass;
    
    // Residual invariant mass limit check
    const double diffSquaredMass(maxSquaredMass - minSquaredMass);
    double o_randomSquaredMass;
    if(diffSquaredMass<=0.)
    {
        if(-diffSquaredMass<phaseSpaceSquaredMassTolerance)
        {
            // The error is within defined tolerance, continue the sampling procedure
            o_randomSquaredMass = maxSquaredMass;
        }
        else
        {
            // The error is too large, throw exception
            throw HadronizationException("Error during phase space configuration sampling, invariant mass limits violation",
                                         __FUNCTION__,125);
        }
    }
    else
    {
        // Generate residual squared invariant mass
        // TODO: perform generation of the full set of needed random numbers in one call
        o_randomSquaredMass = m_randomGenerator->getUniformRandom(minSquaredMass,maxSquaredMass);
    }
    
    // Update phase space sampling weight
    // TODO: do we need abs??? We have to set it to zero instead!
    m_weight *= std::abs(diffSquaredMass);
    
    // Return sampled squared mass
    return o_randomSquaredMass;
}
