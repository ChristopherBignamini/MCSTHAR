#include "../Include/MCSTHARSetup.h"
#include "../Include/HadronizationSetupLoader.h"
#include "../../../Utilities/Include/HadronizationException.h"

HadronizationSetupLoader::HadronizationSetupLoader(const string& i_hadronizationSetupFileName)
                                                  :m_hadronizationSetupAvailability(false)
{    
    unique_ptr<::MCSTHARSetup> hadronizationSetupFile;
    try
    {
        // Parse setup file
        hadronizationSetupFile = MCSTHARSetup_(i_hadronizationSetupFileName,
                                               xml_schema::flags::dont_validate);
    }
    catch(const xml_schema::exception& exception)
    {
        // Error during setup file parsing, throw exception
        string errorMessage("Unable to parse the specified setup file. ");
        errorMessage += exception.what();
        throw HadronizationException(errorMessage,__FUNCTION__,401);
    }
    
    // Check if provided string is empty
    m_hadronizationSetup.partitionFunctionDataSetPath = hadronizationSetupFile->partitionFunctionDataSetPath();
    if(m_hadronizationSetup.partitionFunctionDataSetPath.empty())
    {
        // Empty data set parh provided, throw exception
        throw HadronizationException("Error during hadronization setup, empty partition function data set path provided",
                                     __FUNCTION__,402);
    }
    
    // Check if provided string is empty
    m_hadronizationSetup.hadronDataSetFileName = hadronizationSetupFile->hadronDataFileName();
    if(m_hadronizationSetup.hadronDataSetFileName.empty())
    {
        // Empty data file name provided, throw exception
        throw HadronizationException("Error during hadronization setup, empty hadron data file name provided",
                                     __FUNCTION__,403);
    }
    
    // Check energy density value
    m_hadronizationSetup.energyDensity = hadronizationSetupFile->clusterEnergyDensity();
    if(m_hadronizationSetup.energyDensity<=0.)
    {
        // Non positive energy density provided, throw exception
        throw HadronizationException("Error during hadronization setup, non positive energy density parameter provided",
                                     __FUNCTION__,404);
    }
    
    // Check strangeness suppression parameter value
    m_hadronizationSetup.gammaS = hadronizationSetupFile->strangenessSuppressionParameter();
    if((m_hadronizationSetup.gammaS<0.) || (m_hadronizationSetup.gammaS>1.))
    {
        // Strangeness suppression parameter outside [0,1] range provided, throw exception
        throw HadronizationException("Error during hadronization setup, strangeness suppression parameter outside [0,1] range provided",
                                     __FUNCTION__,405);
    }
    
    m_hadronizationSetup.clusterMergingFlag = false;
    // Retrieve cluster merging (light and heavy flavored) cluster minimum mass
    if(hadronizationSetupFile->clusterMergingMinimumMass().present())
    {
        m_hadronizationSetup.clusterMergingMinMass = hadronizationSetupFile->clusterMergingMinimumMass().get();
        
        if(m_hadronizationSetup.clusterMergingMinMass<0.)
        {
            // Negative minimum mass provided, throw exception
            throw HadronizationException("Error during hadronization setup, negative minimum cluster merging (light and heavy flavored) mass provided",
                                         __FUNCTION__,406);
        }
        
        if(m_hadronizationSetup.clusterMergingMinMass>m_hadronizationSetup.maxClusterMass)
        {
            // Minimum mass outside available range provided, throw exception
            throw HadronizationException("Error during hadronization setup, cluster merging (light and heavy flavored) mass outside available range provided",
                                         __FUNCTION__,406);
        }
        
        m_hadronizationSetup.clusterMergingFlag = true;
    }
    else
    {
        m_hadronizationSetup.clusterMergingMinMass = -1.;
    }

    // Retrieve chamr cluster merging minimum mass
    if(hadronizationSetupFile->charmClusterMergingMinimumMass().present())
    {
        m_hadronizationSetup.charmClusterMergingMinMass = hadronizationSetupFile->charmClusterMergingMinimumMass().get();

        if(m_hadronizationSetup.charmClusterMergingMinMass<0.)
        {
            // Negative minimum mass provided, throw exception
            throw HadronizationException("Error during hadronization setup, negative minimum charm cluster merging mass provided",
                                         __FUNCTION__,407);
        }
        
        if(m_hadronizationSetup.charmClusterMergingMinMass>m_hadronizationSetup.maxClusterMass)
        {
            // Minimum mass outside available range provided, throw exception
            throw HadronizationException("Error during hadronization setup, charm cluster merging (mass outside available range provided",
                                         __FUNCTION__,407);
        }
        
        if(m_hadronizationSetup.clusterMergingFlag==false)
        {
            m_hadronizationSetup.clusterMergingFlag = true;
        }
    }
    

    // Retrieve bottom cluster merging minimum mass
    if(hadronizationSetupFile->bottomClusterMergingMinimumMass().present())
    {
        m_hadronizationSetup.bottomClusterMergingMinMass = hadronizationSetupFile->bottomClusterMergingMinimumMass().get();
        
        if(m_hadronizationSetup.bottomClusterMergingMinMass<0.)
        {
            // Negative minimum mass provided, throw exception
            throw HadronizationException("Error during hadronization setup, negative minimum bottom cluster merging mass provided",
                                         __FUNCTION__,408);
        }
        
        if(m_hadronizationSetup.bottomClusterMergingMinMass>m_hadronizationSetup.maxClusterMass)
        {
            // Minimum mass outside available range provided, throw exception
            throw HadronizationException("Error during hadronization setup, bottom cluster merging (mass outside available range provided",
                                         __FUNCTION__,408);
        }
        if(m_hadronizationSetup.clusterMergingFlag==false)
        {
            m_hadronizationSetup.clusterMergingFlag = true;
        }

    }

    // Retrieve random number generator seed
    if(hadronizationSetupFile->randomNumberGeneratorSeed().present())
    {
        m_hadronizationSetup.randomNumberGeneratorSeed = hadronizationSetupFile->randomNumberGeneratorSeed().get();
    }
    
    m_hadronizationSetupAvailability = true;
}

HadronizationSetupLoader::~HadronizationSetupLoader(void)
{
}

const HadronizationSetup& HadronizationSetupLoader::getHadronizationSetup(void)
{
    if(!m_hadronizationSetupAvailability)
    {
        // Hadronization parameters not correctly set, throw exception
        throw HadronizationException("Error during hadronization parameter retrieveal, parameters not correctly set",
                                     __FUNCTION__,402);
    }
    else
    {
        return m_hadronizationSetup;
    }
}
