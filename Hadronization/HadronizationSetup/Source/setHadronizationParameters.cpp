#include <iostream>
#include "../Include/setHadronizationParameters.h"
#include "../Include/runSetupInterface.h"
#include "../../../Utilities/Include/Constants.h"
#include "../../../Utilities/Include/HadronizationException.h"

void setHadronizationParameters(HadronizationSetup &io_setupParameters)
{
    
    // Set fit parameters
	io_setupParameters.energyDensity = 0.35;
	io_setupParameters.gammaS = 0.65;	

    // Set partition function interpolation parameters
    // TODO: to be used together with the final handling of partition
    // function data structure and likely moved to another place
	io_setupParameters.minClusterMass = 0.26996;
	io_setupParameters.maxClusterMass = 10.;
	io_setupParameters.minEnergyDensity = 0.2;
	io_setupParameters.maxEnergyDensity = 0.5;
	io_setupParameters.minGammaS = 0.5;	
	io_setupParameters.maxGammaS = 1.;	
	
    // Set channel sampling parameters
	io_setupParameters.samplingTemperature = 0.160;
	io_setupParameters.samplingEnergyDensity = io_setupParameters.minEnergyDensity;
	
    // Set cluster merging parameters
	io_setupParameters.clusterMergingFlag = 1;
	io_setupParameters.clusterMergingMinMass = 1.8;
    io_setupParameters.charmClusterMergingMinMass = 3.5;
    io_setupParameters.bottomClusterMergingMinMass = 5.7;

    // TODO: add check about heavy and global cluster mass
    // limits! charmClusterMergingMinMass>clusterMergingFlag and
    // bottomClusterMergingMinMass>clusterMergingFlag otherwise
    // the cluster merging procedure could fail (code has not been
    // tested in case of limits not ordered as above)
        
    // Set random number generator seed and initialize the generation algorithm
	io_setupParameters.randomNumberGeneratorSeed = 9019984;
	
    // Set partition function grid path 
	io_setupParameters.partitionFunctionDataSetPath =
        "/Users/Christopher/Development/WorkingDirs/MCSTHAR++CPC/MCSTHAR++/Resources/Herwig6510/PartitionFunctionSet";
    
    // Set hadron list file
	io_setupParameters.hadronDataSetFileName = "../../Resources/Herwig6510/HadronList.ascii";
    
    // Set light flavored hadrons maximum mass value
    io_setupParameters.lightHadronMaxMass = 1.8;
	
    // Run parameter setup user interface
	runSetupInterface(io_setupParameters);
}
