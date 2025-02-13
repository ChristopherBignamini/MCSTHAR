############### General compiler settings ###############

CPP = /opt/homebrew/opt/llvm/bin/clang++
ifeq ($DEBUG,0)
#ifnotdebug
    CPPFLAGS = -O3 -Wall -Wreorder -std=c++17
    FFLAGS = -O3 -Wall
else
    CPPFLAGS = -g -Wall -Wreorder -std=c++17
    FFLAGS = -g -Wall
endif
#F = gfortran
F = /opt/homebrew/bin/gfortran
#else
#endif
 
#########################################################

############### Compiler external include and link flags#
# TODO: switch to autotools..
#ROOT_INCLUDE = $(shell root-config --cflags)
ROOT_INCLUDE = -D_REENTRANT -pthread -m64 -I/opt/homebrew/include/root
#ROOT_LINK = $(shell root-config --libs)
ROOT_LINK = -L/opt/homebrew/lib/root -lCore -lRIO -lNet -lHist \
	-lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl
XSD_INCLUDE = -I/opt/homebrew/include
XSD_LINK = -L/opt/homebrew/lib/ -lxerces-c
XERCES_INCLUDE = -I/opt/homebrew/include/xercesc/
HERWIG_INCLUDE = -I/Users/bignamic/Development/HEP/MCSTHAR/Main/Herwig6510/Include
HEPMC_INCLUDE = -I/opt/homebrew/Cellar/hepmc2/2.06.11/include/
HEPMC_LINK = -L/opt/homebrew/Cellar/hepmc2/2.06.11/lib/ -lHepMC -lHepMCfio
GSL_INCLUDE = -I/opt/homebrew/include/
GSL_LIB = -L/opt/homebrew/lib/ -lgsl
#########################################################

############### Compiler fortran libraries link flag ####

#FLINKFLAG = -L/usr/local/lib/gcc/x86_64-apple-darwin10.4.0/4.6.0 \
        -lgfortran -lgcc_ext.10.5 -lgcc -no_compact_unwind -lSystem
FLINKFLAG = -L/opt/homebrew/lib/gcc/current/ -lgfortran 

#########################################################

############### Targets and requirements ################

VPATH = ../../Hadronization/HadronizationChannelGenerator/Include:\
	../../Hadronization/HadronizationChannelGenerator/Source:\
	../../Hadronization/HadronizationHadronSet/Include:\
	../../Hadronization/HadronizationHadronSet/Source:\
	../../Hadronization/HadronizationObjects/Include:\
	../../Hadronization/HadronizationObjects/Source:\
	../../Hadronization/HadronizationSetup/Include:\
	../../Hadronization/HadronizationSetup/Source:\
	../../Hadronization/HadronizationSteps/Include:\
	../../Hadronization/HadronizationSteps/Source:\
	../../Hadronization/HadronizationPartitionFunctionStorage/Include:\
	../../Hadronization/HadronizationPartitionFunctionStorage/Source:\
	../../Interfaces/Herwig6510/Include:\
	../../Interfaces/Herwig6510/Source:\
	../../Utilities/Include:\
	../../Utilities/Source:\
	../../Analysis/HepMCEventStorage/Include:\
	../../Analysis/HepMCEventStorage/Source:\
	Source:\
	Include:\

MCSTHAROBJECTS = \
	HadronizationObject.o \
	Particle.o \
	Cluster.o \
	ClusterMerging.o \
	HadronizationHandler.o \
	HadronSet.o \
	HadronData.o \
	HadronDataSet.o \
	HadronDataSetLoader.o \
	HadronSamplingGroups.o \
	computeWaveFunctionStrangeComponent.o \
	PhaseSpaceSampling.o \
	HadronSampling.o \
	ClusterHadronization.o \
	HadronizationEventRecord.o \
	HadronizationException.o \
	RandomNumberGenerator.o \
	GSLRandomNumberGenerator.o \
	ROOTRandomNumberGenerator.o \
	setHadronizationParameters.o \
	runSetupInterface.o \
	setParameters.o \
	MCSTHARSetup.o \
	HadronizationSetupLoader.o \
	PartitionFunctionHandling.o \
	PartitionFunctionDataFile.o \
	PartitionFunctionDataSummaryFile.o \
	PartitionFunctionArchiveFile.o \
	PartitionFunctionArchiveBaseTypes.o \
	Vector3D.o \
	LorentzVector.o \
	fileUtils.o 

INTERFACEOBJECTS = \
	hwcforMCSTHAR.o \
	loadClusters.o \
	convertParticleId.o \
	updateEventRecord.o

# ANALYSISOBJECTS = \
# 	buildHepMCEventRecord.o

HERWIGOBJECTS = herwig6510.o

OBJECTS = $(MCSTHAROBJECTS) $(INTERFACEOBJECTS) $(HERWIGOBJECTS) $(ANALYSISOBJECTS) herwig.o

#########################################################

############### Compile main and build executable #######

EXE = herwig.exe

default: $(EXE)

herwig.exe: $(OBJECTS)
	$(CPP) -o $(EXE) $(OBJECTS) $(FLINKFLAG) $(ROOT_LINK) $(GSL_LIB) $(XSD_LINK) $(HEPMC_LINK) $(EFENCE_LINK) 

#buildHepMCEventRecord.h 
#$(HEPMC_INCLUDE)
herwig.o: herwig.cpp \
	hwcforMCSTHAR.h loadClusters.h \
	updateEventRecord.h Herwig6510Wrapper.h \
	HadronizationSetupLoader.h \
	HadronizationHandler.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $(ROOT_INCLUDE) $<

############### Compile MCSTHAR++ code ##################

HadronizationObject.o: HadronizationObject.cpp HadronizationObject.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $< 

Particle.o: Particle.cpp Particle.h HadronizationObject.h HadronData.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $< 

Cluster.o: Cluster.cpp Cluster.h HadronizationObject.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

ClusterMerging.o: ClusterMerging.cpp ClusterMerging.h Cluster.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

HadronizationHandler.o: HadronizationHandler.cpp HadronizationHandler.h \
	HadronizationEventRecord.h HadronizationException.h \
	ClusterMerging.h HadronizationSetup.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $(ROOT_INCLUDE) $<

HadronSet.o: HadronSet.cpp HadronSet.h HadronizationSetup.h HadronData.h
	$(CPP) -c $(CPPFLAGS) $<

HadronData.o: HadronData.cpp HadronData.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $<

HadronDataSet.o: HadronDataSet.cpp HadronDataSet.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

HadronDataSetLoader.o: HadronDataSetLoader.cpp HadronDataSetLoader.h HadronDataSet.h \
	HadronData.h HadronizationException.h 
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

HadronSamplingGroups.o: HadronSamplingGroups.cpp HadronSamplingGroups.h HadronSet.h
	$(CPP) -c $(CPPFLAGS) $<

computeWaveFunctionStrangeComponent.o: \
	computeWaveFunctionStrangeComponent.cpp computeWaveFunctionStrangeComponent.h
	$(CPP) -c $(CPPFLAGS) $<

PhaseSpaceSampling.o: PhaseSpaceSampling.cpp \
	PhaseSpaceSampling.h Cluster.h HadronData.h Constants.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

HadronSampling.o: HadronSampling.cpp \
	HadronSampling.h ChargeConfiguration.h HadronData.h Constants.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(GSL_INCLUDE) $<

ClusterHadronization.o: ClusterHadronization.cpp ClusterHadronization.h HadronSamplingGroups.h\
	PhaseSpaceSampling.h Cluster.h HadronData.h HadronizationChannel.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $(ROOT_INCLUDE) $<

HadronizationEventRecord.o: HadronizationEventRecord.cpp HadronizationEventRecord.h \
	HadronizationChannel.h HadronizationObject.h Cluster.h Particle.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

HadronizationException.o: HadronizationException.cpp HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $<

RandomNumberGenerator.o: RandomNumberGenerator.cpp RandomNumberGenerator.h
	$(CPP) -c $(CPPFLAGS) $<

GSLRandomNumberGenerator.o: GSLRandomNumberGenerator.cpp GSLRandomNumberGenerator.h \
	RandomNumberGenerator.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(GSL_INCLUDE) $<

ROOTRandomNumberGenerator.o: ROOTRandomNumberGenerator.cpp ROOTRandomNumberGenerator.h \
	ROOTRandomNumberGenerator.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

setHadronizationParameters.o: \
	setHadronizationParameters.cpp setHadronizationParameters.h \
	runSetupInterface.h Constants.h \
	HadronizationSetup.h HadronSet.h \
	computeWaveFunctionStrangeComponent.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

runSetupInterface.o: runSetupInterface.cpp runSetupInterface.h setParameters.h HadronizationSetup.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

setParameters.o: setParameters.cpp setParameters.h HadronizationSetup.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

MCSTHARSetup.o: MCSTHARSetup.cpp MCSTHARSetup.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

HadronizationSetupLoader.o: HadronizationSetupLoader.cpp HadronizationSetupLoader.h \
	MCSTHARSetup.h HadronizationSetup.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

PartitionFunctionHandling.o: PartitionFunctionHandling.cpp PartitionFunctionHandling.h \
	PartitionFunctionDataSummaryFile.h PartitionFunctionDataSummaryFile.h \
	PartitionFunctionDataFile.h partitionFunctionDataConstants.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

PartitionFunctionDataFile.o: PartitionFunctionDataFile.cpp PartitionFunctionDataFile.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

PartitionFunctionDataSummaryFile.o: PartitionFunctionDataSummaryFile.cpp PartitionFunctionDataSummaryFile.h \
	PartitionFunctionArchiveFile.h 
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

PartitionFunctionArchiveFile.o: PartitionFunctionArchiveFile.cpp PartitionFunctionArchiveFile.h \
	PartitionFunctionArchiveBaseTypes.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

PartitionFunctionArchiveBaseTypes.o: PartitionFunctionArchiveBaseTypes.cpp PartitionFunctionArchiveBaseTypes.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $<

Vector3D.o: Vector3D.cpp Vector3D.h
	$(CPP) -c $(CPPFLAGS) $<

LorentzVector.o: LorentzVector.cpp LorentzVector.h \
	Vector3D.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $<

fileUtils.o: fileUtils.cpp fileUtils.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $<

#partfunct.o: partfunct.cpp partfunct.h HadronSamplingGroupHandler.h \
#	HadronizationChannelGenerator.h kingen_data.h Particle.h \
#	PhaseSpaceIntegrationWrapper.h
#	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $<

#omega.o: omega.f 
#	$(F) -c $(FFLAGS) $<

#########################################################

############### Compile Analysis code ############

# buildHepMCEventRecord.o: buildHepMCEventRecord.cpp buildHepMCEventRecord.h
# 	$(CPP) -c $(CPPFLAGS) $(HEPMC_INCLUDE) $< 

#########################################################

############### Compile Herwig specific code ############

herwig6510.o: herwig6510.f herwig6510.inc HERWIG65.INC
	$(F) -c $(FFLAGS) $(HERWIG_INCLUDE) $< 

#########################################################

############### Compile Herwig specific interface code ############

hwcforMCSTHAR.o: \
	hwcforMCSTHAR.f HERWIG65.INC
	$(F) -c $(FFLAGS) $(HERWIG_INCLUDE) $(EFENCE_LINK) $< 

loadClusters.o: loadClusters.cpp loadClusters.h \
	Herwig6510Constants.h Herwig6510Wrapper.h \
	HadronizationException.h ChargeData.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(HEPMC_INCLUDE) $(EFENCE_LINK) $<

convertParticleId.o: \
	convertParticleId.f herwig6510.inc HERWIG65.INC
	$(F) -c $(FFLAGS) $(HERWIG_INCLUDE) $(EFENCE_LINK) $<

updateEventRecord.o: updateEventRecord.cpp updateEventRecord.h \
	HadronizationEventRecord.h Herwig6510Constants.h \
	convertParticleId.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(EFENCE_LINK) $<

#########################################################

############### Clean compilation file ##################

cleanall:
	rm -rf *.exe $(OBJECTS)

cleanherwig:
	rm -rf $(HERWIGOBJECTS)

cleanmcsthar:
	rm -rf $(MCSTHAROBJECTS) herwig.o herwig.exe

cleaninterface:
	rm -rf $(INTERFACEOBJECTS)

cleananalysis:
	rm -rf $(ANALYSISOBJECTS)


#########################################################
