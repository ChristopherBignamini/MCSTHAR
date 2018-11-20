############### General compiler settings ###############
CPP = mpic++
CPPFLAGS = -O3 -Wall -Wreorder
#CPPFLAGS = -g -Wall -Wreorder

#########################################################

############### Compiler external include and link flags#
# TODO: switch to autotools..
#ROOT_INCLUDE = $(shell root-config --cflags)
ROOT_INCLUDE = -D_REENTRANT -pthread -I/Users/Christopher/Development/Libraries/root_5_28//include/root
#ROOT_LINK = $(shell root-config --libs)
ROOT_LINK = -L/Users/Christopher/Development/Libraries/root_5_28//lib/root -lCore -lCint -lRIO -lNet -lHist \
	-lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl
XSD_INCLUDE = -I/Users/Christopher/Development/Libraries/xsd-3.3.0-i686-macosx/libxsd
XSD_LINK = -L/Users/Christopher/Development/Libraries/xerces-c_3_1_1/lib -lxerces-c
XERCES_INCLUDE = -I/Users/Christopher/Development/Libraries/xerces-c_3_1_1/include
OMP_LINK = -lgomp
GSL_INCLUDE = -I/Users/Christopher/Development/Libraries/gsl_1_14/include 
GSL_LIB = -L/Users/Christopher/Development/Libraries/gsl_1_14/lib -lgsl
#RANDOMGEN = -D_USEROOTRANDOM 
RANDOMGEN = 
LOGLEV = -D_FULLLOG
#LOGLEV = 
#########################################################

############### Targets and requirements ################

VPATH = ../../Hadronization/HadronizationChannelGenerator/Include:\
	../../Hadronization/HadronizationChannelGenerator/Source:\
	../../Hadronization/HadronizationHadronSet/Include:\
	../../Hadronization/HadronizationHadronSet/Source:\
	../../Hadronization/HadronizationPartitionFunction/Include:\
	../../Hadronization/HadronizationPartitionFunction/Source:\
	../../Hadronization/HadronizationPartitionFunctionStorage/Include:\
	../../Hadronization/HadronizationPartitionFunctionStorage/Source:\
	../../Hadronization/HadronizationSetup/Include:\
	../../Hadronization/HadronizationSetup/Source:\
	../../Utilities/Include:\
	../../Utilities/Source:\
	Source:\
	Include:\

MCSTHAROBJECTS = \
	HadronSet.o \
	HadronData.o \
	HadronDataSet.o \
	HadronDataSetLoader.o \
	HadronSamplingGroups.o \
	HadronSampling.o \
	HadronSamplingNew.o \
	PartitionFunctionGenerator.o \
	PartitionFunctionGeneratorMaster.o \
	PartitionFunctionGeneratorSlave.o \
	PartitionFunctionGenerationHandler.o \
	PartitionFunctionGenerationHandlerMaster.o \
	PartitionFunctionGenerationHandlerSlave.o \
	PartitionFunctionArchiveWriter.o \
	PartitionFunctionDataFile.o \
	PartitionFunctionDataSummaryFile.o \
	PartitionFunctionArchiveFile.o \
	PartitionFunctionArchiveBaseTypes.o \
	MCSTHARPartitionFunctionCalculationSetup.o \
	PartitionFunctionCalculationSetupLoader.o \
	PhaseSpaceIntegrator.o \
	PartitionFunctionPhaseSpaceSampling.o \
	PhaseSpaceIntegrandFunction.o \
	PhaseSpaceExactIntegrandFunction.o \
	PhaseSpaceApproximatedIntegrandFunction.o \
	HadronizationException.o \
	RandomNumberGenerator.o \
	GSLRandomNumberGenerator.o \
	ROOTRandomNumberGenerator.o \
	logMessage.o \
	computeWaveFunctionStrangeComponent.o \
	fileUtils.o \
	buildParameterRegularGrid.o \
	mathFunctions.o


OBJECTS = $(MCSTHAROBJECTS) partitionFunctionSetBuilder.o

#########################################################

############### Compile main and build executable #######

EXE = partitionFunctionSetBuilder

default: $(EXE)

partitionFunctionSetBuilder: $(OBJECTS)
	$(CPP) $(CPPFLAGS) -o $(EXE) $(OBJECTS) $(OMP_LINK) $(ROOT_LINK) $(XSD_LINK) $(GSL_LIB) $(RANDOMGEN)

partitionFunctionSetBuilder.o: partitionFunctionSetBuilder.cpp logMessage.h \
	PartitionFunctionCalculationSetupLoader.h PartitionFunctionGenerationHandlerMaster.h \
	PartitionFunctionGenerationHandlerSlave.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $(RANDOMGEN) $(ROOT_INCLUDE) $< 

############### Compile MCSTHAR++ code ##################

MCSTHARPartitionFunctionCalculationSetup.o: MCSTHARPartitionFunctionCalculationSetup.cpp \
	MCSTHARPartitionFunctionCalculationSetup.h PartitionFunctionArchiveBaseTypes.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $< 

PartitionFunctionCalculationSetupLoader.o: PartitionFunctionCalculationSetupLoader.cpp \
	PartitionFunctionCalculationSetupLoader.h HadronizationException.h \
	MCSTHARPartitionFunctionCalculationSetup.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $< 

PartitionFunctionGenerationHandler.o: PartitionFunctionGenerationHandler.cpp PartitionFunctionGenerationHandler.h \
	PartitionFunctionGenerator.h HadronSamplingGroups.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $< 

PartitionFunctionGenerationHandlerMaster.o: PartitionFunctionGenerationHandlerMaster.cpp PartitionFunctionGenerationHandlerMaster.h \
	PartitionFunctionGenerationHandler.h PartitionFunctionArchiveWriter.h PartitionFunctionCalculationSetup.h \
	PartitionFunctionGeneratorMaster.h logMessage.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(XSD_INCLUDE) $(XERCES_INCLUDE) $(GSL_INCLUDE) $< 

PartitionFunctionGenerationHandlerSlave.o: PartitionFunctionGenerationHandlerSlave.cpp PartitionFunctionGenerationHandlerSlave.h \
	PartitionFunctionGenerationHandler.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(GSL_INCLUDE) $< 

buildParameterRegularGrid.o: buildParameterRegularGrid.cpp buildParameterRegularGrid.h \
	HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $< 

PartitionFunctionGenerator.o: PartitionFunctionGenerator.cpp PartitionFunctionGenerator.h \
	ChargeConfiguration.h HadronSamplingGroups.h HadronSampling.h HadronSamplingNew.h \
	HadronizationException.h RandomNumberGenerator.h PartitionFunctionPhaseSpaceSampling.h buildParameterRegularGrid.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(LOGLEV) $(GSL_INCLUDE) $(RANDOMGEN) $< 

PartitionFunctionGeneratorMaster.o: PartitionFunctionGeneratorMaster.cpp PartitionFunctionGeneratorMaster.h \
	PartitionFunctionGenerator.h PartitionFunctionGenerationData.h logMessage.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(LOGLEV) $(GSL_INCLUDE) $< 

PartitionFunctionGeneratorSlave.o: PartitionFunctionGeneratorSlave.cpp PartitionFunctionGeneratorSlave.h \
	PartitionFunctionGenerator.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(LOGLEV) $(GSL_INCLUDE) $< 

PartitionFunctionArchiveWriter.o: PartitionFunctionArchiveWriter.cpp PartitionFunctionArchiveWriter.h \
	PartitionFunctionDataFile.h PartitionFunctionArchiveFile.h \
	partitionFunctionDataConstants.h HadronizationException.h
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

PhaseSpaceIntegrator.o: PhaseSpaceIntegrator.cpp PhaseSpaceIntegrator.h \
	PhaseSpaceIntegrationData.h PhaseSpaceIntegrandFunction.h \
	PhaseSpaceExactIntegrandFunction.h PhaseSpaceApproximatedIntegrandFunction.h \
	HadronizationException.h logMessage.h GSLRandomNumberGenerator.h ROOTRandomNumberGenerator.h \
	RandomNumberGenerator.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(GSL_INCLUDE) $(RANDOMGEN) $< 

PartitionFunctionPhaseSpaceSampling.o: PartitionFunctionPhaseSpaceSampling.cpp PartitionFunctionPhaseSpaceSampling.h \
	PhaseSpaceIntegrandFunction.h PhaseSpaceExactIntegrandFunction.h PhaseSpaceApproximatedIntegrandFunction.h \
	HadronizationException.h GSLRandomNumberGenerator.h ROOTRandomNumberGenerator.h \
	RandomNumberGenerator.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(GSL_INCLUDE) $(RANDOMGEN) $<

PhaseSpaceIntegrandFunction.o: PhaseSpaceIntegrandFunction.cpp PhaseSpaceIntegrandFunction.h 
	$(CPP) -c $(CPPFLAGS) $< 

PhaseSpaceExactIntegrandFunction.o: PhaseSpaceExactIntegrandFunction.cpp PhaseSpaceExactIntegrandFunction.h 
	$(CPP) -c $(CPPFLAGS) $< 

PhaseSpaceApproximatedIntegrandFunction.o: PhaseSpaceApproximatedIntegrandFunction.cpp \
	PhaseSpaceApproximatedIntegrandFunction.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $(LOGLEV) $< 

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

HadronSampling.o: HadronSampling.cpp HadronSampling.h ChargeConfiguration.h \
	HadronData.h Constants.h HadronizationException.h RandomNumberGenerator.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(GSL_INCLUDE) $<

HadronSamplingNew.o: HadronSamplingNew.cpp HadronSamplingNew.h ChargeConfiguration.h \
	HadronData.h Constants.h HadronizationException.h RandomNumberGenerator.h
	$(CPP) -c $(CPPFLAGS) $(ROOT_INCLUDE) $(GSL_INCLUDE) $<

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

logMessage.o: logMessage.cpp logMessage.h
	$(CPP) -c $(CPPFLAGS) $<

computeWaveFunctionStrangeComponent.o: \
	computeWaveFunctionStrangeComponent.cpp computeWaveFunctionStrangeComponent.h
	$(CPP) -c $(CPPFLAGS) $<

fileUtils.o: \
	fileUtils.cpp fileUtils.h
	$(CPP) -c $(CPPFLAGS) $<

mathFunctions.o: mathFunctions.cpp mathFunctions.h HadronizationException.h
	$(CPP) -c $(CPPFLAGS) $<

#########################################################


############### Clean compilation file ##################

cleanall:
	rm -rf partitionFunctionSetBuilder $(OBJECTS)

cleanmcsthar:
	rm -rf $(OBJECTS)

#########################################################
