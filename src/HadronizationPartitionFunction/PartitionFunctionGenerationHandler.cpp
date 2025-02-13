#include "../Include/PartitionFunctionGenerationHandler.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/buildParameterRegularGrid.h"

PartitionFunctionGenerationHandler::PartitionFunctionGenerationHandler(const string& i_hadronDataSetFileName,
                                                                       const double i_lightHadronMaxMass,
                                                                       const double i_minEnergyDensity,
                                                                       const double i_maxEnergyDensity,
                                                                       const double i_numberOfEnergyDensityValues,
                                                                       const vector<double>& i_massValueList,
                                                                       const MCSTHAR::Utilities::ChargeConfiguration& i_chargeConfiguration)
                                                                      :m_isHandlerAvailable(false)
                                                                      ,m_hadronSamplingGroups(i_hadronDataSetFileName,
                                                                                              i_lightHadronMaxMass)
                                                                      ,m_massValues(i_massValueList)
                                                                      ,m_chargeConfiguration(i_chargeConfiguration)
{
    // TODO: check if all calculation setup data are being used, the remaining can be removed
    try
    {
        // TODO: add missing checks
        if(m_massValues.size()==0)
        {
            throw HadronizationException("Error during partition function calculation run, empty mass value set provided",
                                         __FUNCTION__,831);
        }
        
        // TODO: this operation is probably partially done in m_partitionFunctionWriter builder, factorize!
        // Build energy density value set
        buildParameterRegularGrid(i_minEnergyDensity,
                                  i_maxEnergyDensity,
                                  i_numberOfEnergyDensityValues,
                                  m_energyDensityValues);
    }
    catch(HadronizationException& ex)
    {
        throw ex;
    }
}

// TODO: all these empty destructors can be avoided!
PartitionFunctionGenerationHandler::~PartitionFunctionGenerationHandler(void)
{
}
