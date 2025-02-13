#include "../Include/buildParameterRegularGrid.h"
#include "../Include/HadronizationException.h"
// TODO: check and update all method/function name
void buildParameterRegularGrid(const double i_minimumParameterValue,
                               const double i_maximumParameterValue,
                               const unsigned int i_numberParameterValues,
                               vector<double>& io_parameterRegularGrid)
{
    // Check grid structure
    if(i_minimumParameterValue>i_maximumParameterValue)
    {
        throw HadronizationException("Error during one dimensional parameter regular grid building, minimum/maximum value ordering error",
                                     __FUNCTION__,731);
    }
    if(i_numberParameterValues==0)
    {
        throw HadronizationException("Error during one dimensional parameter regular grid building, number of required points in grid is zero",
                                     __FUNCTION__,732);
    }
    if((i_minimumParameterValue==i_maximumParameterValue) && (i_numberParameterValues>1))
    {
        throw HadronizationException("Error during one dimensional parameter regular grid building, parameter range is zero and more than one value has been required",
                                     __FUNCTION__,733);
    }
    
    
    // Build grid
    double parameterStep(0.);
    io_parameterRegularGrid.resize(i_numberParameterValues,0.);
    if(i_numberParameterValues>1)
    {
        parameterStep = (i_maximumParameterValue - i_minimumParameterValue)/(i_numberParameterValues-1);
    }
    double parameterValue(i_minimumParameterValue);
    io_parameterRegularGrid[0] = parameterValue;
    for(unsigned int gridPointIndex=1;gridPointIndex<i_numberParameterValues;++gridPointIndex)
    {
        parameterValue += parameterStep;
        io_parameterRegularGrid[gridPointIndex] = parameterValue;
    }
    
    return;
}
