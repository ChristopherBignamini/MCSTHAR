#ifndef BUILDPARAMETERREGULARGRID_H
#define BUILDPARAMETERREGULARGRID_H

#include <vector>

using namespace std;

/**
* @brief Single parameter regular (one dimensional) grid building method
*
* @param i_minimumParameterValue Minimum parameter value
* @param i_maximumParameterValue Maximum parameter value
* @param i_numberParameterValues Number of parameter values in grid
* @param io_parameterRegularGrid Parameter regular (one dimensional) grid (I/O parameter)
* @throw HadronizationException in case of min/max parameter value inconsistency or zero number of parameter values
*
* @author Christopher Bignamini
*/
void buildParameterRegularGrid(double i_minimumParameterValue,
                               double i_maximumParameterValue,
                               unsigned int i_numberParameterValues,
                               vector<double>& io_parameterRegularGrid);
#endif