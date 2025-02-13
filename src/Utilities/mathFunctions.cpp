#include "../Include/mathFunctions.h"
#include "../Include/HadronizationException.h"

double factorial(const unsigned int i_n)// TODO: switch to look up table method
{
    // Check if factorial can be represented as IEEE double precision (max is 170!)
    if(i_n<171)
    {
        double o_factorialValue(1.);
        for(unsigned int integerIndex=2;integerIndex<=i_n;++integerIndex)
        {
            o_factorialValue *= integerIndex;
        }
        
        return o_factorialValue;
    }
    else
    {
        throw HadronizationException("Error during factorial calculation, provided input integer larger than 170",
                                     __FUNCTION__,751);
    }
}
