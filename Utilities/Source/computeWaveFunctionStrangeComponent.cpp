#include <cmath>
#include "../Include/Constants.h"
#include "../Include/computeWaveFunctionStrangeComponent.h"
#include "../Include/HadronizationException.h"

double computeWaveFunctionStrangeComponent(const MesonSeries i_mesonSeries,const double i_mixingAngle)
{
	const double radMixingAngle(i_mixingAngle*piOn180);
	
	if(i_mesonSeries == Singlet)
    {
		return cos(radMixingAngle + constantMixingAngle)*cos(radMixingAngle + constantMixingAngle);
	}
	else if(i_mesonSeries == Octect)
    {
		return sin(radMixingAngle + constantMixingAngle)*sin(radMixingAngle + constantMixingAngle);
	}
    else
    {
        throw HadronizationException("Error in wave function strange component calculation, not identified meson series",
                                     __FUNCTION__,701);
    }
}
