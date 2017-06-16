#include "Common.h"
#include <cstdlib>
double NeuroNet::Common::getRand(double vmin, double vmax, bool integer)
{
	if (integer)
		return (int)(vmin + rand() * (vmax - vmin) / RAND_MAX);
	return vmin + rand() * (vmax - vmin) / RAND_MAX;//
}
