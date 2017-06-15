#include "Common.h"
#include <cstdlib>
double NeuroNet::Common::getRand(double vmin, double vmax)
{
	return vmin + rand() * (vmax - vmin) / RAND_MAX;//
	//return 0.12345;
}
