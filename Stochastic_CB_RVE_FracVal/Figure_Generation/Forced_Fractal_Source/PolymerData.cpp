#include <cmath>
#include "PolymerData.h"

PolymerData::PolymerData(){};

PolymerData::~PolymerData(){};

void PolymerData::initialize(double rho, long double pbh, int TunnDist, float p, double d)
{
    Density = rho*1000.; //*d*d*d/pow(pow(10,-9),3);
    PolymerBarrierHeight = pbh;
    TunnelingDistance = TunnDist;
    nu = p;
}

double PolymerData::GetDensity()
{
    return Density;
}

float PolymerData::GetPoissonsRatio()
{
    return nu;
}
