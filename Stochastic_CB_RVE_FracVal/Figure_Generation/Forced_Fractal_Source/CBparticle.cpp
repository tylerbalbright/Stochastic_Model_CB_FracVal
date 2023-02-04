#include "RVEdata.h"
#include "PolymerData.h"
#include "CBparticle.h"
#include <cmath>
#include <random>
#include <iostream>

CBparticle::CBparticle()
{
    Pos = false;
    Neg = false;
    //radius = 15;
}

void CBparticle::Initialize(long int i, int type, double r, int m)
{
    ID = i;
    Pos = false;
    Neg = false;
    radius = r;
	MatType = type;
    mult = m;
}

CBparticle::~CBparticle()
{
    neighbors.clear();
}

void CBparticle::SetRadius(double r)
{
    radius = r;
}

void CBparticle::Set_Index(long int i)
{
    index = i;
}

void CBparticle::SetXYZ(long int x, long int y, long int z, long int NumX, int TunnelDist)
{
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    
    //Calc distance to terminals
    double DistPos = 0, DistNeg = 0;
    DistNeg = (NumX - static_cast<double> (xyz[0]))*mult;
    DistPos = static_cast<double> (xyz[0])*mult;
    
    if (DistPos <= ((radius*mult)+TunnelDist))
    {
        Pos = true;
    }
    if (DistNeg <= ((radius*mult)+TunnelDist))
    {
        Neg = true;
    }
}

void CBparticle::AddNeighbor(long int i, double d)
{
    neighbors.push_back(i);
    NeighborDist.push_back(d);
}

int CBparticle::GetXCoordinate()
{
    return xyz[0];
}

int CBparticle::GetYCoordinate()
{
    return xyz[1];
}

int CBparticle::GetZCoordinate()
{
    return xyz[2];
}

double CBparticle::FindDistanceToNeighbor(CBparticle *i)
{
    long int iX, iY, iZ;
    iX = i->GetXCoordinate();
    iY = i->GetYCoordinate();
    iZ = i->GetZCoordinate();
    
    double distance;
    distance = mult * sqrt(std::abs(((iX-xyz[0])*(iX-xyz[0]))+((iY-xyz[1])*(iY-xyz[1]))+((iZ-xyz[2])*(iZ-xyz[2]))));
    return distance;
}

double CBparticle::FindDistanceToNeighbor(CBparticle *i, long int x, long int y, long int z)
{
    long int iX, iY, iZ;
    iX = i->GetXCoordinate();
    iY = i->GetYCoordinate();
    iZ = i->GetZCoordinate();
    
    double distance;
    distance = mult * sqrt(std::abs(((iX-x)*(iX-x))+((iY-y)*(iY-y))+((iZ-z)*(iZ-z))));
    return distance;
}

double CBparticle::FindDistanceToNeighbor(CBparticle *i, float eps1, float eps2, float eps3)
{
    int iX, iY, iZ;
    iX = i->GetXCoordinate();
    iY = i->GetYCoordinate();
    iZ = i->GetZCoordinate();
    
    double distance;
    distance = mult * sqrt(std::abs(((iX-xyz[0])*eps1*(iX-xyz[0])*eps1)+((iY-xyz[1])*eps2*(iY-xyz[1])*eps2)+((iZ-xyz[2])*eps3*(iZ-xyz[2])*eps3)));
    return distance;
}

void CBparticle::SetMatType(int i, RVEdata *dat)
{
    MatType = i;
    SetRadius(dat->GetCBMatType(i)->GetRandomRadius());
}

int CBparticle::GetMatType()
{
    return MatType;
}

std::vector<long int>& CBparticle::GetNeighborsVector()
{
    return neighbors;
}

std::vector<double>& CBparticle::GetNeighborsDistVector()
{
    return NeighborDist;
}


double CBparticle::GetDiameter()
{
    //return 30;
    return radius*2.0;
}

bool CBparticle::IsPositive()
{
    return Pos;
}

bool CBparticle::IsNegative()
{
    return Neg;
}
