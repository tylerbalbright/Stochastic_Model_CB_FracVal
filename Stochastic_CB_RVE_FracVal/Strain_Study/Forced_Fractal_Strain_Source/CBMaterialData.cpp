#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <iostream>
#include "CBMaterialData.h"
#include "time.h"

CBMaterialData::CBMaterialData()
{
	ID = 0;
	Density = 0;
	diamLowerBound = 0;
	diamUpperBound = 0;
	PSoftness = 0;
	NumberOfParticlesInRVE = 0;
};

CBMaterialData::~CBMaterialData(){};

void CBMaterialData::initialize(int i, double rho, int l, int u, int ps, double d)
{
    ID = i;
    Density = rho*1000.;
    diamLowerBound = l;
    diamUpperBound = u;
    PSoftness = ps;
}

double CBMaterialData::GetDensity()
{
    return Density;
}

long double CBMaterialData::GetAvgParticleWeight(double d)
{
    //Mass of 1 particle based on a spherical shape and bulk density properties from manuf.
    long double M1p = (4./3.)*M_PI*(pow((diamLowerBound+diamUpperBound)*d/4.,3))*Density;
    //return average particle weight
    return M1p;
}

void CBMaterialData::SetNumParticles(int i)
{
    NumberOfParticlesInRVE = i;
}

int CBMaterialData::GetNumParticles()
{
    return NumberOfParticlesInRVE;
}

int CBMaterialData::GetPenetrationAllowance()
{
    return PSoftness;
}

double CBMaterialData::GetRandomRadius()
{
    double d = 0;
    std::uniform_int_distribution<unsigned int> distribution(diamLowerBound, diamUpperBound);
    srand(time(NULL));
    std::mt19937_64 gen(rand());
    d = distribution(gen)/2.0;
    
    return d;
}

int CBMaterialData::GetMaxDiameter()
{
    return diamUpperBound;
}

int CBMaterialData::GetMinDiameter()
{
    return diamLowerBound;
}
