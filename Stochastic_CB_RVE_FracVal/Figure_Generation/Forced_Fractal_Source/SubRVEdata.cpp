#include <random>
#include <iostream>
#include <algorithm>
#include "RVEdata.h"
#include "SubRVEdata.h"

SubRVEdata::SubRVEdata(){}

SubRVEdata::~SubRVEdata(){}

void SubRVEdata::initialize(long int Lx, long int Ux, long int Ly, long int Uy, long int Lz, long int Uz)
{
    UpperX = Ux;
    LowerX = Lx;
    UpperY = Uy;
    LowerY = Ly;
    UpperZ = Uz;
    LowerZ = Lz;
    
    std::uniform_int_distribution<long int> X_Coord_List(LowerX, UpperX);
    std::uniform_int_distribution<long int> Y_Coord_List(LowerY, UpperY);
    std::uniform_int_distribution<long int> Z_Coord_List(LowerZ, UpperZ);
    Coord_List.push_back(X_Coord_List);
    Coord_List.push_back(Y_Coord_List);
    Coord_List.push_back(Z_Coord_List);
}

void SubRVEdata::Add_Boundary_Particle(long int i)
{
    Sub_Particles.push_back(i);
}

//Return integer distribution list of coordinates within the subRVE
std::vector<std::uniform_int_distribution<long int> > &SubRVEdata::Get_Coord_List()
{
    return Coord_List;
}

std::vector<long int> SubRVEdata::Get_Particles()
{
    return Sub_Particles;
}

void SubRVEdata::PlaceCNTParticles(){}

long int SubRVEdata::GetNumDispFail()
{
    return NumDispFail;
}
