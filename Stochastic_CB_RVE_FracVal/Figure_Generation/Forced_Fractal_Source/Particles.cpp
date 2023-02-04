#include "Particles.h"

Particles::Particles()
{
    Type = -1;
    TypeID = -1;
    GlobalID = -1;
    PosTerminalResistance = 0;
}
Particles::~Particles(){}

void Particles::Initialize(long int t, long int tID, long int gID)
{
    Type = t;
    TypeID = tID;
    GlobalID = gID;
    PosTerminalResistance = 0;
}

long int Particles::GetType()
{
    return Type;
}

long int Particles::GetTypeID()
{
    return TypeID;
}

long int Particles::GetGlobalID()
{
    return GlobalID;
}

double Particles::GetPosTerminalResistance()
{
    return PosTerminalResistance;
}

void Particles::SetPosTerminalResistance(double r)
{
    PosTerminalResistance = r;
}
