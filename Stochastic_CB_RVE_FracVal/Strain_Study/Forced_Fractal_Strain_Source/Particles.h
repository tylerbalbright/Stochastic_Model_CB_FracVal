//This class holds the particle information such as particle type, and subsequent particle ID based on that type specifier.

#pragma once

class Particles
{
private:
    long int Type;
    long int TypeID;
    long int GlobalID;
    double PosTerminalResistance;
    
public:
    //Constructor and Destructor
    Particles();
    ~Particles();
    
    void Initialize(long int t, long int tID, long int gID);
    long int GetType();
    long int GetTypeID();
    long int GetGlobalID();
    void SetPosTerminalResistance(double r);
    double GetPosTerminalResistance();
};
