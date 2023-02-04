/* Define material properties for the RVE */

#pragma once

class CBMaterialData
{
private:
    int ID;
    double Density;
    int diamLowerBound, diamUpperBound;
    int NumberOfParticlesInRVE;
    int PSoftness;
    
public:
    //Constructors and Destructors
    CBMaterialData();
    ~CBMaterialData();
    
    //Initialize material properties
    void initialize(int i, double rho, int l, int u, int ps, double d);
    
    //Return the density of the material
    double GetDensity();
    
    //Return average particle weight
    long double GetAvgParticleWeight(double d);
    
    //Set number of particles based on weight fraction
    void SetNumParticles(int i);
    
    //Get the number of particles of that type
    int GetNumParticles();
    
    //Get particle softness (penetration allowance)
    int GetPenetrationAllowance();
    
    //Return random radius for CB particle
    double GetRandomRadius();
    
    //Return max diameter of particle
    int GetMaxDiameter();
    
    //Return min diameter of particle
    int GetMinDiameter();
};
