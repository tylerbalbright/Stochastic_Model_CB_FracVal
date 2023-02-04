/* The CB particle class holds information about the particle size, neighbors, and location within the RVE structure. */
#include <vector>

#pragma once
class RVEdata;
class CBparticle
{
private:
    long int ID;     //particle ID from the global system
    long int index;  //index of the particle according to the global system
    long int xyz[3];     //XYZ coordinates of the particle
    double radius;  //spherical particle radius
    std::vector<long int> neighbors; //holds ID's of particles neighboring this particle
    std::vector<double> NeighborDist; //Holds dist between neigboring particles
    bool Pos;   //connected to positive terminal? (true/false)
    bool Neg;   //connected to negative terminal? (true/false)
    int MatType;    //Material type integer ID
    int mult;    //magnitude of which the step distance needs multiplied to determine true distance (discretization)
    
public:
    //Constructor and Destructor
    CBparticle();
    ~CBparticle();
    
    //Initialize id of particle
    void Initialize(long int i, int type, double r, int m);
    
    //Set radius to a double value r
    void SetRadius(double r);
    
    //Set index of particle placement
    void Set_Index(long int i);
    
    //Get particle diameter
    double GetDiameter();
    
    //Set x, y, and z location in RVE and check closeness to terminals
    void SetXYZ(long int x, long int y, long int z, long int NumX, int TunnelDist);
    
    //Add neighbor to vectors
    void AddNeighbor(long int i, double d);
    
    //Clear neighbors vectors for strain test
    void Clear_Neighbors();

    //Return x, y, and z coordinates
    int GetXCoordinate();
    int GetYCoordinate();
    int GetZCoordinate();
    
    //Find distance to neighbor i (pre-placement)
    double FindDistanceToNeighbor(CBparticle *i, long int x, long int y, long int z);
    
    //Using pre-set particle info (post-placement)
    double FindDistanceToNeighbor(CBparticle *i);
    
    //Using strained particle state
    double FindStrainDistanceToNeighbor(CBparticle *i, float eps1, float eps2, float eps3);
    
    //Return positively or negatively connected status
    bool IsPositive();
    bool IsNegative();
    
    //Return positively or negatively connected status for strain
    bool Strain_IsPositive(float eps1,long int NumX, int TunnelDist);
    bool Strain_IsNegative(float eps1,long int NumX, int TunnelDist);
    
    //Set Mat Type
    void SetMatType(int i, RVEdata *dat);
    
    //Return vector of neighboring particles
    std::vector<long int> &GetNeighborsVector();
    
    //Return vector of neighboring particles distances
    std::vector<double> &GetNeighborsDistVector();
    
    //Return Mat Type
    int GetMatType();
    
    //Return CB shell and potential maps based on diameter
    void Get_Maps(std::vector<long int> &Shell, std::vector<long int> Potential);
    
};
