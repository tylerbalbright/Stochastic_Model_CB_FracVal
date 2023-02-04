/* The CB particle class holds information about the particle size, neighbors, and location within the RVE structure. */
#include <vector>

#pragma once
class RVEdata;
class CBparticle
{
private:
    long int ID;     //particle number
    long int index;
    long int xyz[3];     //XYZ coordinates
    double radius;  //defines particle size
    std::vector<long int> neighbors; //holds ID's of particles neighboring this particle
    std::vector<double> NeighborDist; //Holds dist between neigboring particles
    bool Pos;   //connected to positive terminal? (t/f)
    bool Neg;   //connected to negative terminal? (t/f)
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
    
    //Return x, y, and z coordinates
    int GetXCoordinate();
    int GetYCoordinate();
    int GetZCoordinate();
    
    //Find distance to neighbor i (pre-placement)
    double FindDistanceToNeighbor(CBparticle *i, long int x, long int y, long int z);
    
    //Using pre-set particle info (post-placement)
    double FindDistanceToNeighbor(CBparticle *i);
    
    //Using strained particle state
    double FindDistanceToNeighbor(CBparticle *i, float eps1, float eps2, float eps3);
    
    //Return positively or negatively connected status
    bool IsPositive();
    bool IsNegative();
    
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
