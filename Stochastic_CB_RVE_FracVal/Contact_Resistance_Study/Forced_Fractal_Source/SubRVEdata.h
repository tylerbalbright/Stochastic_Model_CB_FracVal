/* Sub RVE data storage object. */
#include <vector>
#pragma once

class RVEdata;
class SubRVEdata
{
private:
    int ID;
    long int NumDispFail;
    std::vector<long int> Sub_Particles;
    long int UpperX, LowerX, UpperY, LowerY, UpperZ, LowerZ;
    std::vector<std::uniform_int_distribution<long int> > Coord_List;
    
public:
    //Constructor and destructor
    SubRVEdata();
    ~SubRVEdata();
    
    //Initialize instance of Sub RVE
    void initialize(long int Ux, long int Lx, long int Uy, long int Ly, long int Uz, long int Lz);
    
    //function to add boundary particle to domain particle list
    void Add_Boundary_Particle(long int i);
    
    //Return integer distribution list of coordinates within the subRVE
    std::vector<std::uniform_int_distribution<long int> > &Get_Coord_List();
    
    //Return vector of particles that are in the sub domain
    std::vector<long int> Get_Particles();
    
    //Place CNT particles in the sub domain
    void PlaceCNTParticles();
    
    long int GetNumDispFail();
};
