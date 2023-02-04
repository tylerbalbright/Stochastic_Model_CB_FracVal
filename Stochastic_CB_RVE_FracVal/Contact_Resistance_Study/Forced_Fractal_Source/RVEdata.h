/* The RVE data class holds all of the data and major functions required to store the runtime data. */
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <mutex>
#include <deque>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include "Particles.h"
#include "CBparticle.h"
#include "CNTparticle.h"
#include "CBMaterialData.h"
#include "PolymerData.h"
#pragma once

class SubRVEdata;
class RVEdata
{
private:
    long int Xbound, Ybound, Zbound;   //RVE size definition
    
    //Defining the volume of the "internal" portion of the sub-RVE's
    float Wcb;                      //temp holder of weight percentage filler
    int NumCBFillers;             //Number of CB filler types in RVE
    int NumCNTFillers;           //Number of CNT filler types in RVE
    int MaxParticleDiameter;                      //Max particle diameter used for determing boundary volume
    long int Total_Num_Particles;
    int mult;                                     //magnitude by which the step distance needs to be multiplied to acheive true distance in nanometers
    int Td;                                       //Tunneling distance in nanometers
    double d;                                     //unit size in nanometers
    float DispQuality;                            //Float value between zero and one
    float Pct_Rand;                               //Percentage of particles to be placed purely randomly
    float Pct_Fill;
    int sumInt = 0;                               //Total interferring particles in dataset
    double sumCurrent = 0;                        //Summing the current flowing into the positively connected particles
    std::vector<bool> network;
    std::vector<std::vector<bool> > posbank;
    std::vector<std::vector<bool> > gndbank;
    double RVEresistance;                         //Resistance across the RVE network
    double RVEconductivity;                       //Conductivity of the RVE network
    double VIN = 5;                               //Voltage differential applied across the RVE
    bool NetExists;                               //Boolean to flag if a network exists in the RVE
    long int running_index;
    int num_forced_agg = 0;                       // Number of forced agglomerates to be placed
    int num_rand_agg = 0;
    
    double Rc = 0;
    
    //Mapping and containers for particle shell and potential indices
    std::vector<std::vector<int> > CB_Potential_Maps_X;
    std::vector<std::vector<int> > CB_Potential_Maps_Y;
    std::vector<std::vector<int> > CB_Potential_Maps_Z;
    std::vector<std::vector<int> > CB_Shell_Maps_X;
    std::vector<std::vector<int> > CB_Shell_Maps_Y;
    std::vector<std::vector<int> > CB_Shell_Maps_Z;
    std::vector<std::vector<long int> > CB_Potential_Maps_Indices;
    std::vector<std::vector<long int> > CB_Shell_Maps_Indices;
    
    //Mapping and containers for strain calc stuff
    std::vector<std::vector<int> > CB_Strain_Potential_Maps_X;
    std::vector<std::vector<int> > CB_Strain_Potential_Maps_Y;
    std::vector<std::vector<int> > CB_Strain_Potential_Maps_Z;
    std::vector<std::vector<long int> > CB_Strain_Potential_Maps_Indices;
    std::vector<std::vector<int> > CB_Strain_Shell_Maps_X;
    std::vector<std::vector<int> > CB_Strain_Shell_Maps_Y;
    std::vector<std::vector<int> > CB_Strain_Shell_Maps_Z;
    std::vector<std::vector<long int> > CB_Strain_Shell_Maps_Indices;

    //Bitfields for indexing particles and shells
    boost::dynamic_bitset<> Potential;
    boost::dynamic_bitset<> Shell;
    std::deque<Eigen::SparseMatrix<double, Eigen::ColMajor> > All_Indices;
    
    //Agglomerate Banks
    std::deque<std::vector<std::vector<int> > > Agglomerates;
    std::deque<std::vector<long int> > Agglomerates_Indices;
    std::deque<std::vector<int> > Agglomerates_Bounds;
    
    std::vector<float> MassFractionOfFillers;     //Define volume fractions of filler(s)
    std::vector<CBMaterialData> CB_data;          //Array of CBMaterialData storage objects
    std::deque<CBparticle> CB_particles;          //2-D Array of CBparticle types storage objects
    std::vector<CNTparticle> CNT_particles;       //2-D Array of CNT particle types storage objects
    PolymerData p_data;                           ///Object storing info about polymer
    std::vector<long int> NumParticles;           //vector storing number of each type of particle in RVE
    std::vector<Particles> particle_data;         //Vector of particle objects holding details about the particle types and their local/global IDs
    std::vector<double> Nodal_Voltages;
    std::vector<std::vector<double> > ConductivityMatrix;
    int no_p_per_agg;
    
    

public:
    //Constructor and Destructor
    RVEdata();
    ~RVEdata();
    
    //Function to read input file
    void ReadInputData(std::ifstream &fin, std::ofstream &ferr, int cmd_PctRandLower, int cmd_PctRandUpper, int cmd_XY, int cmd_thickness, int cmd_wcb, int cmd_fill, double cmd_Rc);
    
    //Function to read agglomerate files
    void Read_Agg_Files();
    
    //Return positive-ground distance (limit)
    long int GetXBoundaryLimit();
    
    //Function that checks vector N for particle O
    bool CheckInVector(std::vector<long int> &ParticleN, long int ParticleIndexO);
    
    //Main function for building the RVE
    void Build_RVE();
    
    //Place the remainder of the particles next to the seed particles
    void Place_Remainder(long int num_rand);
    
    //Place agglomerates in the RVE
    void Place_Rand_Agglomerates(long int num_rand);
    
    //Place agglomerates in the RVE
    void Place_Forced_Agglomerates(long int num_rand);
    
    //Find the next available placement location
    void Get_Next_Location(long int &next);
    
    //Initialize a new particle at the given index
    void Initialize_Index(long int next_index, long int ID);
    
    //Populate the SubRVEs with the prescribed particles defined in the input file
    void Place_Seeds();
    
    //Add the new CB particles neighboring indices to the placement list
    void Map_CB(long int CB);
    
    int Get_CB_Map_Index(long int ID);
    
    int Get_Max_CB_Map_Index(long int ID);
    
    void Initialize_CB_Potential_Shell_Maps();
    void Initialize_CB_Strain_Potential_Maps();
    
    //Return number of particles
    int GetNumParticles();
    
    //Return number of CB fillers
    int GetNumCBFillers();
    
    //Check if sparse matrix position is null or filled
    bool isNull(const Eigen::SparseMatrix<double, Eigen::ColMajor>& mat, int row, int col);
    
    //Return point to particle i
    Particles* GetParticle(long int i);

    //Return pointer to CB particle #i
    CBparticle* GetCBParticle(long int num);
    
    //Return pointer to CNT particle #i
    CNTparticle* GetCNTParticle(long int num);
    
    //Return pointer to material type object
    CBMaterialData* GetCBMatType(int i);
    
    //Return pointer to polymer data object
    PolymerData* GetPolymerData();
    
    //Find the neighboring particles
    void Find_Neighbors();
    
    //Find new neighbor relations with applied strains
    void Find_Strained_Neighbors(float strain_elong);
    
    //Clear neighbors lists for strain tests
    void Reset_Neighbors();
    
    //Check adjacent indices for particles and add as neighbors
    void Check_Adjacent_Indices(long int i);
    bool Seed_Check_Adjacent_Indices(long int x, long int y, long int z, long int i);
    
    //Check adjacent particles for strain analysis
    void Check_Adjacent_Strain_Indices(long int i, float eps1, float eps2, float eps3);
    
    //Function to calculate the conductivity of the RVE
    bool SolveForConductivity();
    
    //Function to calculate the conductivity of the RVE
    bool SolveForStrainConductivity(float strain_elong, std::ofstream &cond_out);
    
    //Function to solve for the nodal voltages of the networked RVE particles
    void SolveForNodalVoltages();
    
    //Function to thread the filling of the conductivity and current matrices
    std::vector<Eigen::Triplet<double> > FillCurrAndG(std::vector<double> &ValuesCurr, std::vector<long int> &Indices, std::vector<long int> &NeighborIndices, double Rc, int VIN, long int start, long int finish);

    //Find the network via positive -> negative algorithm
    bool FindNetwork(long int N, std::deque<std::mutex> &PosLock, std::deque<std::mutex> &GndLocks);
    
    //Find the strained network via positive -> negative algorithm
    bool FindStrainNetwork(long int N, std::deque<std::mutex> &PosLock, std::deque<std::mutex> &GndLocks, float eps1);
    
    //Branch searching functions for parrallel find net algorithm
    void NegativeBranch(int NumProcessors, std::deque<std::mutex> &GndLocks);
    void PositiveBranch(int NumProcessors, std::deque<std::mutex> &PosLocks);
    void PositiveSearch(std::deque<std::mutex> &PosLocks);
    void NegativeSearch(std::deque<std::mutex> &GndLocks);

    //Return instance of neighbors vector for particle i
    std::vector<long int> &GetNeighbors(long int i);
    
    //Return instance of neighbors distance vector for particle i
    std::vector<double> &GetNeighborsDist(long int i);

    //Return tunneling threshold between two particles
    double CalcTunnelingThreshold(long int i, long int j);
    
    //Return tunneling threshold between particle and terminal
    double CalcTunnelingThresholdToTerminal(long int i);
    
    //Return cross-sectional area of tunneling path
    long double FindTunnelingArea(long int i, long int j);

    //Return tunneling resistance between two particles
    double CalcTunnelingResistance(double dist, double TunnelingThresh, double u, long double At);
    
    //Return distance from particle to terminal (pos or gnd defined by boolean)
    double  FindDistanceToTerminal(bool terminal, long int i);
    
    //Check if particle is close to pos or neg terminal
    bool CheckPosNeg(long int i, bool j);
    
    //Error reporting
    void CheckNeighbors();
    
    void OutputResults(std::ofstream &fout);
    
    void Output_Ovito_File(std::ofstream &fplot, std::string uid, std::string bid, std::ofstream &cond_out);
};
