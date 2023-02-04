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
    //RVE size definition
    long int Xbound, Ybound, Zbound;
    
    //temp holder of weight percentage filler
    float Wcb;
    
    //Number of CB filler types in RVE
    int NumCBFillers;
    
    //Number of CNT filler types in RVE
    int NumCNTFillers;
    
    //Max particle diameter used for determing boundary volume
    int MaxParticleDiameter;
    
    //The Total_Num_Particles variable is used in multiple functions. For a single particle type system,
    //the usage of this variable is for determining how many particles are to be placed as random aggregates,
    //as forced aggregates, and as forced agglomeraion individual particles. For multi-particle systems, this
    //usage will need to be modified in multiple locations throughout this file specifically.
    long int Total_Num_Particles;
    
    //magnitude by which the step distance needs to be multiplied to acheive true distance in nanometers
    int mult;
    
    //Tunneling distance in nanometers
    int Td;
    
    //unit length in nanometers
    double d;
    
    //Float value between zero and one for agglomeration algorithm
    float DispQuality;
    
    //Percentage of particles to be placed purely randomly and as individaul fill particles
    float Pct_Rand;
    float Pct_Fill;
    
    //Total interferring particles in dataset
    int sumInt = 0;
    
    //Variable for summing the current flowing into the positively connected particles
    double sumCurrent = 0;
    
    //Vectors for determing whether a conductive network exists between the electrodes
    std::vector<bool> network;
    std::vector<std::vector<bool> > posbank;
    std::vector<std::vector<bool> > gndbank;
    
    //Resistance across the RVE network
    double RVEresistance;
    
    //Conductivity of the RVE network
    double RVEconductivity;
    
    //Voltage differential applied across the X+ and X- faces of the RVE
    double VIN = 5;
    
    //Boolean to flag if a network exists in the RVE
    bool NetExists;
    
    //Indexing variable to increment ID #'s for each new particle placed.
    long int running_index;
    
    //Number of forced and random aggregates to be placed
    int num_forced_agg = 0;
    int num_rand_agg = 0;
    
    //Index mapping containers for particle shell and potential indices
    std::vector<std::vector<int> > CB_Potential_Maps_X;
    std::vector<std::vector<int> > CB_Potential_Maps_Y;
    std::vector<std::vector<int> > CB_Potential_Maps_Z;
    std::vector<std::vector<int> > CB_Shell_Maps_X;
    std::vector<std::vector<int> > CB_Shell_Maps_Y;
    std::vector<std::vector<int> > CB_Shell_Maps_Z;
    std::vector<std::vector<long int> > CB_Potential_Maps_Indices;
    std::vector<std::vector<long int> > CB_Shell_Maps_Indices;
    
    //Index mapping containers for strain simulations. Maps are enlarged to allow for scaling to
    //create new interparticle relationships that might've been unrealized in unstrainted config.
    //The potential map is used to determine which indices (surrounding a
    //reference index) have the "potential" to create neighboring relationships with the particle
    //located at the reference index. The shell map is used to determine which indices surrounding
    //the reference index are invalid as a consequnece of their near proximity to the reference index.
    //The maps are stored as vectors of 1-D arrays of integers for corresponding coordinates and indices.
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
    
    //The All_Indices deque is an array of 2-D sparse matrices corresponding to all possible RVE
    //coordinates. The array axis (i.e. All_Indices[i]<-) corresponds to the X coordinate. The sparse matrix
    //coordinates (i.e. All_Indices[i].coeffRef(j,k)) j and k correlate to the Z coordinate and Y coordinate
    //respectively.
    std::deque<Eigen::SparseMatrix<double, Eigen::ColMajor> > All_Indices;
    
    //Containers for storing indices and local coordinates of FracVal aggregate input files
    std::deque<std::vector<std::vector<int> > > Agglomerates;
    std::deque<std::vector<long int> > Agglomerates_Indices;
    std::deque<std::vector<int> > Agglomerates_Bounds;

    //Variable used to define the number of primary particles per aggregate. For all 2021-2022 papers, this
    //variable was maintained at 100 particles per aggregate. It set by reading the # of coordinate lines
    //in each FracVal aggregate file
    int no_p_per_agg;
    
    //Vector to store mass fraction(s) of filler(s) defined via input file
    std::vector<float> MassFractionOfFillers;
    
    //Array of CBMaterialData storage objects
    std::vector<CBMaterialData> CB_data;
    
    //Array of CBparticle types storage objects
    std::deque<CBparticle> CB_particles;

    //Array of CNT particle types storage objects
    std::vector<CNTparticle> CNT_particles;
    
    //Object storing info about polymer
    PolymerData p_data;

    //vector storing number of each type of particle in RVE
    std::vector<long int> NumParticles;

    //Vector of particle objects holding details about the particle types and their local/global IDs
    std::vector<Particles> particle_data;

    //Vector of nodal voltages of all particles identified as part of conductive network between electrodes
    std::vector<double> Nodal_Voltages;
    
public:
    //Constructor and Destructor
    RVEdata();
    ~RVEdata();
    
    //Function to read input file
    void ReadInputData(std::ifstream &fin, std::ofstream &ferr, int cmd_PctRandLower, int cmd_PctRandUpper, int cmd_XY, int cmd_thickness, int cmd_wcb, int cmd_fill);
    
    //Function to read agglomerate files
    void Read_Agg_Files();
    
    //Return positive-ground distance (limit)
    long int GetXBoundaryLimit();
    
    //Function that checks vector N for particle O
    bool CheckInVector(std::vector<long int> &ParticleN, long int ParticleIndexO);
    
    //Main function for building the RVE (placing particles)
    void Build_RVE();
    
    //Place the remainder of the particles next to the seed particles
    void Place_Remainder(long int num_rand);
    
    //Place random (seed) aggregates in the RVE
    void Place_Rand_Agglomerates(long int num_rand);
    
    //Place forced agglomeration aggregates in the RVE
    void Place_Forced_Agglomerates(long int num_rand);
    
    //Find the next available placement location
    void Get_Next_Location(long int &next);
    
    //Initialize a new particle at the given index
    void Initialize_Index(long int next_index, long int ID);
    
    //Add the new CB particles neighboring indices to the placement list
    void Map_CB(long int CB);
    
    //Function to determine the index of the correct map for a particle
    int Get_CB_Map_Index(long int ID);

    //Function to return the index of the largest particle index map (used in strain test)
    int Get_Max_CB_Map_Index(long int ID);
    
    //Functions to initialize shell and potential maps
    void Initialize_CB_Potential_Shell_Maps();
    void Initialize_CB_Strain_Potential_Maps();
    
    //Return number of particles in the RVE
    int GetNumParticles();
    
    //Return number of CB fillers
    int GetNumCBFillers();
    
    //Check if sparse matrix position is null or filled
    bool isNull(const Eigen::SparseMatrix<double, Eigen::ColMajor>& mat, int row, int col);
    
    //Return pointer to particle i
    Particles* GetParticle(long int i);

    //Return pointer to CB particle #i
    CBparticle* GetCBParticle(long int num);
    
    //Return pointer to CNT particle #i
    CNTparticle* GetCNTParticle(long int num);
    
    //Return pointer to material type object
    CBMaterialData* GetCBMatType(int i);
    
    //Return pointer to polymer data object
    PolymerData* GetPolymerData();
    
    //Function to find the neighboring pairs of particles in the RVE
    void Find_Neighbors();
    
    //Function to find new neighboring pairs with applied strain
    void Find_Strained_Neighbors(float strain_elong);
    
    //Clear neighboring pair data containers for strain tests
    void Reset_Neighbors();
    
    //Functions to check adjacent indices for particles and add as neighbors
    void Check_Adjacent_Indices(long int i);

    //Function to check adjacent indices for particles during strain testing
    void Check_Adjacent_Strain_Indices(long int i, float eps1, float eps2, float eps3);
    
    //Function to solve for the conductivity of the RVE
    bool SolveForConductivity();
    
    //Function to solve for the conductivity of the strained RVE
    bool SolveForStrainConductivity(float strain_elong, std::ofstream &cond_out);
    
    //Function to solve for the nodal voltages of the networked RVE particles
    void SolveForNodalVoltages();
    
    //Function to parallelize the filling of the conductivity and current matrices
    std::vector<Eigen::Triplet<double> > FillCurrAndG(std::vector<double> &ValuesCurr, std::vector<long int> &Indices, std::vector<long int> &NeighborIndices, double Rc, int VIN, long int start, long int finish);

    //Function to find the network via positive -> negative algorithm
    bool FindNetwork(long int N, std::deque<std::mutex> &PosLock, std::deque<std::mutex> &GndLocks);
    
    //Function to find the strained network via positive -> negative algorithm
    bool FindStrainNetwork(long int N, std::deque<std::mutex> &PosLock, std::deque<std::mutex> &GndLocks, float eps1);
    
    //Searching functions for parrallelized network finding algorithm
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
    
    //Function to output results
    void OutputResults(std::ofstream &fout);

    //Function to output RVE coordinate file in Ovito accessible format
    void Output_Ovito_File(std::ofstream &fplot, std::string uid, std::string bid, std::ofstream &cond_out);
};
