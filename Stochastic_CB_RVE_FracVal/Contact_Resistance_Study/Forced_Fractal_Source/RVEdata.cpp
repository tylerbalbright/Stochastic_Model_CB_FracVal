#define _USE_MATH_DEFINES

//#define WINDOWS

#ifdef WINDOWS
#include <direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#define getcwd getcwd
#endif

#include <iostream>
#include <fstream>
#include <future>
#include <mutex>
#include <deque>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include <algorithm>
#include "Files.h"
#include "PolymerData.h"
#include "CBparticle.h"
#include "CBMaterialData.h"
#include "SubRVEdata.h"
#include "visit_writer.h"
#include <Eigen/Core>
//#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include "RVEdata.h"

RVEdata::RVEdata(){};

RVEdata::~RVEdata(){};

//Check vectors to see if they hold a similar integer
bool RVEdata::CheckInVector(std::vector<long int> &ParticleN, long int ParticleIndexO)
{
    if (std::find(ParticleN.begin(), ParticleN.end(), ParticleIndexO) != ParticleN.end())
    {
        return (true);
    }
    else
    {
        return (false);
    }
}

void RVEdata::ReadInputData(std::ifstream &fin, std::ofstream &ferr, int cmd_PctRandLower, int cmd_PctRandUpper, int cmd_XY, int cmd_thickness, int cmd_wcb, int cmd_fill, double cmd_Rc)
{
    //Define Rc from sim data
    Rc = cmd_Rc;
    
    sumInt = 0;
    double NotNeeded;
    std::string lines, ignore;
    
    //Read in basic problem parameters
    fin >> Xbound >> Ybound >> Zbound >> NumCBFillers >> NumCNTFillers >> Pct_Rand >> d;
    
    //Define the multiplier for converting unit steps to distance
    mult = round(d/pow(10,-9));
    
    //********** Modified Input for Pct Rand and thickness **************
    //Initialize random number generator instance for determining params
    // Convert both the integers to string
    std::string s1 = std::to_string(time(NULL));
    std::string s2 = std::to_string(getpid());
    // Concatenate both strings
    std::string s = s1 + s2;
    // Convert the concatenated string
    // to integer
    long int c = std::stol(s);
    srand(c);
    std::mt19937_64 gen(rand());
    std::uniform_int_distribution<int> Qdist(cmd_PctRandLower, cmd_PctRandUpper);
    Pct_Rand = Qdist(gen)/10000000.;
    Pct_Fill = cmd_fill/100.;
    Zbound = cmd_thickness;
    Xbound = cmd_XY;
    Ybound = cmd_XY;
    float Real_Wcb = cmd_wcb/100.;
    //Read in and create CB filler material(s) object(s)
    CB_data.resize(NumCBFillers);
    
    //std::cout << "mult = " << mult << std::endl;
    
    int ID, diamUpper, diamLower, ps;
    double rho;
    MaxParticleDiameter = 0;
    for (int i = 0; i < NumCBFillers; i++)
    {
        fin >> ID >> Wcb >> rho >> diamLower >> diamUpper >> ps;
        Wcb = Real_Wcb;
        CB_data[i].initialize(ID, rho, diamLower, diamUpper, ps, d);
        MassFractionOfFillers.push_back(Wcb);
        
        //Check if CB filler i has a larger max radius than i-1
        if (diamUpper > MaxParticleDiameter)
        {
            MaxParticleDiameter = diamUpper;
        }
    }
    
    for (int i = 0; i < NumCNTFillers; i++)
    {
        //Read CNT data
    }
    
    long double PolymerBarrierHeight;
    int TunnelingDist;
    float p;
    fin >> rho >> PolymerBarrierHeight >> TunnelingDist >> p;
    p_data.initialize(rho, PolymerBarrierHeight, TunnelingDist, p, d);
    Td = TunnelingDist*mult;
    
    //Determine quantity of each particle type based on weight percentage
    float Wp = 1;       //polymer at 100% mass fraction
    for (int i = 0; i < MassFractionOfFillers.size(); i++)
    {
        Wp -= MassFractionOfFillers[i]/100.0;      //reduced by mass fraction occupied by filler particles
    }
    long double cells = Xbound*Ybound*Zbound; //static_cast<long double>((Xbound)*(Ybound)*(Zbound));
    long double Vrve = cells * (pow(d, 3));     //volume of RVE in nm^3
    //std::cout << "d = " << d << std::endl;
    //std::cout << "cells = " << cells << ", Vrve = " << Vrve << std::endl;
    
    //Define the indices holding deque
    for (int t = 0; t <= Xbound; t++)
    {
        All_Indices.push_back(Eigen::SparseMatrix<double, Eigen::ColMajor>(Ybound+1, Zbound+1));
    }
    
    //Find density of composite using weight fractions https://nptel.ac.in/courses/101106038/mod03lec01.pdf
    double CompRho = Wp/p_data.GetDensity();
    for (int i = 0; i < NumCBFillers; i++)
    {
        CompRho += (MassFractionOfFillers[i]/100.0)/CB_data[i].GetDensity();
    }
    
    //Calculate number of particles of each type and add to CB_data
    long int particle_count = 0;
    long int CB_count = 0;
    long int CNT_count = 0;

    //Calculating the number of particles to be included in the "sub-Rves"
    Total_Num_Particles = 0;
    
    for (int i = 0; i<NumCBFillers; i++)
    {
        NumParticles.push_back(round((MassFractionOfFillers[i]/100.0)*Vrve*(1.0/CompRho)/CB_data[i].GetAvgParticleWeight(d)));
        //std::cout << "Num_p = " << NumParticles[i] << std::endl;
        for (long int j = 0; j < NumParticles[i]; j++)
        {
            particle_data.emplace_back();
            particle_data[particle_count].Initialize(i, CB_count, particle_count);
			CB_particles.emplace_back();
			CB_particles[CB_count].Initialize(CB_count, i, CB_data[i].GetRandomRadius(), round(d/pow(10,-9)));
			CB_count += 1;
            particle_count += 1;
        }
        CB_data[i].SetNumParticles(NumParticles[i]);
        Total_Num_Particles += NumParticles[i];
    }
    
    //std::cout << "Num Particles = " << Total_Num_Particles << std::endl;
    
    //Define bit flag field sizes for shells and potentials
    Potential.resize((Xbound+1)*(Ybound+1)*(Zbound+1));
    Potential.reset();
    Shell.resize((Xbound+1)*(Ybound+1)*(Zbound+1));
    Shell.reset();
    //std::cout << "Shell resized" << std::endl;
    //Resize number of boundary particle locks and initialize CB particles=
    for (int i = 0+NumCBFillers; i<NumCBFillers+NumCNTFillers; i++)
    {
        /*
        NumParticles.push_back(round((MassFractionOfFillers[i]/100.0)*Vrve*(1.0/CompRho)/CNT_data[i].GetAvgParticleWeight()));
        CB_data[i].SetNumParticles(NumParticles[i]);
        
         //Someday we will add CNT to this simulation!!!
        */
    }
}

void RVEdata::Read_Agg_Files()
{
    std::string cwd = getcwd(NULL,0);
    std::string dir = cwd + "/Agg_Files/";
    boost::filesystem::path p(dir);
    std::vector<std::string> file_names;
    
    for (auto i = boost::filesystem::directory_iterator(p); i != boost::filesystem::directory_iterator(); i++)
    {
        //we eliminate directories in a list
        if (boost::filesystem::is_directory(i->path()))
        {
            continue;
        }
        else
        {
            if (i->path().filename().string()[0] == '.')
            {
                continue;
            }
            if (i->path().filename().string()[0] != '_')
            {
                file_names.push_back(i->path().filename().string());
            }
        }
    }
    for (int i = 0; i < file_names.size(); i++)
    {
        //std::cout << "file no. " << i << "\r" << std::flush;
        std::ifstream Agg_file;
        Agg_file.open(dir + file_names[i]);
        std::vector<std::vector<int> > temp_agg_xyz;
        std::vector<std::vector<int> > translated_agg_xyz;
        int num;
        Agg_file >> num;
        std::string x="", y="", z="", r="";
        int xi=0, yi=0, zi=0, ri=0;
        int minx=100, miny=100, minz=100, maxx=-100, maxy=-100, maxz=-100;
        std::string line;
        //Ignore second two lines
        std::getline(Agg_file,line);
        std::getline(Agg_file,line);
        while (std::getline(Agg_file,line))
        {
            //std::cout << line << std::endl;
            std::istringstream s(line);
            s >> x >> std::ws >> y >> std::ws >> z >> std::ws >> r;
            xi = std::stoi(x);
            yi = std::stoi(y);
            zi = std::stoi(z);
            ri = std::stoi(r);
            if (xi > maxx){ maxx = xi;}
            if (yi > maxy){ maxy = yi;}
            if (zi > maxz){ maxz = zi;}
            if (xi < minx){ minx = xi;}
            if (yi < miny){ miny = yi;}
            if (zi < minz){ minz = zi;}
            temp_agg_xyz.push_back(std::vector<int> {xi, yi, zi});
        }
        //std::cout << minx << " " << miny << " " << minz << std::endl;
        if (i == 0){no_p_per_agg = temp_agg_xyz.size();}
        for (int j = 0; j < temp_agg_xyz.size(); j++)
        {
            translated_agg_xyz.push_back(std::vector<int> {temp_agg_xyz[j][0]-minx, temp_agg_xyz[j][1]-miny, temp_agg_xyz[j][2]-minz});
            //std::cout << "j: " << j << " " << temp_agg_xyz[j][0]-minx << " " << temp_agg_xyz[j][1]-miny << " " << temp_agg_xyz[j][2]-minz << std::endl;
        }
        //std::cout << "minx = " << minx << " miny = " << miny <<  " minz = " << minz << " maxx = " << maxx << " maxy = " << maxy << " maxz = " << maxz << std::endl;
        Agglomerates.push_back(translated_agg_xyz);
        std::vector<long int> temp_indices;
        for (int j = 0; j < translated_agg_xyz.size(); j++)
        {
            long int index = ((Ybound+1)*(Zbound+1)*translated_agg_xyz[j][0]) + ((Ybound+1)*translated_agg_xyz[j][2]) + translated_agg_xyz[j][1];
            if (index < 0)
            {
                std::cout << "index<0: " << index << ", xyz - " << translated_agg_xyz[j][0] << "," << translated_agg_xyz[j][1] << "," << translated_agg_xyz[j][2] << std::endl;
            }
            temp_indices.push_back(index);
        }
        Agglomerates_Indices.push_back(temp_indices);
        Agglomerates_Bounds.push_back(std::vector<int> {maxx-minx, maxy-miny, maxz-minz});
    }
}

int RVEdata::GetNumParticles()
{
    return particle_data.size();
}

int RVEdata::GetNumCBFillers()
{
    return NumCBFillers;
}

void RVEdata::Build_RVE()
{
    //Determine number of available processors for threading
    int NumProcessors = std::thread::hardware_concurrency()*2;
	
    //Determine number of random particle seeds
    long int num_rand = ceil(Total_Num_Particles*Pct_Rand);
    
    if (num_rand == 0)
    {
        std::cout << "ZERO RAND ERR: Pct_Rand = " << Pct_Rand << " TotalNum = " << Total_Num_Particles << std::endl;
        exit(1);
    }

    //Initialize the placement mapping variables
    //std::cout << "Shell/Potential Maps Initialized" << std::endl;
    Initialize_CB_Potential_Shell_Maps();
    Initialize_CB_Strain_Potential_Maps();
    //std::cout << "Strain Potential Maps Initialized" << std::endl;
    
    //Place the seed particles
    //Place_Seeds();
    //std::cout << "Placing Rand Agg" << std::endl;
    //Place the rest of the particles
    Place_Rand_Agglomerates(num_rand);
    //std::cout << "Placing Forced Agg" << std::endl;
    Place_Forced_Agglomerates(num_rand);
    //std::cout << "Placing Remainder" << std::endl;
    Place_Remainder(num_rand);

    /*
    //Output potential and shell for checking. Files get really 'uge really quicklike
    std::ofstream check;
    check.open("check.csv");
    check << Total_Num_Particles << std::endl << std::endl;
    for (int i = 0; i < Total_Num_Particles; i++)
    {
        check << GetCBParticle(i)->GetXCoordinate() << " " << GetCBParticle(i)->GetYCoordinate() << " " << GetCBParticle(i)->GetZCoordinate() << " " << 15 << std::endl;
    }
    check.close();*/
    
    //Output coordinates for checking and comparing the things to the things
    /*
    std::ofstream check;
    check.open("check1.csv");
    check << num_rand + 1 + Potential.count() << std::endl << std::endl;
    for (int i = 0; i < num_rand; i++)
    {
        check << GetCBParticle(i)->GetXCoordinate() << " " << GetCBParticle(i)->GetYCoordinate() << " " << GetCBParticle(i)->GetZCoordinate() << " " << 15 << " " << 0 << std::endl;
    }
    
    check << GetCBParticle(num_rand)->GetXCoordinate() << " " << GetCBParticle(num_rand)->GetYCoordinate() << " " << GetCBParticle(num_rand)->GetZCoordinate() << " " << 15 << " " << 100 << std::endl;
    
    long int counter = 0;
    for (int i = 0; i <= Xbound; i++)
    {
        for (int k = 0; k <= Zbound; k++)
        {
            for (int j = 0; j <= Ybound; j++)
            {
                if (Potential[counter] == 1)
                {
                    check << i << " " << j << " " << k << " " << 0.5 << " " << 0 << std::endl;
                }
                counter++;
            }
        }
    }
    check.close();*/
}

void RVEdata::Place_Rand_Agglomerates(long int num_rand)
{
    //Initialize random number generator instance
    
    //********** Modified Input for Pct Rand and thickness **************
    //Initialize random number generator instance for determining params
    // Convert both the integers to string
    std::string s1 = std::to_string(time(NULL));
    std::string s2 = std::to_string(getpid());
    // Concatenate both strings
    std::string s = s1 + s2;
    // Convert the concatenated string
    // to integer
    long int c = std::stol(s);
    srand(c);
    std::mt19937_64 gen(rand());
    std::uniform_int_distribution<long int> XYZ(0, Potential.size());
    std::uniform_int_distribution<int> Num_Agglomerates(0,Agglomerates.size()-1);
    std::vector<long int> Checked;
    
    long int next_index;
    int next_agglomerate;
    
    num_rand_agg = floor((Total_Num_Particles*Pct_Rand)/no_p_per_agg);
    if (num_rand_agg == 0)
    {
        num_rand_agg = 1;
    }
    //std::cout << "num_rand_agg: " << num_rand_agg << std::endl;
    
    running_index = 0;
    
    for (long int i = 0; i < num_rand_agg; i ++)
    {
        //std::cout << "Agglomerate #" << i << "\r" << std::flush;
        bool Good = false;
        int total_iter = 0;
        //Pick aggregate from the list
        next_agglomerate = Num_Agglomerates(gen);
        
        while (Good == false)
        {
            bool Already_Checked = true;
            
            while(Already_Checked == true)
            {
            //Set next index to random location
                next_index = XYZ(gen);
                //Check the bounds of the agg to make sure it won't go beyond the RVE
                int x, y, z, x_max, y_max, z_max;
                x = floor(next_index/((Ybound+1)*(Zbound+1)));
                z = floor((next_index-(x*(Ybound+1)*(Zbound+1)))/(Ybound+1));
                y = next_index - (x*(Ybound+1)*(Zbound+1)) - (z*(Ybound+1));
                x_max = x + Agglomerates_Bounds[next_agglomerate][0];
                y_max = y + Agglomerates_Bounds[next_agglomerate][1];
                z_max = z + Agglomerates_Bounds[next_agglomerate][2];
                if (x_max > Xbound || y_max > Ybound || z_max > Zbound)
                {
                    Checked.push_back(next_index);
                    continue;
                }
                
                for (int k = 0; k < Checked.size(); k++)
                {
                    if (next_index == Checked[k]){ break;}
                    if (k == Checked.size()-1){ Already_Checked = false;}
                }
            }
            
            //Loop through the agglomerate indices and check to see if possible
            for (int j = 0; j < Agglomerates[next_agglomerate].size(); j++)
            {
                if (Shell[next_index + Agglomerates_Indices[next_agglomerate][j]] == 1)
                {
                    Checked.push_back(next_index);
                    break;
                }
                if (j == Agglomerates_Indices[next_agglomerate].size()-1)
                {
                    Good = true;
                }
            }
            if (Good == true)
            {
                Checked.clear();
                for (int z = 0; z < Agglomerates_Indices[next_agglomerate].size(); z++)
                {
                    //Got the next index to map. Getting its coordinates and mapping it:
                    Initialize_Index(next_index + Agglomerates_Indices[next_agglomerate][z], running_index);
                    
                    //Map the particle
                    Map_CB(running_index);
                    
                    //Increment running index
                    running_index += 1;
                }
            }
        }
        //Find neighbors?
        //Find_Neighbors(i);
    }
    //std::cout << std::endl;
}

void RVEdata::Place_Forced_Agglomerates(long int num_rand)
{
    //Initialize random number generator instance
    
    //********** Modified Input for Pct Rand and thickness **************
    //Initialize random number generator instance for determining params
    // Convert both the integers to string
    std::string s1 = std::to_string(time(NULL));
    std::string s2 = std::to_string(getpid());
    // Concatenate both strings
    std::string s = s1 + s2;
    // Convert the concatenated string
    // to integer
    long int c = std::stol(s);
    srand(c);
    std::mt19937_64 gen(rand());
    std::uniform_int_distribution<long int> XYZ(0, Potential.size());
    std::uniform_int_distribution<int> Num_Agglomerates(0, Agglomerates.size()-1);
    std::vector<long int> Checked;
    
    long int next_index;
    int next_agglomerate;
    
    num_forced_agg = floor(((Total_Num_Particles - (floor((Total_Num_Particles*Pct_Rand)/no_p_per_agg)*no_p_per_agg))*(1.-Pct_Fill))/no_p_per_agg);    //std::cout << "num_forced_agg: " << num_forced_agg << std::endl;
        
    for (long int i = 0; i < num_forced_agg; i ++)
    {
        //std::cout << "Agglomerate #" << i << std::endl;
        bool Good = false;
        int total_iter = 0;
        //Pick aggregate from the list
        next_agglomerate = Num_Agglomerates(gen);
        
        int num_attempts = 0;
        while (Good == false)
        {
            //std::cout << "Attempt # = " << num_attempts << "\r" << std::flush;
            bool Already_Checked = true;
            
            while(Already_Checked == true)
            {
                //Set next index to random location
                next_index = XYZ(gen);
                //Check the bounds of the agg to make sure it won't go beyond the RVE
                int x, y, z, x_max, y_max, z_max;
                x = floor(next_index/((Ybound+1)*(Zbound+1)));
                z = floor((next_index-(x*(Ybound+1)*(Zbound+1)))/(Ybound+1));
                y = next_index - (x*(Ybound+1)*(Zbound+1)) - (z*(Ybound+1));
                x_max = x + Agglomerates_Bounds[next_agglomerate][0];
                y_max = y + Agglomerates_Bounds[next_agglomerate][1];
                z_max = z + Agglomerates_Bounds[next_agglomerate][2];
                if (x_max > Xbound || y_max > Ybound || z_max > Zbound)
                {
                    Checked.push_back(next_index);
                    continue;
                }
                
                for (int k = 0; k < Checked.size(); k++)
                {
                    if (next_index == Checked[k]){ break;}
                    if (k == Checked.size()-1){ Already_Checked = false;}
                }
            }
            
            bool Got_Neighbors = false;
            //Loop through the agglomerate indices and check to see if possible
            for (int j = 0; j < Agglomerates[next_agglomerate].size(); j++)
            {
                if (Shell[next_index + Agglomerates_Indices[next_agglomerate][j]] == 1)
                {
                    Checked.push_back(next_index);
                    break;
                }
                if (Potential[next_index + Agglomerates_Indices[next_agglomerate][j]] == 1)
                {
                    //std::cout << "ever?" << std::endl;
                    Got_Neighbors = true;
                }
                if (j == Agglomerates_Indices[next_agglomerate].size()-1 && Got_Neighbors == true)
                {
                    Good = true;
                }
            }

            if (Got_Neighbors == true && Good == true)
            {
                Checked.clear();
                for (int z = 0; z < Agglomerates_Indices[next_agglomerate].size(); z++)
                {
                    //Got the next index to map. Getting its coordinates and mapping it:
                    Initialize_Index(next_index + Agglomerates_Indices[next_agglomerate][z], running_index);
                    
                    //Map the particle
                    Map_CB(running_index);
                    
                    //Increment running index
                    running_index += 1;
                }
            }
            num_attempts++;
        }
        //Find neighbors?
        //Find_Neighbors(i);
    }
    //std::cout << std::endl;
}

void RVEdata::Place_Remainder(long int num_rand)
{
    //Initialize random number generator instance
    
    //********** Modified Input for Pct Rand and thickness **************
    //Initialize random number generator instance for determining params
    // Convert both the integers to string
    std::string s1 = std::to_string(time(NULL));
    std::string s2 = std::to_string(getpid());
    // Concatenate both strings
    std::string s = s1 + s2;
    // Convert the concatenated string
    // to integer
    long int c = std::stol(s);
    srand(c);
    std::mt19937_64 gen(rand());
    std::uniform_int_distribution<long int> XYZ(0, Potential.size());
    
    //int num_rand_agg = floor((Total_Num_Particles*Pct_Rand)/100.);
    //int num_forced_agg = floor(((Total_Num_Particles - (floor((Total_Num_Particles*Pct_Rand)/100.)*100.))*(1.-Pct_Fill))/100.);
    int num_remaining = Total_Num_Particles - ((num_rand_agg+num_forced_agg)*no_p_per_agg);
    long int next_index;
    
    for (long int i = Total_Num_Particles-num_remaining; i < Total_Num_Particles; i ++)
    {
        //std::cout << "Particle #" << i << "\r" << std::flush;
        bool Good = false;
        int total_iter = 0;
        while (Good == false)
        {
            //Set next index to random location
            next_index = XYZ(gen);
            if (Potential[next_index] == 1)
            {
                Good = true;
            }
            total_iter ++;
            if (total_iter > 10000000)
            {
                std::cout << "Max Iter: " << total_iter << std::endl;
                long int sum_ones = 0;
                for (long int z=0; z < Potential.size(); z++)
                {
                    if (Potential[z] ==1)
                    {
                        sum_ones += 1;
                    }
                }
                std::cout << "sum_ones in Potential = " << sum_ones << std::endl;
                exit(1);
            }
        }
        //Got the next index to map. Getting its coordinates and mapping it:
        Initialize_Index(next_index, i);
        
        //Map the particle
        Map_CB(i);
        
        //Find neighbors?
        //Find_Neighbors(i);
    }
    //std::cout << std::endl;
}

void RVEdata::Get_Next_Location(long int &next)
{
    //Check if "next" is a valid placement
    if (Potential[next] == 1) { return;}
    
    //Check up to last index and return once a valid placement is found
    for (long int i = next; i < Potential.size(); i++)
    {
        if (Potential[i] == 1) { next = i; return;}
    }
    
    //If arrived here, need to start over at first (0th) index
    for (long int i = 0; i < next; i++)
    {
        if (Potential[i] == 1) { next = i; return;}
    }
    
    //If got all the way here, either the thing is full or something's wrong
    std::cout << "RVE FULL or WORSE... THERE'S A BUG!\n";
    exit(1);
}

void RVEdata::Initialize_Index(long int next_index, long int ID)
{
    int x, y, z;
    x = floor(next_index/((Ybound+1)*(Zbound+1)));
    z = floor((next_index-(x*(Ybound+1)*(Zbound+1)))/(Ybound+1));
    y = next_index - (x*(Ybound+1)*(Zbound+1)) - (z*(Ybound+1));
    
    //long int CB_index = ((Ybound+1)*(Zbound+1)*x) + ((Ybound+1)*z) + y;
    
    GetCBParticle(GetParticle(ID)->GetTypeID())->SetXYZ(x, y, z, Xbound, Td);
}

void RVEdata::Place_Seeds()
{
    int NumProcessors = std::thread::hardware_concurrency()*2;
	
    //Start by placing the seed particles
    long int num_rand = ceil(Total_Num_Particles*Pct_Rand);
    //std::cout << "num_rand = " << num_rand << std::endl;
	
    //Initialize random number generator instance
    //********** Modified Input for Pct Rand and thickness **************
    //Initialize random number generator instance for determining params
    // Convert both the integers to string
    std::string s1 = std::to_string(time(NULL));
    std::string s2 = std::to_string(getpid());
    // Concatenate both strings
    std::string s = s1 + s2;
    // Convert the concatenated string
    // to integer
    long int c = std::stol(s);
    srand(c);
    std::mt19937_64 gen(rand());
    
    //Define uniform standard distribution of potential coordinates
    std::uniform_int_distribution<long int> X_coords(0, Xbound);
    std::uniform_int_distribution<long int> Y_coords(0, Ybound);
    std::uniform_int_distribution<long int> Z_coords(0, Zbound);
    
    //std::cout << Xbound << " " << Ybound << " " << Zbound << std::endl;
    
    //Defining seed placement variables
    bool Is_Placed = false;
    bool Valid_Place = true;
    int iPSoftness, jPSoftness, dpi, dpj;
    
    //Temporary Neighbor data
    std::vector<long int> neighIndex;
    std::vector<double> neighDist;

    //std::cout << "Num Rand: " << num_rand << std::endl;
    //Loop through and place all seed particles while checking placement fidelity
	for (long int i = 0; i < num_rand; i++)
	{
        //std::cout << "Particle #" << i << "\r" << std::flush;
        Is_Placed = false;
        while (Is_Placed == false)
        {
            //Generate random coordinate for seed particle
            long int x = X_coords(gen);
            long int y = Y_coords(gen);
            long int z = Z_coords(gen);
            
            if (Seed_Check_Adjacent_Indices(x,y,z,i) == true)
            {
                Is_Placed = true;
                GetCBParticle(GetParticle(i)->GetTypeID())->SetXYZ(x, y, z, Xbound, Td);
                Map_CB(i);
            }
        }
	}
}

void RVEdata::Map_CB(long int CB)
{
    //Get CB Particle Diameter
    int diam = GetCBParticle(CB)->GetDiameter();
    int x = GetCBParticle(CB)->GetXCoordinate();
    int y = GetCBParticle(CB)->GetYCoordinate();
    int z = GetCBParticle(CB)->GetZCoordinate();
    
    //Add index to the deque of indices for future placement checking
    All_Indices[x].coeffRef(y,z) = CB;
    
    long int CB_index = ((Ybound+1)*(Zbound+1)*x) + ((Ybound+1)*z) + y;
    
    GetCBParticle(CB)->Set_Index(CB_index);
    
    //Loop through the lists of potential and shells to add/remove to/from the potenial and shell bitfields
    int which_set = Get_CB_Map_Index(CB);
    long int max_index = Potential.size();
    long int which_index;
    
    for (int i = 0; i < CB_Potential_Maps_Indices[which_set].size(); i++)
    {
        which_index = CB_index + CB_Potential_Maps_Indices[which_set][i];
        if (x + CB_Potential_Maps_X[which_set][i] > Xbound || x + CB_Potential_Maps_X[which_set][i] < 0) {continue;}
        if (y + CB_Potential_Maps_Y[which_set][i] > Ybound || y + CB_Potential_Maps_Y[which_set][i] < 0) {continue;}
        if (z + CB_Potential_Maps_Z[which_set][i] > Zbound || z + CB_Potential_Maps_Z[which_set][i] < 0) {continue;}
        if (which_index > max_index) {continue;}
        if (Shell[which_index] == 1) {continue;}
        
        //if made it past the "if gauntlet" then the placement of a 1 in the potential column is the play
        Potential[which_index] = 1;
    }
    for (int i = 0; i < CB_Shell_Maps_Indices[which_set].size(); i++)
    {
        which_index = CB_index + CB_Shell_Maps_Indices[which_set][i];
        if (x + CB_Shell_Maps_X[which_set][i] > Xbound || x + CB_Shell_Maps_X[which_set][i] < 0) {continue;}
        if (y + CB_Shell_Maps_Y[which_set][i] > Ybound || y + CB_Shell_Maps_Y[which_set][i] < 0) {continue;}
        if (z + CB_Shell_Maps_Z[which_set][i] > Zbound || z + CB_Shell_Maps_Z[which_set][i] < 0) {continue;}
        if (which_index > max_index) {continue;}
        
        //if made it past the "if gauntlet" then the placement of a 1 in the shell column is the play
        Shell[which_index] = 1;
        Potential[which_index] = 0;
    }
}

int RVEdata::Get_CB_Map_Index(long int ID)
{
    //return the index to the correct map for the particular particle corresponding to its size
    int type = CB_particles[particle_data[ID].GetTypeID()].GetMatType();
    int diam = CB_particles[particle_data[ID].GetTypeID()].GetDiameter();
    
    int count = 0;
    for (int i = 0; i < type; i++)
    {
        count += CB_data[i].GetMaxDiameter() - CB_data[i].GetMinDiameter();
        //std::cout << "count = " << count << std::endl;
    }
    count += diam - CB_data[type].GetMinDiameter();
    //std::cout << "count final = " << count << std::endl;
    //std::cout << "diam = " << diam << ", MinD = " << CB_data[type].GetMinDiameter() << std::endl;
    
    return count;
}

int RVEdata::Get_Max_CB_Map_Index(long int ID)
{
    //return the index to the correct map for the particular particle corresponding to its size
    int type = CB_particles[particle_data[ID].GetTypeID()].GetMatType();
    int diam = CB_data[type].GetMaxDiameter();
    
    int count = 0;
    for (int i = 0; i < type; i++)
    {
        count += CB_data[i].GetMaxDiameter() - CB_data[i].GetMinDiameter();
        //std::cout << "count = " << count << std::endl;
    }
    count += diam - CB_data[type].GetMinDiameter();
    //std::cout << "count final = " << count << std::endl;
    //std::cout << "diam = " << diam << ", MinD = " << CB_data[type].GetMinDiameter() << std::endl;
    
    return count;
}

void RVEdata::Initialize_CB_Potential_Shell_Maps()
{
    //Loop through each filler and add the mapping scheme for each particle size for each type
    int Map_Count = 0;
    for (int x = 0; x < NumCBFillers; x++)
    {
        for (int y = 0; y <= GetCBMatType(x)->GetMaxDiameter() - GetCBMatType(x)->GetMinDiameter(); y++)
        {
            double diam = (GetCBMatType(x)->GetMinDiameter() + y)*mult;
            
            int PS = GetCBMatType(x)->GetPenetrationAllowance()*mult;

            int xb = diam*2 ;   //+ (2*Td);
            int yb = diam*2 ;   //+ (2*Td);
            int zb = diam*2 ;   //+ (2*Td);
            
            CB_Potential_Maps_X.push_back(std::vector<int>(0));
            CB_Potential_Maps_Y.push_back(std::vector<int>(0));
            CB_Potential_Maps_Z.push_back(std::vector<int>(0));
            CB_Shell_Maps_X.push_back(std::vector<int>(0));
            CB_Shell_Maps_Y.push_back(std::vector<int>(0));
            CB_Shell_Maps_Z.push_back(std::vector<int>(0));
            CB_Potential_Maps_Indices.push_back(std::vector<long int>(0));
            CB_Shell_Maps_Indices.push_back(std::vector<long int>(0));

            int x0 = diam;  //+Td,
            int y0 = diam;  //+Td,
            int z0 = diam;  //+Td;
            double dist;
            
            long int index_0 = 0 - ((Ybound+1)*(Zbound+1)*x0) - ((Ybound+1)*z0) - (y0);
            long int zer = 0 - index_0;
            
            for (int i = 0; i <= xb; i++)
            {
                for (int j = 0; j <= zb; j++)
                {
                    for (int k = 0; k <= yb; k++)
                    {
                        dist = mult*pow(((abs(y0-k)*abs(y0-k)) + (abs(z0-j)*abs(z0-j)) + (abs(x0-i)*abs(x0-i))), 0.5);
                        if (dist < (diam + Td) || std::fabs(dist - diam - Td) < 0.000001)
                        {
                            if (dist > (diam-PS) || std::fabs(dist-diam+PS) < 0.000001)
                            {
                                CB_Potential_Maps_X[Map_Count].push_back(i-x0);
                                CB_Potential_Maps_Z[Map_Count].push_back(j-z0);
                                CB_Potential_Maps_Y[Map_Count].push_back(k-y0);
                                CB_Potential_Maps_Indices[Map_Count].push_back(index_0);
                            }
                            else
                            {
                                CB_Shell_Maps_X[Map_Count].push_back(i-x0);
                                CB_Shell_Maps_Z[Map_Count].push_back(j-z0);
                                CB_Shell_Maps_Y[Map_Count].push_back(k-y0);
                                CB_Shell_Maps_Indices[Map_Count].push_back(index_0);
                            }
                        }
                        /*else     //Too far outside the bounds of the particle to reach
                        {
                            CB_Shell_Maps[y].push_back(index_0);
                        }*/
                        index_0++;
                    }
                    index_0 += ((Ybound+1) - yb - 1);
                }
                index_0 += ((Ybound+1)*(Zbound+1) - ((zb+1)*(Ybound+1)));
            }
            //***********************************************************************************
             /******************* FOR IMAGE OUTPUT/DEBUGGING OF MAP ****************************
             //**********************************************************************************
            std::ofstream output;
            std::cout << "Map_Count = " << Map_Count << std::endl;
            output.open("map_" + std::to_string(Map_Count) + ".csv");
            output << CB_Potential_Maps_Indices[Map_Count].size() + CB_Shell_Maps_Indices[Map_Count].size() << std::endl << std::endl;
            for (int p = 0; p < CB_Potential_Maps_Indices[Map_Count].size(); p++)
            {
                output << x0+CB_Potential_Maps_X[Map_Count][p] << " " << y0+CB_Potential_Maps_Y[Map_Count][p] << " " << z0+CB_Potential_Maps_Z[Map_Count][p] << " 0.5 100 0 " << CB_Potential_Maps_Indices[Map_Count][p] + zer <<std::endl;
            }
            
            for (int p = 0; p < CB_Shell_Maps_Indices[Map_Count].size(); p++)
            {
                output << x0+CB_Shell_Maps_X[Map_Count][p] << " " << y0+CB_Shell_Maps_Y[Map_Count][p] << " " << z0+CB_Shell_Maps_Z[Map_Count][p] << " 0.5 0 100 " << CB_Shell_Maps_Indices[Map_Count][p] << std::endl;
            }
            
            output.close();
             //**************************************************************************
            */
            Map_Count++;
        }
    }
}

void RVEdata::Initialize_CB_Strain_Potential_Maps()
{
    //Loop through each filler and add the mapping scheme for each particle size for each type
    int Map_Count = 0;
    for (int x = 0; x < NumCBFillers; x++)
    {
        for (int y = 0; y <= GetCBMatType(x)->GetMaxDiameter() - GetCBMatType(x)->GetMinDiameter(); y++)
        {
            double diam = (GetCBMatType(x)->GetMinDiameter() + y)*mult;
            
            int PS = GetCBMatType(x)->GetPenetrationAllowance()*mult;
            
            int max_strain_dist = 6;          //must be an even integer
            int xb = diam*2 + max_strain_dist + (2*Td);
            int yb = diam*2 + max_strain_dist + (2*Td);
            int zb = diam*2 + max_strain_dist + (2*Td);
            
            CB_Strain_Potential_Maps_X.push_back(std::vector<int>(0));
            CB_Strain_Potential_Maps_Y.push_back(std::vector<int>(0));
            CB_Strain_Potential_Maps_Z.push_back(std::vector<int>(0));
            CB_Strain_Potential_Maps_Indices.push_back(std::vector<long int>(0));
            
            CB_Strain_Shell_Maps_X.push_back(std::vector<int>(0));
            CB_Strain_Shell_Maps_Y.push_back(std::vector<int>(0));
            CB_Strain_Shell_Maps_Z.push_back(std::vector<int>(0));
            CB_Strain_Shell_Maps_Indices.push_back(std::vector<long int>(0));
            
            int x0 = xb/2;  //+Td,
            int y0 = xb/2;  //+Td,
            int z0 = xb/2;  //+Td;
            double dist;
            
            long int index_0 = 0 - ((Ybound+1)*(Zbound+1)*x0) - ((Ybound+1)*z0) - (y0);
            long int zer = 0 - index_0;
            
            for (int i = 0; i <= xb; i++)
            {
                for (int j = 0; j <= zb; j++)
                {
                    for (int k = 0; k <= yb; k++)
                    {
                        dist = mult*pow(((abs(y0-k)*abs(y0-k)) + (abs(z0-j)*abs(z0-j)) + (abs(x0-i)*abs(x0-i))), 0.5);
                        if (dist < (diam + Td + max_strain_dist) || std::fabs(dist - diam - Td - max_strain_dist) < 0.000001)
                        {
                            if (dist > (diam-PS) || std::fabs(dist-diam+PS) < 0.000001)
                            {
                                CB_Strain_Potential_Maps_X[Map_Count].push_back(i-x0);
                                CB_Strain_Potential_Maps_Z[Map_Count].push_back(j-z0);
                                CB_Strain_Potential_Maps_Y[Map_Count].push_back(k-y0);
                                CB_Strain_Potential_Maps_Indices[Map_Count].push_back(index_0);
                            }
                            else
                            {
                                continue;
                            }
                        }
                        else     //Too far outside the bounds of the particle to reach
                        {
                            /*
                             CB_Strain_Shell_Maps_X[Map_Count].push_back(i-x0);
                             CB_Strain_Shell_Maps_Z[Map_Count].push_back(j-z0);
                             CB_Strain_Shell_Maps_Y[Map_Count].push_back(k-y0);
                             CB_Strain_Shell_Maps_Indices[Map_Count].push_back(index_0);
                             */
                        }
                        index_0++;
                    }
                    index_0 += ((Ybound+1) - yb - 1);
                }
                index_0 += ((Ybound+1)*(Zbound+1) - ((zb+1)*(Ybound+1)));
            }
            //***********************************************************************************
            /* ******************* FOR IMAGE OUTPUT/DEBUGGING OF MAP ****************************
             //*********************************************************************************
            std::ofstream output;
            std::cout << "Map_Count = " << Map_Count << std::endl;
            output.open("strain-map_" + std::to_string(Map_Count) + ".csv");
            output << CB_Strain_Potential_Maps_Indices[Map_Count].size()+CB_Strain_Shell_Maps_Indices[Map_Count].size() << std::endl << std::endl;
            for (int p = 0; p < CB_Strain_Potential_Maps_Indices[Map_Count].size(); p++)
            {
                output << x0+CB_Strain_Potential_Maps_X[Map_Count][p] << " " << y0+CB_Strain_Potential_Maps_Y[Map_Count][p] << " " << z0+CB_Strain_Potential_Maps_Z[Map_Count][p] << " 0.5 100 0 " << CB_Strain_Potential_Maps_Indices[Map_Count][p] + zer <<std::endl;
            }
            
            for (int p = 0; p < CB_Strain_Shell_Maps_Indices[Map_Count].size(); p++)
            {
                output << x0+CB_Strain_Shell_Maps_X[Map_Count][p] << " " << y0+CB_Strain_Shell_Maps_Y[Map_Count][p] << " " << z0+CB_Strain_Shell_Maps_Z[Map_Count][p] << " 0.5 0 100 " << CB_Strain_Shell_Maps_Indices[Map_Count][p] + zer <<std::endl;
            }
            
            output.close();
            Map_Count++;
            */
        }
    }
}

void RVEdata::Find_Neighbors()
{
    //Count the # of processors for parallelization
    int num_processors = std::thread::hardware_concurrency();
    
    //Started at the bottom now we're here. Make the futures for parralelizatiomn.
    std::vector<std::future<void> > Search;
    long int num_rand = ceil(Total_Num_Particles*Pct_Rand);

    long int index = 0;
    int Active_Threads = 0;
    
    while (index < Total_Num_Particles)
    {
        if (Search.size() < num_processors)
        {
            Search.push_back(std::async(std::launch::async, &RVEdata::Check_Adjacent_Indices, this, index));
            index++;
            //std::cout << "Neighbors: " << index << "\r";
        }
        else
        {
            for (int i = 0; i < num_processors; i++)
            {
                if (index == Total_Num_Particles) { break;}
                auto status = Search[i].wait_for(std::chrono::nanoseconds(10));
                if (status == std::future_status::ready)
                {
                    Search[i] = std::async(std::launch::async, &RVEdata::Check_Adjacent_Indices, this, index);
                    index++;
                    //std::cout << "Neighbors: " << index << "\r";
                }
            }
        }
    }
    //std::cout << std::endl;
    
    //Wait for all the threads to finish up!
    for (int i = 0; i < num_processors; i++)
    {
        Search[i].wait();
    }
}

void RVEdata::Find_Strained_Neighbors(float strain_elong)
{
    //Count the # of processors for parallelization
    int num_processors = std::thread::hardware_concurrency();
    
    float Poissons = p_data.GetPoissonsRatio();
    float eps1 = 1 + strain_elong;
    float eps2, eps3;
    if (Poissons != 0)
    {
        eps2 = 1 - (strain_elong*Poissons);
        eps3 = 1 - (strain_elong*Poissons);
    }
    else
    {
        eps2 = 1;
        eps3 = 1;
    }
    
    //Started at the bottom now we're here. Make the futures for parralelizatiomn.
    std::vector<std::future<void> > Search;
    long int num_rand = ceil(Total_Num_Particles*Pct_Rand);
    
    long int index = 0;
    int Active_Threads = 0;
    
    while (index < Total_Num_Particles)
    {
        if (Search.size() < num_processors)
        {
            Search.push_back(std::async(std::launch::async, &RVEdata::Check_Adjacent_Strain_Indices, this, index, eps1, eps2, eps3));
            index++;
            //std::cout << "Neighbors: " << index << "\r";
        }
        else
        {
            for (int i = 0; i < num_processors; i++)
            {
                if (index == Total_Num_Particles) { break;}
                auto status = Search[i].wait_for(std::chrono::nanoseconds(10));
                if (status == std::future_status::ready)
                {
                    Search[i] = std::async(std::launch::async, &RVEdata::Check_Adjacent_Strain_Indices, this, index, eps1, eps2, eps3);
                    index++;
                    //std::cout << "Neighbors: " << index << "\r";
                }
            }
        }
    }
    //std::cout << std::endl;
    
    //Wait for all the threads to finish up!
    for (int i = 0; i < num_processors; i++)
    {
        Search[i].wait();
    }
}

void RVEdata::Reset_Neighbors()
{
    for (long int i = 0; i < Total_Num_Particles; i++)
    {
        GetCBParticle(GetParticle(i)->GetTypeID())->Clear_Neighbors();
    }
}

void RVEdata::Check_Adjacent_Strain_Indices(long int i, float eps1, float eps2, float eps3)
{
    //Compare the distance of adjacent indices. Loop through the set of potential neighbor spots
    long int x = GetCBParticle(GetParticle(i)->GetTypeID())->GetXCoordinate();
    long int y = GetCBParticle(GetParticle(i)->GetTypeID())->GetYCoordinate();
    long int z = GetCBParticle(GetParticle(i)->GetTypeID())->GetZCoordinate();
    
    long int CB_index = ((Ybound+1)*(Zbound+1)*x) + ((Ybound+1)*z) + y;
    
    int iPSoftness, jPSoftness, dpi, dpj;
    double dist;
    
    int which_set = Get_Max_CB_Map_Index(i);
    long int which_index;
    long int max_index = Potential.size();
    
    for (int j = 0; j < CB_Strain_Potential_Maps_Indices[which_set].size(); j++)
    {
        //if (j == i) { continue;}
        long int p_ID;
        which_index = CB_index + CB_Strain_Potential_Maps_Indices[which_set][j];
        
        //Check the bounds of the potential map index
        if (x + CB_Strain_Potential_Maps_X[which_set][j] > Xbound || x + CB_Strain_Potential_Maps_X[which_set][j] < 0) {continue;}
        if (y + CB_Strain_Potential_Maps_Y[which_set][j] > Ybound || y + CB_Strain_Potential_Maps_Y[which_set][j] < 0) {continue;}
        if (z + CB_Strain_Potential_Maps_Z[which_set][j] > Zbound || z + CB_Strain_Potential_Maps_Z[which_set][j] < 0) {continue;}
        if (which_index > max_index) {continue;}
        //If passed the "if" statements, the index is within the RVE. Check it for a particle
        
        if(isNull(All_Indices[x+CB_Strain_Potential_Maps_X[which_set][j]], y+CB_Strain_Potential_Maps_Y[which_set][j], z+CB_Strain_Potential_Maps_Z[which_set][j]) == true)
        {
            //If == 0, then there is no particle at the index. Onto the next one!
            continue;
        }
        else
        {
            p_ID = All_Indices[x+CB_Strain_Potential_Maps_X[which_set][j]].coeffRef(y+CB_Strain_Potential_Maps_Y[which_set][j], z+CB_Strain_Potential_Maps_Z[which_set][j]);
        }
        
        int Particle_i_Type = GetParticle(i)->GetType();
        int Particle_j_Type = GetParticle(p_ID)->GetType();
        
        //Set particle softness of particle j
        if (Particle_i_Type <= NumCBFillers)
        {
            iPSoftness = GetCBMatType(GetCBParticle(GetParticle(i)->GetTypeID())->GetMatType())->GetPenetrationAllowance()*mult;
            dpi = GetCBParticle(GetParticle(i)->GetTypeID())->GetDiameter()*mult;
        }
        
        if (Particle_j_Type <= NumCBFillers)
        {
            jPSoftness = GetCBMatType(GetCBParticle(GetParticle(p_ID)->GetTypeID())->GetMatType())->GetPenetrationAllowance()*mult;
            dpj = GetCBParticle(GetParticle(p_ID)->GetTypeID())->GetDiameter()*mult;
        }
        
        //Define the most restrictive penetration allowance
        double PS;
        if (iPSoftness < jPSoftness)
        {
            PS = iPSoftness;
        }
        else
        {
            PS = jPSoftness;
        }
        
        dist = GetCBParticle(GetParticle(i)->GetTypeID())->FindStrainDistanceToNeighbor(GetCBParticle(GetParticle(p_ID)->GetTypeID()),eps1, eps2, eps3);
        
        // checking distance w.r.t the particle sizes and tunneling allowance
        if (dist < ((dpi+dpj)/2.0) + (Td) || std::fabs(dist - ((dpi+dpj)/2.0) - (Td)) < 0.000001)
        {
            //Checking interference allowance of most restrictive particle
            if (dist > (((dpi+dpj)/2.0)-(mult*PS)) || std::fabs(dist-(((dpi+dpj)/2.0)+(mult*PS))) < 0.000001)
            {
                GetCBParticle(GetParticle(i)->GetTypeID())->AddNeighbor(p_ID, dist);
            }
            else
            {
                std::cout << "TOO CLOSE WHAT?" << std::endl;
            }
        }
    }
}

void RVEdata::Check_Adjacent_Indices(long int i)
{
    //Compare the distance of adjacent indices. Loop through the set of potential neighbor spots
    long int x = GetCBParticle(GetParticle(i)->GetTypeID())->GetXCoordinate();
    long int y = GetCBParticle(GetParticle(i)->GetTypeID())->GetYCoordinate();
    long int z = GetCBParticle(GetParticle(i)->GetTypeID())->GetZCoordinate();

    long int CB_index = ((Ybound+1)*(Zbound+1)*x) + ((Ybound+1)*z) + y;

    int iPSoftness, jPSoftness, dpi, dpj;
    double dist;
    
    int which_set = Get_Max_CB_Map_Index(i);
    long int which_index;
    long int max_index = Potential.size();

    for (int j = 0; j < CB_Potential_Maps_Indices[which_set].size(); j++)
    {
        //if (j == i) { continue;}
        long int p_ID;
        which_index = CB_index + CB_Potential_Maps_Indices[which_set][j];

        //Check the bounds of the potential map index
        if (x + CB_Potential_Maps_X[which_set][j] > Xbound || x + CB_Potential_Maps_X[which_set][j] < 0) {continue;}
        if (y + CB_Potential_Maps_Y[which_set][j] > Ybound || y + CB_Potential_Maps_Y[which_set][j] < 0) {continue;}
        if (z + CB_Potential_Maps_Z[which_set][j] > Zbound || z + CB_Potential_Maps_Z[which_set][j] < 0) {continue;}
        if (which_index > max_index) {continue;}
        //If passed the "if" statements, the index is within the RVE. Check it for a particle
        
        if(isNull(All_Indices[x+CB_Potential_Maps_X[which_set][j]], y+CB_Potential_Maps_Y[which_set][j], z+CB_Potential_Maps_Z[which_set][j]) == true)
        {
            //If == 0, then there is no particle at the index. Onto the next one!
            continue;
        }
        else
        {
            p_ID = All_Indices[x+CB_Potential_Maps_X[which_set][j]].coeffRef(y+CB_Potential_Maps_Y[which_set][j], z+CB_Potential_Maps_Z[which_set][j]);
        }
        
        int Particle_i_Type = GetParticle(i)->GetType();
        int Particle_j_Type = GetParticle(p_ID)->GetType();
        
        //Set particle softness of particle j
        if (Particle_i_Type <= NumCBFillers)
        {
            iPSoftness = GetCBMatType(GetCBParticle(GetParticle(i)->GetTypeID())->GetMatType())->GetPenetrationAllowance()*mult;
            dpi = GetCBParticle(GetParticle(i)->GetTypeID())->GetDiameter()*mult;
        }
        
        if (Particle_j_Type <= NumCBFillers)
        {
            jPSoftness = GetCBMatType(GetCBParticle(GetParticle(p_ID)->GetTypeID())->GetMatType())->GetPenetrationAllowance()*mult;
            dpj = GetCBParticle(GetParticle(p_ID)->GetTypeID())->GetDiameter()*mult;
        }
        
        //Define the most restrictive penetration allowance
        double PS;
        if (iPSoftness < jPSoftness)
        {
            PS = iPSoftness;
        }
        else
        {
            PS = jPSoftness;
        }
        
        dist = GetCBParticle(GetParticle(i)->GetTypeID())->FindDistanceToNeighbor(GetCBParticle(GetParticle(p_ID)->GetTypeID()));

        // checking distance w.r.t the particle sizes and tunneling allowance
        if (dist < ((dpi+dpj)/2.0) + (Td) || std::fabs(dist - ((dpi+dpj)/2.0) - (Td)) < 0.000001)
        {
            //Checking interference allowance of most restrictive particle
            if (dist > (((dpi+dpj)/2.0)-(mult*PS)) || std::fabs(dist-(((dpi+dpj)/2.0)+(mult*PS))) < 0.000001)
            {
                GetCBParticle(GetParticle(i)->GetTypeID())->AddNeighbor(p_ID, dist);
            }
        }
    }
}

bool RVEdata::Seed_Check_Adjacent_Indices(long int x, long int y, long int z, long int i)
{
    //Compare the distance of adjacent indices. Loop through the set of potential neighbor spots
    long int CB_index = ((Ybound+1)*(Zbound+1)*x) + ((Ybound+1)*z) + y;
    
    int iPSoftness, jPSoftness, dpi, dpj;
    double dist;
    
    int which_set = Get_Max_CB_Map_Index(i);
    long int which_index;
    long int max_index = Potential.size();
    
    for (int j = 0; j < CB_Potential_Maps_Indices[which_set].size(); j++)
    {
        //if (j == i) { continue;}
        long int p_ID;
        which_index = CB_index + CB_Potential_Maps_Indices[which_set][j];
        
        //Check the bounds of the potential map index
        if (x + CB_Potential_Maps_X[which_set][j] > Xbound || x + CB_Potential_Maps_X[which_set][j] < 0) {continue;}
        if (y + CB_Potential_Maps_Y[which_set][j] > Ybound || y + CB_Potential_Maps_Y[which_set][j] < 0) {continue;}
        if (z + CB_Potential_Maps_Z[which_set][j] > Zbound || z + CB_Potential_Maps_Z[which_set][j] < 0) {continue;}
        if (which_index > max_index) {continue;}
        //If passed the "if" statements, the index is within the RVE. Check it for a particle
        
        if(isNull(All_Indices[x+CB_Potential_Maps_X[which_set][j]], y+CB_Potential_Maps_Y[which_set][j], z+CB_Potential_Maps_Z[which_set][j]) == true)
        {
            //If == 0, then there is no particle at the index. Onto the next one!
            continue;
        }
        else
        {
            p_ID = All_Indices[x+CB_Potential_Maps_X[which_set][j]].coeffRef(y+CB_Potential_Maps_Y[which_set][j], z+CB_Potential_Maps_Z[which_set][j]);
        }
        
        int Particle_i_Type = GetParticle(i)->GetType();
        int Particle_j_Type = GetParticle(p_ID)->GetType();
        
        //Set particle softness of particle j
        if (Particle_i_Type <= NumCBFillers)
        {
            iPSoftness = GetCBMatType(GetCBParticle(GetParticle(i)->GetTypeID())->GetMatType())->GetPenetrationAllowance()*mult;
            dpi = GetCBParticle(GetParticle(i)->GetTypeID())->GetDiameter()*mult;
        }
        
        if (Particle_j_Type <= NumCBFillers)
        {
            jPSoftness = GetCBMatType(GetCBParticle(GetParticle(p_ID)->GetTypeID())->GetMatType())->GetPenetrationAllowance()*mult;
            dpj = GetCBParticle(GetParticle(p_ID)->GetTypeID())->GetDiameter()*mult;
        }
        
        //Define the most restrictive penetration allowance
        double PS;
        if (iPSoftness < jPSoftness)
        {
            PS = iPSoftness;
        }
        else
        {
            PS = jPSoftness;
        }
        
        dist = GetCBParticle(GetParticle(i)->GetTypeID())->FindDistanceToNeighbor(GetCBParticle(GetParticle(p_ID)->GetTypeID()));
        
        // checking distance w.r.t the particle sizes and tunneling allowance
        if (dist < ((dpi+dpj)/2.0) + (Td) || std::fabs(dist - ((dpi+dpj)/2.0) - (Td)) < 0.000001)
        {
            //Checking interference allowance of most restrictive particle
            if (dist > (((dpi+dpj)/2.0)-(mult*PS)) || std::fabs(dist-(((dpi+dpj)/2.0)+(mult*PS))) < 0.000001)
            {
                //Do nothing, it is valid
            }
            else
            {
                return false;
            }
        }
    }
    return true;
}

bool RVEdata::isNull(const Eigen::SparseMatrix<double, Eigen::ColMajor>& mat, int row, int col)
{
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(mat, col); it; ++it)
    {
        if (it.row() == row) {return false;}
    }
    return true;
}

Particles* RVEdata::GetParticle(long int i)
{
    return &particle_data[i];
}

CBparticle* RVEdata::GetCBParticle(long int num)
{
    return &CB_particles[num];
}

CNTparticle* RVEdata::GetCNTParticle(long int num)
{
    return &CNT_particles[num];
}

CBMaterialData* RVEdata::GetCBMatType(int i)
{
    return &CB_data[i];
}

PolymerData* RVEdata::GetPolymerData()
{
    return &p_data;
}

long int RVEdata::GetXBoundaryLimit()
{
    return Xbound;
}

bool RVEdata::SolveForStrainConductivity(float strain_elong, std::ofstream &cond_out)
{
    float eps1 = 1 + strain_elong;
    float Poissons = p_data.GetPoissonsRatio();
    float eps2, eps3;
    if (Poissons != 0)
    {
        eps2 = 1 - (strain_elong*Poissons);
        eps3 = 1 - (strain_elong*Poissons);
    }
    else
    {
        eps2 = 1;
        eps3 = 1;
    }
    
    //Vectors passed to the find network algorithm
    long int N = particle_data.size();
    
    std::deque<std::mutex> PosLocks(N);
    std::deque<std::mutex> GndLocks(N);
    
    posbank.resize(N, std::vector<bool> (2));
    gndbank.resize(N, std::vector<bool> (2));
    network.resize(N);
    
    //Initialize vector lock objects
    for (long int i = 0; i<N; i++)
    {
        posbank[i][0] = false;
        gndbank[i][0] = false;
        posbank[i][1] = false;
        gndbank[i][1] = false;
    }
    
    //Find the networked particles
    NetExists = FindStrainNetwork(N, PosLocks, GndLocks, eps1);
    
    if (NetExists == false)
    {
        //std::cout << "No Net" << std::endl;
        RVEconductivity = 0;
        RVEresistance = 0;
        cond_out << strain_elong << "," << RVEresistance << "," << RVEconductivity << "," << std::endl;
        return false;
    }
    //std::cout << "Net Exists" << std::endl;
    
    //If arrived at this point, solve for conducitiviy of network particles
    SolveForNodalVoltages();
    
    RVEresistance = VIN/sumCurrent;
    sumCurrent = 0;
    
    RVEconductivity = (Xbound*d*eps1)/(Ybound*Zbound*d*d*eps2*eps3*RVEresistance);
    //std::cout << "Resistance = " << RVEresistance << ", strain=" << strain_elong << std::endl;
    cond_out << strain_elong << "," << RVEresistance << "," << RVEconductivity << "," << std::endl;
    return true;
}

bool RVEdata::SolveForConductivity()
{
    //Vectors passed to the find network algorithm
    long int N = particle_data.size();
    
    std::deque<std::mutex> PosLocks(N);
    std::deque<std::mutex> GndLocks(N);
    
    posbank.resize(N, std::vector<bool> (2));
    gndbank.resize(N, std::vector<bool> (2));
    network.resize(N);
    
    //Initialize vector lock objects
    for (long int i = 0; i<N; i++)
    {
        posbank[i][0] = false;
        gndbank[i][0] = false;
        posbank[i][1] = false;
        gndbank[i][1] = false;
    }
    
    //Find the networked particles
    NetExists = FindNetwork(N, PosLocks, GndLocks);
    
    if (NetExists == false)
    {
        //std::cout << "No Net" << std::endl;
        RVEconductivity = 0;
        RVEresistance = 0;
        return false;
    }
    //std::cout << "Net Exists" << std::endl;
    
    //If arrived at this point, solve for conducitiviy of network particles
    SolveForNodalVoltages();

    RVEresistance = VIN/sumCurrent;
    RVEconductivity = (Xbound*d)/(Ybound*Zbound*d*d*RVEresistance);
    sumCurrent = 0;
    //std::cout << "Resistance = " << RVEresistance << std::endl;
    return true;
}

bool RVEdata::FindNetwork(long int N, std::deque<std::mutex> &PosLocks, std::deque<std::mutex> &GndLocks)
{
    //Fill initial posbank and gndbank first rows
    for (long int i = 0; i<particle_data.size(); i++)
    {
        if (particle_data[i].GetType() < NumCBFillers)
        {
            //It's a cb particle
            posbank[i][0] = CB_particles[particle_data[i].GetTypeID()].IsPositive();
            gndbank[i][0] = CB_particles[particle_data[i].GetTypeID()].IsNegative();
        }
        else
        {
            //It's a cnt particle
        }
    }
    
    //Count num of processesors available
    int NumProcessors = std::thread::hardware_concurrency();

    //Half of available processors do positive search, half do negative search
    std::vector<std::future<void> > BranchSearchSplit;
    
    BranchSearchSplit.push_back(std::async(std::launch::async, &RVEdata::PositiveBranch, this, NumProcessors/2, std::ref(PosLocks)));
    BranchSearchSplit.push_back(std::async(std::launch::async, &RVEdata::NegativeBranch, this, NumProcessors/2, std::ref(GndLocks)));
    BranchSearchSplit[0].wait();
    BranchSearchSplit[1].wait();
    
    //Check for Network
    int NumNetwork = 0;
    
    for (unsigned long int i = 0; i<N; i++)
    {
        if (posbank[i][1] == true && gndbank[i][1] == true)
        {
            network[i] = 1;
            NumNetwork += 1;
        }
        else
        {
            network[i] = 0;
        }
    }
    
    if (NumNetwork > 0)
    {
        return(true);
    }
    else
    {
        return(false);
    }
}

bool RVEdata::FindStrainNetwork(long int N, std::deque<std::mutex> &PosLocks, std::deque<std::mutex> &GndLocks, float eps1)
{
    //Fill initial posbank and gndbank first rows
    for (long int i = 0; i<particle_data.size(); i++)
    {
        if (particle_data[i].GetType() < NumCBFillers)
        {
            //It's a cb particle
            if (CB_particles[particle_data[i].GetTypeID()].Strain_IsPositive(eps1,Xbound,Td) == true)
            {
                posbank[i][0] = true;
            }
            else
            {
                posbank[i][0] = false;
                particle_data[i].SetPosTerminalResistance(0);
            }
            if (CB_particles[particle_data[i].GetTypeID()].Strain_IsNegative(eps1,Xbound,Td) == true)
            {
                gndbank[i][0] = true;
            }
            else
            {
                gndbank[i][0] = false;
            }
        }
        else
        {
            //It's a cnt particle
        }
    }
    
    //Count num of processesors available
    int NumProcessors = std::thread::hardware_concurrency();
    
    //Half of available processors do positive search, half do negative search
    std::vector<std::future<void> > BranchSearchSplit;
    
    BranchSearchSplit.push_back(std::async(std::launch::async, &RVEdata::PositiveBranch, this, NumProcessors/2, std::ref(PosLocks)));
    BranchSearchSplit.push_back(std::async(std::launch::async, &RVEdata::NegativeBranch, this, NumProcessors/2, std::ref(GndLocks)));
    BranchSearchSplit[0].wait();
    BranchSearchSplit[1].wait();
    
    //Check for Network
    int NumNetwork = 0;
    
    for (unsigned long int i = 0; i<N; i++)
    {
        if (posbank[i][1] == true && gndbank[i][1] == true)
        {
            network[i] = 1;
            NumNetwork += 1;
        }
        else
        {
            network[i] = 0;
        }
    }
    
    if (NumNetwork > 0)
    {
        return(true);
    }
    else
    {
        return(false);
    }
}

void RVEdata::PositiveBranch(int NumProcessors, std::deque<std::mutex> &PosLocks)
{
    //initializing variables that help identify when all pos network has been found
    int NumPos2 = 0;
    int NumPos = 1;
    
    // This while loop is dependent on the two positive network identifiers. In each iteration
    // the code will find positive network particles identified by a 1 in the first row of the
    // posbank matrix. If the particle has a "0" in the second row and a 1 in the first row,
    // this signals that the particle has not yet been processed. The processing involves searching
    // the particle for its neighbors, and changing those neighbor indexes in the 1st row of the
    // posbank to 1's. This process is repeated until there are no more changes in the number of
    // 1's in the posbank from the beginning of the loop to the end of the loop.
    
    while (NumPos != NumPos2)
    {
        // Finding the total number of particles in the positive network at the loop start
        NumPos = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (posbank[w][1] == true)
            {
                NumPos += 1;
            }
        }
        
        //Split the search party into futures
        std::vector<std::future<void> > SearchersPos;
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersPos.push_back(std::async(std::launch::async, &RVEdata::PositiveSearch, this, std::ref(PosLocks)));
        }
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersPos[i].wait();
        }
        
        // Finding the total number of particles in the pos network at the end of the loop
        NumPos2 = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (posbank[w][1] == true)
            {
                NumPos2 += 1;
            }
        }
    }
    //std::cout << "NumPos2 = " << NumPos2 << std::endl;
}

//Algorithm that searches the posbank for particles that need to be searched and have their neighbors added to posbank
void RVEdata::PositiveSearch(std::deque<std::mutex> &PosLocks)
{
    //Neighbors vector placeholder
    std::vector<long int> neighbors;
    unsigned long int NumNeigh = 0;
    
    for (unsigned long int i = 0; i<particle_data.size(); i++)
    {
        //If try_lock is successful, this thread has ownership of the index i search
        if (PosLocks[i].try_lock() == true)
        {
            //If statement conditions below holds true, the particle neighbors need to be added to the bank
            if (posbank[i][0] == true && posbank[i][1] == false)
            {
                neighbors = GetNeighbors(i);
                NumNeigh = neighbors.size();
                if (NumNeigh != 0)
                {
                    //Loop through neighbors and them to the bank
                    for (unsigned int b = 0; b<NumNeigh; b++)
                    {
                        bool Success = false;
                        int FailCount = 0;
                        while (Success == false)
                        {
                            //If try_lock succeeds, add neighbor to bank. If fails, loop again and wait for successful lock (max of 5 times)
                            if (PosLocks[neighbors[b]].try_lock() == true)
                            {
                                posbank[neighbors[b]][0] = true;
                                PosLocks[neighbors[b]].unlock();
                                Success = true;
                            }
                            else
                            {
                                FailCount += 1;
                            }
                            
                            //If try_lock fails 5 times, move on to the next while marking this index as needing further attention
                            if (FailCount > 5)
                            {
                                PosLocks[i].unlock();
                                goto TryLaterNeg;
                            }
                        }
                    }
                }
                posbank[i][1] = 1;
            }
            PosLocks[i].unlock();
        TryLaterNeg:
            continue;
        }
        else
        {
            continue;
        }
    }
}

void RVEdata::NegativeBranch(int NumProcessors, std::deque<std::mutex> &GndLocks)
{
    // Finding the total number of particles direclty connected to positive electrode
    int NumGnd = 1;
    int NumGnd2 = 0;
    // This while loop is dependent on the two negative network identifiers. In each iteration
    // the code will find negative network particles identified by a 1 in the first row of the
    // gndbank matrix. If the particle has a "0" in the second row and a 1 in the first row,
    // this signals that the particle has not yet been processed. The processing involves searching
    // the particle for its neighbors, and changing those neighbor indexes in the 1st row of the
    // gndbank to 1's. This process is repeated until there are no more changes in the number of
    // 1's in the gndbank from the beginning of the loop to the end of the loop.

    while (NumGnd != NumGnd2)
    {
        // Finding the total number of particles in the ground network at the loop start
        NumGnd = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (gndbank[w][1] == true)
            {
                NumGnd += 1;
            }
        }
        
        //Split the search party into futures
        std::vector<std::future<void> > SearchersNeg;
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersNeg.push_back(std::async(std::launch::async, &RVEdata::NegativeSearch, this, std::ref(GndLocks)));
        }
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersNeg[i].wait();
        }
        
        // Finding the total number of particles in the ground network at the loop end
        NumGnd2 = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (gndbank[w][1] == true)
            {
                NumGnd2 += 1;
            }
        }
    }
    //std::cout << "NumGnd2 = " << NumGnd2 << std::endl;
}

//Algorithm that searches the gndbank for particles that need to be searched and have their neighbors added to gndbank
void RVEdata::NegativeSearch(std::deque<std::mutex> &GndLocks)
{
    //Neighbors vector placeholder
    std::vector<long int> neighbors;
    unsigned long int NumNeigh = 0;

    for (unsigned long int i = 0; i<particle_data.size(); i++)
    {
        //If try_lock is successful, this thread has ownership of the index i search
        if (GndLocks[i].try_lock() == true)
        {
            //If statement conditions below holds true, the particle neighbors need to be added to the bank
            if (gndbank[i][0] == true && gndbank[i][1] == false)
            {
                neighbors = GetNeighbors(i);
                NumNeigh = neighbors.size();
                if (NumNeigh != 0)
                {
                    //Loop through neighbors and them to the bank
                    for (unsigned int b = 0; b<NumNeigh; b++)
                    {
                        bool Success = false;
                        int FailCount = 0;
                        while (Success == false)
                        {
                            //If try_lock succeeds, add neighbor to bank. If fails, loop again and wait for successful lock (max of 5 times)
                            if (GndLocks[neighbors[b]].try_lock() == true)
                            {
                                gndbank[neighbors[b]][0] = true;
                                GndLocks[neighbors[b]].unlock();
                                Success = true;
                            }
                            else
                            {
                                FailCount += 1;
                            }
                            
                            //If try_lock fails 5 times, move on to the next while marking this index as needing further attention (i.e. don't mark true in the [i][1] position

                            if (FailCount > 5)
                            {
                                GndLocks[i].unlock();
                                goto TryLaterNeg;
                            }
                        }
                    }
                }
                gndbank[i][1] = 1;
            }
            GndLocks[i].unlock();
        TryLaterNeg:
            continue;
        }
        else
        {
            continue;
        }
    }
}

void RVEdata::OutputResults(std::ofstream &fout)
{
    fout << RVEresistance << "," << RVEconductivity << "," << std::endl;
}

//return vector of neighors to particle i
std::vector<long int> &RVEdata::GetNeighbors(long int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        return CB_particles[particle_data[i].GetTypeID()].GetNeighborsVector();
    }
    else
    {
        //Do nothing
        std::cout << "ERROR: CNT particles don't exist yet" << std::endl;
        exit(1);
    }
}

std::vector<double> &RVEdata::GetNeighborsDist(long int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        return CB_particles[particle_data[i].GetTypeID()].GetNeighborsDistVector();
    }
    else
    {
        //Do nothing
        std::cout << "ERROR: CNT particles don't exist yet" << std::endl;
        exit(1);
    }
}

//Function used to solve for the nodal voltages of the particles in the rve.
void RVEdata::SolveForNodalVoltages()
{
    //Reset all the storage
    //Shell.reset();
    //Potential.reset();
    //All_Indices.resize(0);
    
    int VIN = 5;
    
    //Count num particles in network and assign row indices to network particles
    std::vector<long int> Indices;
    std::vector<long int> NeighborIndices;
    long int SumNet = 0;
    for (long int i = 0; i < network.size(); i++)
    {
        if (network[i] == true)
        {
            Indices.push_back(i);
            NeighborIndices.push_back(SumNet);
            SumNet += 1;
        }
        else
        {
            NeighborIndices.push_back(-1);
        }
    }
    
    //double Rc = 10;          //Contact resistance between carbon allotropes in ohms
    
    Eigen::SparseMatrix<double, Eigen::ColMajor> g (SumNet, SumNet);
    g.reserve(10);
    Eigen::VectorXd curr = Eigen::VectorXd::Zero(SumNet);
    
    //Count num of processesors available
    int NumProcessors = std::thread::hardware_concurrency();
    
    //Find the particle distribution between cores
    long int NumP = SumNet / NumProcessors;
    long int Left = SumNet % NumProcessors;
    
    //Place particles in the SubRVEs using threads
    std::vector<std::future<std::vector<Eigen::Triplet<double> > > > Waiters;
    typedef Eigen::Triplet<double> T;
    std::vector<T> GtripletList;
    
    std::vector<std::vector<long int> > IndicesG;
    std::vector<std::vector<double> > ValuesG;
    std::vector<double> ValuesCurr;
    IndicesG.resize(SumNet);
    ValuesG.resize(SumNet);
    ValuesCurr.resize(SumNet);
    
    long int start = 0;
    long int finish;
    
    auto FillStart = std::chrono::high_resolution_clock::now();
    //std::cout << "Filling Cond, NumProc = " << NumProcessors << std::endl;
    
    //Loop to fill G and Curr using multiple threads
    for (int i = 0; i<NumProcessors; i++)
    {
        if (i == NumProcessors - 1)
        {
            finish = SumNet;
        }
        else
        {
            finish = start + NumP;
        }
        Waiters.push_back(std::async(std::launch::async, &RVEdata::FillCurrAndG, this, std::ref(ValuesCurr), std::ref(Indices), std::ref(NeighborIndices), Rc, VIN, start, finish));
        start += NumP;
        //std::cout << "start = " << start << "\r" << std::flush;
    }
    
    std::vector<Eigen::Triplet<double> > Gtemp;
    
    //Fill the triplet list created by the threads as they are returned
   for (int i = 0; i<Waiters.size(); i++)
    {
        Waiters[i].wait();
        Gtemp = Waiters[i].get();
        GtripletList.insert(GtripletList.end(), Gtemp.begin(), Gtemp.end());
        Gtemp.clear();
    }
    
    //Fill G and Curr Eigen objects using triplet list
    for (long int i = 0; i<SumNet; i++)
    {
        //Current vector portion
        if (ValuesCurr[i] > 0)
        {
            curr(i) = ValuesCurr[i];
        }
    }
    
    //Fill the cond matrix using the triplets created by the thread
    g.setFromTriplets(GtripletList.begin(), GtripletList.end());
    //std::cout << "G set = " << std::endl;
    
    auto FillFinish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> FilledElapsed = FillFinish - FillStart;
    //std::cout << FilledElapsed.count() << "..Solving Linear Alg " << std::flush;

    //Perform linear algebra to obtain nodal voltages and output info for debugging
    g.finalize();
    
    /*
    //DEBUG G AND I
    std::ofstream checkg;
    checkg.open("checkg.csv");
    for (int q = 0; q < g.rows(); q++)
    {
        for (int r = 0; r < g.cols(); r++)
        {
            checkg << g.coeffRef(q,r) << ",";
        }
        checkg << std::endl;
    }
    checkg << std::endl;
    
    for (int q = 0; q < curr.rows(); q++)
    {
        checkg << curr.coeffRef(q) << std::endl;
    }
    checkg.close();
    */
    
    //Eigen::SimplicialLDLT <Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > solver;
    Eigen::SimplicialLDLT <Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;

    g.makeCompressed();
    //std::cout << "G Compressed " << std::endl;
    solver.analyzePattern(g);
    //std::cout << "G Analyzed " << std::endl;
    solver.compute(g);
    //std::cout << "G Computed " << std::endl;
    Eigen::VectorXd V = Eigen::VectorXd::Zero(SumNet,1);
    solver.factorize(g);
    //std::cout << "G Factorized " << std::endl;

    if (solver.info() != 0)
    {
        std::cout << std::endl << "Matrix Factorization Failed!" << std::endl;//<< solver.lastErrorMessage() << std::endl;
        exit (1);
    }
    V = solver.solve(curr);
    
    //std::cout << "V Solved = " << std::endl;

    //CheckNeighbors();
    auto Solved = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> SolvedFinal = Solved - FillFinish;
    //std::cout << SolvedFinal.count() << ".." << std::flush;
    
    //Fill Voltage vector in RVEdata class with nodal voltages found
    Nodal_Voltages.resize(particle_data.size());
    int ph = 0;
    for (long int i = 0; i < particle_data.size(); i++)
    {
        if (network[i] == false)
        {
            Nodal_Voltages[i] = 0;
            continue;
        }
        Nodal_Voltages[i] = V.coeffRef(ph);
        ph += 1;
        
        if (Nodal_Voltages[i]-VIN > 0.001)
        {
            std::cout << "i: " << i << " V: " << Nodal_Voltages[i] << std::endl;
        }
    }
    
    //Sum the current flowing into the positively connected particles and calculate resistance/conductivity
    for (long int i = 0; i < SumNet; i++)
    {
        if (particle_data[Indices[i]].GetPosTerminalResistance() != 0)
        {
            sumCurrent += (VIN - Nodal_Voltages[Indices[i]])/particle_data[Indices[i]].GetPosTerminalResistance();
        }
    }
}

//Function that is threaded to fill temporary current and conductivity entries for the matrices required to solve for nodal voltages via KVL
std::vector<Eigen::Triplet<double> > RVEdata::FillCurrAndG(std::vector<double> &ValuesCurr, std::vector<long int> &Indices, std::vector<long int> &NeighborIndices, double Rc, int VIN, long int start, long int finish)
{
    //Neighbor vector placeholder and distance placeholder and tunneling threshold
    std::vector<long int> neighbors;
    std::vector<double> neighborsDist;
    std::vector<Eigen::Triplet<double> > G;
    double TunnelingThresh;
    double TunnelingResistance;
    double ValuesG;
    
    for (long int i = start; i<finish; i++)
    {        
        //Fill Conductivity matrix using data from neighbor vector
        neighbors = GetNeighbors(Indices[i]);
        neighborsDist = GetNeighborsDist(Indices[i]);
        for (int j = 0; j < neighbors.size(); j++)
        {
            TunnelingThresh = CalcTunnelingThreshold(Indices[i],neighbors[j]);
            if (neighborsDist[j] < TunnelingThresh || std::abs(neighborsDist[j] - TunnelingThresh) < 0.0001)
            {
                //Particles are in CONTACT
                ValuesG = -1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,NeighborIndices[neighbors[j]],ValuesG));
                ValuesG = 1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
            else
            {
                long double At = FindTunnelingArea(Indices[i],neighbors[j]);
                TunnelingResistance = CalcTunnelingResistance(neighborsDist[j], TunnelingThresh, d, At);
                if (TunnelingResistance < Rc)
                {
                    TunnelingResistance = Rc;
                }
                //std::cout << "TR = " << TunnelingResistance << std::endl;
                ValuesG = -1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,NeighborIndices[neighbors[j]],ValuesG));
                ValuesG = 1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
        }
        //Check particle i's relation to positive and negative terminals
        if (CheckPosNeg(Indices[i],false) == true)
        {
            //Find distance to positive terminal and place conductivity in g matrix and fill current vector with corresponding value
            double dist = FindDistanceToTerminal(false,Indices[i]);
            TunnelingThresh = CalcTunnelingThresholdToTerminal(Indices[i]);
            if (dist < TunnelingThresh || std::abs(dist - TunnelingThresh) < 0.0001)
            {
                //Particles are in CONTACT
                particle_data[Indices[i]].SetPosTerminalResistance(Rc);
                ValuesCurr[i] = (VIN/Rc);
                ValuesG = 1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
            else
            {
                //Particles are tunneling
                long double At = FindTunnelingArea(Indices[i],-1);
                TunnelingResistance = CalcTunnelingResistance(dist, TunnelingThresh, d, At);
                if (TunnelingResistance < Rc)
                {
                    TunnelingResistance = Rc;
                }
                
                particle_data[Indices[i]].SetPosTerminalResistance(TunnelingResistance);
                ValuesCurr[i] = (VIN/TunnelingResistance);
                ValuesG = 1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
        }
        if (CheckPosNeg(Indices[i],true) == true)
        {
            //Find distance to negative terminal and place conductivity in g matrix
            double dist = FindDistanceToTerminal(true,Indices[i]);
            TunnelingThresh = CalcTunnelingThresholdToTerminal(Indices[i]);
            if (dist < TunnelingThresh || std::abs(dist - TunnelingThresh) < 0.0001)
            {
                //Particles are in CONTACT
                ValuesG = 1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
            else
            {
                //particles are tunneling
                long double At = FindTunnelingArea(Indices[i],-1);
                TunnelingResistance = CalcTunnelingResistance(dist, TunnelingThresh, d, At);
                if (TunnelingResistance < Rc)
                {
                    TunnelingResistance = Rc;
                }
                ValuesG = 1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
        }
    }
    return G;
}

double RVEdata::CalcTunnelingThreshold(long int i, long int j)
{
    int iType, jType;
    iType = particle_data[i].GetType();
    jType = particle_data[j].GetType();

    if (iType < NumCBFillers && jType < NumCBFillers)
    {
        //Both are CB particles
        int dpi, dpj;
        dpi = CB_particles[particle_data[i].GetTypeID()].GetDiameter()*mult;
        dpj = CB_particles[particle_data[j].GetTypeID()].GetDiameter()*mult;
        return ((dpi+dpj)/2.0);
    }
    else if (iType > NumCBFillers && jType > NumCBFillers)
    {
        //Both i and j are CNT particles
        return 0;
    }
    else if (iType > NumCBFillers && jType < NumCBFillers)
    {
        //i is a cnt and j is a cb
        return 0;
    }
    else if (iType < NumCBFillers && jType > NumCBFillers)
    {
        //j is a cnt and i is a cb
        return 0;
    }
    else
    {
        return 0;
    }
}

double RVEdata::CalcTunnelingThresholdToTerminal(long int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        //i is a cb particle
        return mult*(CB_particles[particle_data[i].GetTypeID()].GetDiameter()/2.0);
    }
    else
    {
        return 0;
    }
}

double RVEdata::CalcTunnelingResistance(double dist, double TunnelingThresh, double u, long double At)
{
    //Tunneling Conductivity Parameter definition
    long double h = 6.626*pow(10,-34);           //Planck's constant (m^2kg/s)
    long double e = 1.602176634*pow(10,-19);     //Elementary charge in units of Coulombs
    long double tau = 0.5*1.60218*pow(10,-19);   //Average barrier height of Epoxy
    long double Em = 9.10938356*pow(10,-31);     //electron mass
    double d;                                   //tunneling distance
    
    d = (dist - TunnelingThresh)*u/mult;
    
    if (d < 0)
    {
        std::cout << "Can't have negative tunneling distance.." << std::endl;
        exit(1);
    }
    
    long double RtunnNumerator = exp(4.0*M_PI*d*sqrt(2.0 * Em*tau) / h)*(h*h*d);
    long double RtunnDenom = (At*e*e*sqrt(2.0 * Em*tau));
    
    double TR = RtunnNumerator/RtunnDenom;
    
    if (dist - TunnelingThresh > Td)
    {
        std::cout << "Why? How? Dist = " << dist << " TT = " << TunnelingThresh << std::endl;
    }
    
    if (TR == 0)
    {
        std::cout << "Tunn R can't be 0: Rnumerator = " << RtunnNumerator << " dist= " << dist << " thresh= " << TunnelingThresh << " d= " << d << " u = " << u << " 2.0*EM*tau= " << 2.0*Em*tau << std::endl;
    }
    else if (TR < 0)
    {
        std::cout << "Tunn R can't be less than 0: Rnumer = " << RtunnNumerator << " Rdenom= " << RtunnDenom << std::endl;
    }
    if (At == 0)
    {
        std::cout << " At: " << At;
    }
    return TR;
}

double  RVEdata::FindDistanceToTerminal(bool terminal, long int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        //i is a cb particle
        if (terminal == false)
        {
            //positive terminal
            return mult*CB_particles[particle_data[i].GetTypeID()].GetXCoordinate();
        }
        else
        {
            //ground terminal
            return mult*(Xbound - CB_particles[particle_data[i].GetTypeID()].GetXCoordinate());
        }
    }
    else
    {
        return 0;
    }
}

long double RVEdata::FindTunnelingArea(long int i, long int j)
{
    if (j < 0)
    {
        //Particle to terminal tunneling
        if (particle_data[i].GetType() <= NumCBFillers)
        {
            //CB particle to terminal
            double dp = mult*CB_particles[particle_data[i].GetTypeID()].GetDiameter()*pow(10,-9);
            return M_PI*pow(dp/2.0,2);           //cross sectional area of tunneling path (m^2)
        }
        else
        {
            //CNT particle to terminal
            return 0;
        }
    }
    if (j >= 0)
    {
        //Particle to particle tunneling
        if (particle_data[i].GetType() <= NumCBFillers && particle_data[j].GetType() <= NumCBFillers)
        {
            //CB to CB tunneling
            double dpi = mult*CB_particles[particle_data[i].GetTypeID()].GetDiameter()*pow(10,-9);
            double dpj = mult*CB_particles[particle_data[j].GetTypeID()].GetDiameter()*pow(10,-9);
            return M_PI*pow((dpi+dpj)/4.0,2);           //cross sectional area of tunneling path (m^2)
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

bool RVEdata::CheckPosNeg(long int i, bool j)
{
    int Type = particle_data[i].GetType();
    
    if (j == true)
    {
        //Check to Negative terminal
        if (Type < NumCBFillers)
        {
            //CB particle
            return CB_particles[particle_data[i].GetTypeID()].IsNegative();
        }
        else
        {
            //Not CB particle
            return 0;
        }
    }
    if (j == false)
    {
        //Check to positive terminal
        if (Type < NumCBFillers)
        {
            //CB particle
            return CB_particles[particle_data[i].GetTypeID()].IsPositive();
        }
        else
        {
            //Not CB particle
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

void RVEdata::CheckNeighbors()
{
    std::string dir = getcwd(NULL, 0);
    std::ofstream MFErr;
    MFErr.open(dir + "/CheckNeigh.csv");
    long int NumNeighbors = 0;
    for (long int i = 0; i < particle_data.size(); i++)
    {
        std::vector<long int> neigh = CB_particles[i].GetNeighborsVector();
        MFErr << "i: " << i << "," << CB_particles[i].GetXCoordinate() << "," << CB_particles[i].GetYCoordinate() << "," << CB_particles[i].GetZCoordinate() << ",";
        for (int j = 0; j < neigh.size(); j++)
        {
            MFErr << neigh[j] << ",";
            NumNeighbors += 1;
        }
        MFErr << std::endl;
    }
    MFErr.close();
    std::cout << "NumNeighors = " << NumNeighbors << std::endl;
    
    std::ofstream XYZ;
    XYZ.open(dir + "/check_xyz.csv");
    XYZ << Total_Num_Particles << "," << Xbound*mult << "," << Ybound*mult << "," << Zbound*mult << "," << 5 << "," << 1 << "," << 1000 << std::endl;
    
    for (int i = 0; i < Total_Num_Particles; i++)
    {
        XYZ << CB_particles[i].GetXCoordinate()*mult << "," << CB_particles[i].GetYCoordinate()*mult << "," << CB_particles[i].GetZCoordinate()*mult << std::endl;
    }
    XYZ.close();
}

void RVEdata::Output_Ovito_File(std::ofstream &fplot, std::string uid, std::string bashid, std::ofstream &cond_out)
{
    std::string current_dir = getcwd(NULL,0);
    std::string directory("/RVE_Data/" + std::to_string(Xbound*mult) + "_" + std::to_string(Wcb));
    std::string cond_directory("/Output/" + std::to_string(Xbound*mult) + "_" + std::to_string(Wcb));
    std::string New_Directory(current_dir + directory);
    std::string New_Cond_Directory(current_dir + cond_directory);

    //std::cout << "New directory name = " << New_Directory << std::endl;
    if (boost::filesystem::exists(New_Directory) == false)
    {
        boost::filesystem::create_directory(New_Directory);
        //std::cout << "New directory name = " << New_Directory << std::endl;
    }
    if (boost::filesystem::exists(New_Cond_Directory) == false)
    {
        boost::filesystem::create_directory(New_Cond_Directory);
        //std::cout << "New directory name = " << New_Directory << std::endl;
    }
    
    bool opened_file = false;
    int file_num = 0;
    while (opened_file == false)
    {
        if (boost::filesystem::exists(New_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(Pct_Fill) + "_" + std::to_string(file_num) + ".xyz") == true)
        {
            //Increment file number and try to save again
            file_num += 1;
        }
        else
        {
            fplot.open(New_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(Pct_Fill) + "_" + std::to_string(file_num) + ".xyz");
            cond_out.open(New_Cond_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(Pct_Fill) + "_" + std::to_string(file_num) + ".csv");
            opened_file = true;
        }
        if (file_num > 100)
        {
            std::cout << "File Inc: " << file_num << "\r" << std::flush;
        }
    }
    cond_out.close();
    fplot.close();
    fplot.open(New_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(Pct_Fill) + "_" + std::to_string(file_num) + ".xyz");
    fplot << Total_Num_Particles << std::endl << "X Y Z Diameter ColorScale" << std::endl;
    for (long int i = 0; i < Total_Num_Particles; i++)
    {
        fplot << GetCBParticle(i)->GetXCoordinate()*mult << " " << GetCBParticle(i)->GetYCoordinate()*mult << " " << GetCBParticle(i)->GetZCoordinate()*mult << " " << GetCBParticle(i)->GetDiameter()*mult/2 << std::endl;
    }
    fplot.close();
    
    cond_out.open(New_Cond_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(Pct_Fill) + "_" + std::to_string(file_num) + ".csv");
    cond_out << 0 << "," << RVEresistance << "," << RVEconductivity << "," << std::endl;

    /*
    while (opened_file == false)
    {
        if (boost::filesystem::exists(New_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(file_num) + ".vtk") == true)
        {
            //Increment file number and try to save again
            file_num += 1;
        }
        else
        {
            fplot.open(New_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(file_num) + ".vtk");
            opened_file = true;
        }
        if (file_num > 100)
        {
            std::cout << "File Inc: " << file_num << "\r" << std::flush;
        }
    }
    fplot.close();
    std::cout << "While done" << std::endl;
    float* new_arr = new float[Total_Num_Particles*3]();
    int* vardim = new int[2]();
    vardim[0] = 1;
    vardim[1] = 1;
    float **vars;
    vars = new float*[2];
    for (int i = 0; i < 2; i++)
    {
        vars[i] = new float[Total_Num_Particles];
    }
    const char *varnames[] = {"id", "diameter"};
    
    std::vector<int> p_coords(Total_Num_Particles*3);
    for (long int i = 0; i < Total_Num_Particles; i++)
    {
        new_arr[i*3] = GetCBParticle(i)->GetXCoordinate()*mult;
        new_arr[i*3 + 1] = GetCBParticle(i)->GetYCoordinate()*mult;
        new_arr[i*3 + 2] = GetCBParticle(i)->GetZCoordinate()*mult;
        vars[0][i] = 1;
        vars[1][i] = GetCBParticle(i)->GetDiameter()*mult;
        //fplot << GetCBParticle(i)->GetXCoordinate() << " " << GetCBParticle(i)->GetYCoordinate() << " " << GetCBParticle(i)->GetZCoordinate() << " 1 " << GetCBParticle(i)->GetDiameter()/2. << std::endl;
    }
    std::string vtk_name = New_Directory + "/" + bashid + uid + "_" + std::to_string(Zbound*mult) + "_" + std::to_string(Pct_Rand*1000) + "_" + std::to_string(file_num) + ".vtk";
    write_point_mesh(vtk_name.c_str(), 1, Total_Num_Particles, new_arr, 2, vardim, varnames, vars);*/
}

