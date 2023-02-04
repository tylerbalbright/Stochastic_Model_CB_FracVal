#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <future>
#include <chrono>
#include "Files.h"
#include "RVEdata.h"
#include "CBMaterialData.h"
#include "CBparticle.h"
#include "PolymerData.h"

int cmd_PctRandLower;
int cmd_PctRandUpper;
int cmd_XY;
int cmd_thickness;
int cmd_wcb;
int cmd_fill;
double cmd_Rc;
std::string uid;
std::string bid;

int main(int argc, char *argv[])
{
    cmd_PctRandLower = atoi(argv[1]);
    cmd_PctRandUpper = atoi(argv[2]);
    cmd_XY = atoi(argv[3]);
    cmd_thickness = atoi(argv[4]);
    uid = argv[5];
    bid = argv[6];
    cmd_wcb = atoi(argv[7]);
    cmd_fill = atoi(argv[8]);
    cmd_Rc = atof(argv[9]);
    
    std::vector<float> strain_elong{0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004};
    
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    
    std::ifstream fin;      //stream for problem input
    std::ofstream fout;      //stream for problem output
    std::ofstream ferr;      //stream for error output
    std::ofstream fplot;     //stream for plotting RVE
    std::ofstream cond_out;  //stream for reporting conductivities before/after strain
    
    RVEdata dat;        //object to store problem data
    Files o_files;      //object to manipulate files
    
    //Open all file streams
    o_files.OpenFiles(fin, fout, fplot, ferr);
    
    //Read in data for problem
    dat.ReadInputData(fin, ferr, cmd_PctRandLower, cmd_PctRandUpper, cmd_XY, cmd_thickness, cmd_wcb, cmd_fill, cmd_Rc);
    fin.close();
    //std::cout << "Reading Agg Files" << std::endl;
    dat.Read_Agg_Files();
    
    //Build the RVE w/forced agglomeration algorithm
    dat.Build_RVE();
    std::cout << uid << bid << ": RVE built" << std::endl;
    
    //Find the neighboring particles
    dat.Find_Neighbors();
    
    //FOR DEBUGGING
    //dat.CheckNeighbors();
    
    //Solve for conductivity
    bool is_net = dat.SolveForConductivity();
    
    if (is_net == false)
    {
        //std::cout << "No net found" << std::endl;
        return(0);
    }
    //Output results (OLD)
    //dat.OutputResults(fout);
    
    //Get final time
    auto end = std::chrono::high_resolution_clock::now();
    //std::cout << "Time elapsed 1st = " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;
    
    //Output Ovito File
    dat.Output_Ovito_File(fplot, uid, bid, cond_out);
    
    /*
    dat.Reset_Neighbors();
    //Loop through simulated strains and find strained conductivities
    for (int i = 0; i < strain_elong.size(); i++)
    {
        dat.Find_Strained_Neighbors(strain_elong[i]);
        bool is_net_s = dat.SolveForStrainConductivity(strain_elong[i], cond_out);
        if (is_net_s == false)
        {
            std::cout << "No net found, strain=" << strain_elong[i] << std::endl;
        }
        dat.Reset_Neighbors();
        std::cout << uid << bid << ": " << i << "/" << strain_elong.size() << std::endl;
    }
    */
    
    cond_out.close();
    //std::cout << "Placement Finished" << std::endl;
    //Close files
    //o_files.CloseFiles(fout, fplot, ferr);
}
