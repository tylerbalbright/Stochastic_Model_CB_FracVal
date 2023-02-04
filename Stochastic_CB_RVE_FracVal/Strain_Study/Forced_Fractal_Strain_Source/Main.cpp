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

//Definition of global simulation variables to be read via command line input
int cmd_PctRandLower;
int cmd_PctRandUpper;
int cmd_XY;
int cmd_thickness;
int cmd_wcb;
int cmd_fill;
std::string uid;
std::string bid;

//Main function
int main(int argc, char *argv[])
{
    //Read command line input variables into global containers
    cmd_PctRandLower = atoi(argv[1]);
    cmd_PctRandUpper = atoi(argv[2]);
    cmd_XY = atoi(argv[3]);
    cmd_thickness = atoi(argv[4]);
    uid = argv[5];
    bid = argv[6];
    cmd_wcb = atoi(argv[7]);
    cmd_fill = atoi(argv[8]);
    
    //Define vector of floats of uniaxial strain elongations to simulate on generated RVE
    std::vector<float> strain_elong{};
    
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    //Define file streams, some are used for error checking and validation.
    //Some streams may not be used until sections of code are uncommented.
    std::ifstream fin;      //stream for problem input
    std::ofstream fout;      //stream for problem output
    std::ofstream ferr;      //stream for error output
    std::ofstream fplot;     //stream for plotting RVE
    std::ofstream cond_out;  //stream for reporting conductivities before/after strain

    //Definition of RVEdata object. This object houses a majority of the functions needed to populate
    //and test RVEs for conductivity/piezoresistivity
    RVEdata dat;        //object to store problem data
    
    //Defintion of Files object. This object houses functions pertaining to the opening/saving of streamed files.
    Files o_files;      //object to manipulate files
    
    //Open all file streams
    o_files.OpenFiles(fin, fout, fplot, ferr);
    
    //Read in data for problem
    dat.ReadInputData(fin, ferr, cmd_PctRandLower, cmd_PctRandUpper, cmd_XY, cmd_thickness, cmd_wcb, cmd_fill);
    //Close input file
    fin.close();

    //Read in FracVal aggregate files in the execution folder
    dat.Read_Agg_Files();
    
    //Build the RVE w/forced agglomeration algorithm
    dat.Build_RVE();
    std::cout << uid << bid << ": RVE built" << std::endl;
    
    //Find the neighboring particles
    dat.Find_Neighbors();
    
    //Function for outputting coordinates and algorithmically identified neighboring pairs for
    //validating the code is behaving as expected.
    //dat.CheckNeighbors();
    
    //Solve for conductivity and return boolean to note whether a conductive networks exists or not.
    bool is_net = dat.SolveForConductivity();
    
    //If conductive network exists, run the strain tests. If not, exit the program.
    if (is_net == false)
    {
        return(0);
    }
    
    //Function to output conductivity data to save file (OUTDATED)
    //dat.OutputResults(fout);
    
    //Get final time for computing total simulation time
    auto end = std::chrono::high_resolution_clock::now();
    
    //Function to output coordinate file of RVE elements
    dat.Output_Ovito_File(fplot, uid, bid, cond_out);
    
    /*****************************************************************/
    /************************ Strain tests ***************************/
    /*****************************************************************/

    //Call function to clear containers storing neighboring relationship data
    dat.Reset_Neighbors();
    //Loop through strain_elong vector while computing strained conductivities at each strain state
    for (int i = 0; i < strain_elong.size(); i++)
    {
        //Find neighboring relationships under new strain condition applied via directional scaling factors
        dat.Find_Strained_Neighbors(strain_elong[i]);
        //Function to determine if a conductive network exists at the simulated strained state
        bool is_net_s = dat.SolveForStrainConductivity(strain_elong[i], cond_out);
        //If no network found, say so via terminal. Don't exit loop because it complicates post-processing
        if (is_net_s == false)
        {
            std::cout << "No net found, strain=" << strain_elong[i] << std::endl;
        }
        //After conductivity at strain state has been found and recorded, reset containers containing neighbor
        //relationship data and repeat the process for the next strain state.
        dat.Reset_Neighbors();
        std::cout << uid << bid << ": " << i << "/" << strain_elong.size() << std::endl;
    }
    cond_out.close();
    //Function to close all the opened files (OUTDATED)
    //o_files.CloseFiles(fout, fplot, ferr);
}
