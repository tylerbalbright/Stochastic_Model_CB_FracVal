#include <stdio.h>
#include <string.h>
#include <random>
#include <ctime>
#include <iostream>
#include <boost/filesystem.hpp>
#include <unistd.h>
#include <vector>
#include <chrono>
#include <thread>

int main( int argc, char *argv[] )
{
    /*
    int  numprocs, myrank, namelen, i;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    char greeting[MPI_MAX_PROCESSOR_NAME + 80];
    char temp[MPI_MAX_PROCESSOR_NAME + 80];
    MPI_Status status;
    
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Get_processor_name( processor_name, &namelen );
     */
    
    std::vector<std::string> file_names;
    std::vector<std::string> folder_names;
    std::string cwd = getcwd(NULL,0);
    std::string dir = cwd + "/Training_Data/";
    boost::filesystem::path p(dir);
    if (boost::filesystem::is_directory(p) == false)
    {
        std::cout << "ERROR: TRAINING DATA FOLDER NOT FOUND, exiting" << std::endl;
        return(1);
    }
    for (auto i = boost::filesystem::directory_iterator(p); i != boost::filesystem::directory_iterator(); i++)
    {
        //we eliminate directories in a list
        if (!boost::filesystem::is_directory(i->path()))
        {
            continue;
        }
        else
        {
            if (i->path().filename().string()[0] != '_')
            {
                std::cout << i->path().filename().string();
                folder_names.push_back(i->path().filename().string());
            }
        }
    }
}
