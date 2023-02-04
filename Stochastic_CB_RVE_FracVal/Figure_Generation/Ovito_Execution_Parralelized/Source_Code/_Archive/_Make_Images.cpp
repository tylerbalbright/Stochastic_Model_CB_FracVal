#include "mpi.h"
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
    int  numprocs, myrank, namelen, i;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    char greeting[MPI_MAX_PROCESSOR_NAME + 80];
    char temp[MPI_MAX_PROCESSOR_NAME + 80];
    MPI_Status status;
    
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Get_processor_name( processor_name, &namelen );

    sprintf( greeting, "%d", myrank);
    std::string working_folder(argv[1]);
    std::string pixels(argv[2]);
    
    if ( myrank == 0 )
    {
        std::vector<std::string> file_names;
        std::string cwd = getcwd(NULL,0);
        std::string dir = cwd + "/" + working_folder + "/";
        boost::filesystem::path p(dir);
        if (boost::filesystem::is_directory(p) == false)
        {
            std::cout << "ERROR: FOLDER NOT FOUND, exiting" << std::endl;
            return(1);
        }
        for (auto i = boost::filesystem::directory_iterator(p); i != boost::filesystem::directory_iterator(); i++)
        {
            //we eliminate directories in a list
            if (!boost::filesystem::is_directory(i->path()))
            {
                file_names.push_back(i->path().filename().string());
                std::cout << i->path().filename().string() << std::endl;
            }
            else
            {
                continue;
            }
        }
        
        int number_of_files = file_names.size();
        int n = 1;
        while (n <= number_of_files)
        {
            //Fill up the nodes to get started, send inital commands!
            while (n < numprocs && n <= number_of_files)
            {
                bool a = false;
                //First send a bool, false means we're sending another cmd
                MPI_Send(&a, 1, MPI_CXX_BOOL, n, 0, MPI_COMM_WORLD);
                MPI_Send(file_names[n-1].c_str(), file_names[n-1].size(), MPI_CHAR, n, 0, MPI_COMM_WORLD);
                n += 1;
            }
            int next;
            MPI_Recv(&next, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            bool a = false;
            //First send a bool, false means we're sending another cmd
            MPI_Send(&a, 1, MPI_CXX_BOOL, next, 0, MPI_COMM_WORLD);
            MPI_Send(file_names[n-1].c_str(), file_names[n-1].size(), MPI_CHAR, next, 0, MPI_COMM_WORLD);
            n += 1;
        }
        for (int i = 1; i < numprocs; i++)
        {
            bool t = true;
            MPI_Send(&t, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        int i = 0;
        while (i == 0)
        {
            bool type;
            MPI_Recv(&type, sizeof( type ), MPI_CXX_BOOL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (type == true)
            {
                i = 1;
            }
            else
            {
                int num_char = 0;
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_CHAR, &num_char);
                char *temp_char = new char [num_char+1]{};
                temp_char[num_char] = '\0';
                MPI_Recv(temp_char, num_char+1, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                std::string file(temp_char);
                std::cout << "filename: " << file << std::endl;
                std::string command = "/home/talbright/Visit_Source/bin/visit -nowin -verbose -cli -no-icet -s _Image_from_Vtk.py " + working_folder + " " + file + " " + pixels;
                std::system(command.c_str());
                MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                delete [] temp_char;
                continue;
            }
        }
    }
    
    MPI_Finalize( );
    return 0;
}
