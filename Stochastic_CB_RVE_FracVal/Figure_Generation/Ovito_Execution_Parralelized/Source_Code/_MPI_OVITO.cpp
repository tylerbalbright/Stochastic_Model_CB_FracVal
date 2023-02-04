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
    //std::string target_directory(argv[2]);
    std::string pixels(argv[2]);
    std::string XY(argv[3]);
    std::string CB(argv[4]);
    std::string Num_Threads(argv[5]);   
 
    std::string check_home_folder = "";
    std::size_t found; // = str.find_first_of("/");
    int last = 0;
    std::string truncated_folder_check = working_folder;

    for (int i = 0; i < 3; i++)
    {
        if (found!=std::string::npos)
        {
            found = truncated_folder_check.find_first_of("/");
            check_home_folder += truncated_folder_check.substr(0, found+1);
            truncated_folder_check = truncated_folder_check.substr(found + 1, std::string::npos);
        }
    }
    
    bool on_partition;
    if (check_home_folder != "/home/talbright/")
    {
        //std::cout << "Check Home Folder: " << check_home_folder << std::endl;
        //std::cout << "Not in home foolder, transfer first";
        on_partition = true;
    }
    else
    {
        //std::cout << "Not on partition" << std::endl;
        on_partition = false;
    }
    
    if ( myrank == 0 )
    {
        if (on_partition == true)
        {
            std::vector<std::string> file_names;    
	//have to move files intermittently to the partition before beginning to work on them! and then move them back.....
            std::string dir = working_folder + "/";
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
                }
                else
                {
                    continue;
                }
            }
            //Make temp directory. Delete old files.
            if (boost::filesystem::is_directory(boost::filesystem::path("_temp")) == false)
            {
                std::string command_mkdir = "mkdir _temp";
                std::system(command_mkdir.c_str());
            }
            else
            {
                std::string command_rmdir = "rm -rf _temp";
                std::system(command_rmdir.c_str());
                std::string command_mkdir = "mkdir _temp";
                std::system(command_mkdir.c_str());
            }
            //Make image folder in temp directory
            std::string command_mkdir = "mkdir _temp/_" + XY + "_" + CB + "_" + pixels;
            std::system(command_mkdir.c_str());

            std::vector<int> Files_On_Procs(file_names.size());
            
            int number_of_files = file_names.size();
            int n = 1;
            while (n <= number_of_files)
            {
                //Fill up the nodes to get started, send inital commands!
                while (n < numprocs && n <= number_of_files)
                {
                    std::string command = "cp " + working_folder + "/" + file_names[n-1] + " ./_temp/";
                    std::system(command.c_str());
                    bool a = false;
                    //First send a bool, false means we're sending another cmd
                    MPI_Send(&a, 1, MPI_CXX_BOOL, n, 0, MPI_COMM_WORLD);
                    MPI_Send(file_names[n-1].c_str(), file_names[n-1].size(), MPI_CHAR, n, 0, MPI_COMM_WORLD);
                    Files_On_Procs[n-1] = n-1;
                    n += 1;
                }
                int next;
                MPI_Recv(&next, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                std::string remove_file = "rm _temp/" + file_names[Files_On_Procs[next]];
                std::system(remove_file.c_str());
                std::string command = "cp " + working_folder + "/" + file_names[n-1] + " ./_temp/";
                std::system(command.c_str());
                Files_On_Procs[next] = n;
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
            std::vector<std::string> file_names;
            std::string cwd = getcwd(NULL,0);
            std::string dir = working_folder + "/";
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
                }
                else
                {
                    continue;
                }
            }
            
            //Check for picture folder
            std::string pic_folder = dir + "/_" + XY + "_" + CB + "_" + pixels + "/";
            boost::filesystem::path pic(pic_folder);
            if (boost::filesystem::is_directory(pic) == false)
            {
                std::string command_mkdir = "mkdir " + pic_folder;
                std::system(command_mkdir.c_str());
            }
            else
            {
                std::string command_rmdir = "rm -rf " + pic_folder;
                std::system(command_rmdir.c_str());
                std::string command_mkdir = "mkdir " + pic_folder;
                std::system(command_mkdir.c_str());
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
    }
    else
    {
        if (on_partition == true)
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
                    std::string cwd = getcwd(NULL,0);
                    std::string temp_dir = cwd + "/" + "_temp/";
                    std::string unpack = "xz -d -k " + temp_dir + file;
                    std::system(unpack.c_str());
                    std::string command = "./ovito-pro-3.4.3-x86_64/bin/ovitos --nthreads " + Num_Threads + " _OVITO_Image_Gen.py " + file  + " " + temp_dir + " " + pixels + " " + XY + " " + CB;
                    std::system(command.c_str());
                    MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    delete [] temp_char;
                    continue;
                }
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
                    std::string unpack = "xz -d -k " + working_folder + "/" + file;
                    //std::cout << "unpack: " << unpack <<std::endl;
                    std::system(unpack.c_str());
                    std::string command = "./ovito-pro-3.4.3-x86_64/bin/ovitos --nthreads " + Num_Threads + "  _OVITO_Image_Gen.py " + file  + " " + working_folder + " " + pixels + " " + XY + " " + CB;
                    std::system(command.c_str());
                    MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    delete [] temp_char;
                    continue;
                }
            }
        }
    }
    
    MPI_Finalize( );
    return 0;
}
