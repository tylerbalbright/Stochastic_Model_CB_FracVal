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
    //Initialization of MPI variables
    int  numprocs, myrank, namelen, i;
    
    //MPI Initialization stuffs. I borrowed a lot of code for this and it just works
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    char greeting[MPI_MAX_PROCESSOR_NAME + 80];
    char temp[MPI_MAX_PROCESSOR_NAME + 80];
    MPI_Status status;
    
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Get_processor_name( processor_name, &namelen );

    //From hello world example. Each process states its rank
    sprintf( greeting, "%d", myrank);
    
    //Processor of rank zero does all of the task assignment and splitting the load
    //between all other available processes
    if ( myrank == 0 )
    {
        //Define vector of file names and folder names to be processed
        std::vector<std::string> file_names;
        std::vector<std::string> folder_names;
        
        //Get the current working directory
        std::string cwd = getcwd(NULL,0);
        
        //The files to be compressed are located in the RVE_Data folder of the working directory
        std::string dir = cwd + "/RVE_Data/";
        std::cout << "dir = " << dir << std::endl;

        //Initialize a boost object for listing the contents of the directory
        boost::filesystem::path p(dir);
        
        //Check the data folder to see if it exists
        if (boost::filesystem::is_directory(p) == false)
        {
            std::cout << "ERROR: DATA FOLDER NOT FOUND, exiting" << std::endl;
            return(1);
        }
        
        //If the data folder exists, loop through its entries for processing
        for (auto i = boost::filesystem::directory_iterator(p); i != boost::filesystem::directory_iterator(); i++)
        {
            //Check to see if entry is a directory. If it is not a directory, skip it.
            if (!boost::filesystem::is_directory(i->path()))
            {
                continue;
            }
            else
            {
                //If the entry a directory not starting with "_" it needs processed
                if (i->path().filename().string()[0] != '_')
                {
                    //Add foler name to the list of foldernames
                    folder_names.push_back(i->path().filename().string());
                }
            }
        }
        
        //Loop through the each folder to be processed to compress the files it contains
        for (int i = 0; i < folder_names.size(); i++)
        {
            //Set the string of the folder name
            std::string next_dir = dir + folder_names[i] + "/";
            
            //Initialize the boost library using the folder name string
            boost::filesystem::path q(next_dir);
            
            //Loop through the entries in the folder to pull out the files to process
            for (auto j = boost::filesystem::directory_iterator(q); j != boost::filesystem::directory_iterator(); j++)
            {
                //If the entry is a directory, skip it.
                if (!boost::filesystem::is_directory(j->path()))
                {
                    //It is not a directory, so it is assumed to be a file to compress
                    //add the file name to the list to process
                    file_names.push_back(j->path().filename().string());
                }
            }
            
            //This is where the MPI task assigment is performed. In this function, we need to compress
            //all of the files in the file_names vector. To do so, the host process (rank 0) tells
            //another process (let's say process rank 1) to work on a file. It does this by sending
            //process rank 1 the filename of the file to compress. Later, process rank 1 sends a
            //message to process rank 0 when this task is complete, and then the next filename is
            //sent to process rank 1 until all files have been processed.
            int number_of_files = file_names.size();
            int n = 1;
            std::string current_file;
            
            //While loop to make sure all files get processed.
            while (n <= number_of_files)
            {
                //First we send initial tasks to the nodes to get started. Fill em all up!
                //The number of processors available is defined by num_procs
                while (n < numprocs && n <= number_of_files)
                {
                    //Boolean temp variable used to send in MPI commands
                    bool a = false;
                    
                    //Set the current file string
                    current_file = next_dir + file_names[n-1];
                    
                    //First send a boolean to the target process. A false boolean signifies that another
                    //command is going to be sent directly after.
                    MPI_Send(&a, 1, MPI_CXX_BOOL, n, 0, MPI_COMM_WORLD);
                    
                    //Now that the boolean has been sent, send the filename to the process
                    MPI_Send(current_file.c_str(), current_file.size(), MPI_CHAR, n, 0, MPI_COMM_WORLD);
                    
                    //Increment the variable n to move on to the next file to be processed
                    n += 1;
                }
                
                //Now, all of the availble processes have been tasked with compressing a file. In the
                //event that there are still more files to be processed than those sent, we wait until
                //a process is finished, and then we send it another file to work on. The MPI_Recv
                //command waits for any process to attmept to communicate to process rank 0.
                int next;
                MPI_Recv(&next, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                //First send a bool, false means we're sending another cmd. The variable "next" tells
                //us which process sent us the "finished" command.
                bool a = false;
                MPI_Send(&a, 1, MPI_CXX_BOOL, next, 0, MPI_COMM_WORLD);
                
                // We send the next filename to the process rank (next)
                current_file = next_dir + file_names[n-1];
                MPI_Send(current_file.c_str(), current_file.size(), MPI_CHAR, next, 0, MPI_COMM_WORLD);
                
                //Increment the filename counter
                n += 1;
            }
            //Clear the filenames list between iterations of target directories
            file_names.clear();
        }
        //After everything folder and file has been assigned, we need to wait for everything to finish
        //before exiting. To do so, we wait on responses from each process (expect for rank 0) in order.
        //Once all processes have responded, we can send commands to each process killing them, and then
        //exiting finally!
        
        //Loop through the processes
        for (int i = 1; i < numprocs; i++)
        {
            //Boolean set to true triggers the slave process exit command
            bool t = true;
            int not_needed;
            
            //Wait for the response from process i
            MPI_Recv(&not_needed, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            
            //Once the response is received, send the exit command
            MPI_Send(&t, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
        }
    }
    //Processes of rank>0 do the grunt work of compressing the files in a parallelized fashion
    else
    {
        //Set temporary integer that tells us when to exit the program
        int i = 0;
        
        //While exit command has yet to be recieved, wait to receive command from rank 0
        while (i == 0)
        {
            //Set boolean that tells us what the incoming command means
            bool type;
            
            //Recieve the boolean command from rank 0
            MPI_Recv(&type, sizeof( type ), MPI_CXX_BOOL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            
            //If command boolean is true, initiate exit of program
            if (type == true)
            {
                i = 1;
            }
            //If command boolean is false, receive the name of the file to be compressed.
            else
            {
                //Define temporary integer to count the number of characters to receive
                int num_char = 0;
                
                //Get the number of characters from the MPI commands
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_CHAR, &num_char);
                
                //Define an array of characters to read in the command from rank 0
                char *temp_char = new char [num_char+1]{};
                temp_char[num_char] = '\0';
                
                //Recieve the file name to be processed
                MPI_Recv(temp_char, num_char+1, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                //Convert the array of characters to a string and output the file name
                std::string file(temp_char);
                std::cout << "filename: " << file << std::endl;
                
                //Define the command string (called on the command line for compressing the file)
                std::string command = "xz -T0 " + file;
                
                //Execute the compression command via the command line
                std::system(command.c_str());
                
                //After the file is compressed, send the completion signal to process rank 0
                MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                
                //Delete the temporary character array
                delete [] temp_char;
                continue;
            }
        }
    }
    MPI_Finalize( );
    return 0;
}
