/**
 * @file HelloCluster.cpp
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief 
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <omp.h>
#include <mpi.h> 

/**
 * @brief HelloWorld again but in a fancy hybrid OpenMP/MPI version 
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char* argv [])
{
    // Every MPI Programm needs to start with this command, it initializes the MPI-World
    MPI_Init(&argc, &argv);

    int my_proc_ID = 0; 
    
    // Get's processors ID
    MPI_Comm_rank(MPI_COMM_WORLD, &my_proc_ID);
    
    int nr_procs = 0;
    // Get's total number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &nr_procs);
    
    // Get's name of machine/node 
    char my_name[MPI_MAX_PROCESSOR_NAME];
    int name_length = 0;
    MPI_Get_processor_name(my_name, &name_length);
    

    int current_rank = 0;
    while(current_rank < nr_procs)
    {
        if(my_proc_ID == current_rank)
        {
            // you know already this part
            # pragma omp parallel 
            {   
                int nr_total_threads = omp_get_num_threads();
                int my_thread_ID = omp_get_thread_num();

                #pragma omp critical
                {
                    std::cout<<"I'm thread " << my_thread_ID << "/" << nr_total_threads;
                    std::cout<<" of process nr " << my_proc_ID << " on node " << my_name <<std::endl; 
                }
            }           
        }
        current_rank++;
        
        // wait untill all processers arrived at this point -> Synchronisation
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // End MPI Processes
    MPI_Finalize();

    return 0;
}
