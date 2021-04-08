

#include <iostream>
#include <omp.h> // You need to include the OpenMP header
// Rembember to add the "-fopenmp" to your compiler options


int main(int argc, char *argv[])
{
    # pragma omp parallel 
    {   
        int nr_total_threads = omp_get_num_threads();
        int my_ID = omp_get_thread_num();

        #pragma omp critical
        {
            std::cout<<"Hello from " << my_ID << "/" << nr_total_threads << std::endl;
        }
    }  

    return 0;
}