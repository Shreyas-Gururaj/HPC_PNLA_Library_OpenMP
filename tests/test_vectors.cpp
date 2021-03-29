/**
 * @file test_vectors.cpp
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief 
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <limits> //epsilon()
#include "vector_seq.h"
#include "vector_omp.h"
#include <cmath>
#include <string>
#include <omp.h>
#include<chrono>


/**
 * @brief Checks for the correctness of the functions "vector_init_constant_elements" and "vector_euclidean_norm" simultaneously.
 * 
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq".
 * @param test_sucess_count initialized with 0.
 * @param dimension dimension of the vector.
 * @param epsilon tolerable error in computation.
 * @return int int Returns either 0 or 1.
 * 
 */
template<typename Vector>
int test_const_norm(int test_sucess_count, const int dimension, const double epsilon)
{
    Vector const_x;

    auto start_time = std::chrono::high_resolution_clock::now();
    pnla::vector_init_constant_elements(const_x, dimension, 1.0);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time in milliseconds for vector_init_constant_elements is :  " << run_time.count() << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    double norm_function = pnla::vector_euclidean_norm(const_x);
    end_time = std::chrono::high_resolution_clock::now();
    run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time in milliseconds for vector_euclidean_norm is :  " << run_time.count() << std::endl;


    double norm_target = sqrt(dimension);   // Norm of vectors of all elemnts = 1.0 is nothing but the sqrt(dimension of the vector).

    if(abs(norm_function - norm_target) > epsilon)
    {
        std::cout << "Either the constant element initialization or the Eucledian norm does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;

}

/**
 * @brief Checks the correctness of the functions "vector_init_range_elements" and "vector_euclidean_norm" simultaneously.
 * 
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq".
 * @param test_sucess_count initialized with 0.
 * @param dimension dimension of the vector.
 * @param epsilon tolerable error in computation.
 * @return int int Returns either 0 or 1.
 * 
 */
template<typename Vector>
int test_range_norm(int test_sucess_count, const int dimension, const double epsilon)
{
    Vector range_y;

    auto start_time = std::chrono::high_resolution_clock::now();
    pnla::vector_init_range_elements(range_y, dimension);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time in milliseconds for vector_init_range_elements is :  " << run_time.count() << std::endl;

    double norm_function = pnla::vector_euclidean_norm(range_y);
    double norm_target = 0.0;

    for(int i = 0; i < dimension; i++)
    {
        norm_target += static_cast<double>(i * i);  // Norm of vectors containing range of dimension is nothing but the sum(i*i);
    }

    norm_target = sqrt(norm_target);

    if(abs(norm_function - norm_target) > epsilon)
    {
        std::cout << "The constant element initialization does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;
}


/**
 * @brief Checks the correctness of the functions "vector_copy", "vector_scale" and "vector_scaled_addition" simultaneously.
 * 
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq"
 * @param test_sucess_count initialized with 0.
 * @param dimension dimension of the vector.
 * @param epsilon tolerable error in computation.
 * @return int Returns either 0 or 1.
 * 
 */
template<typename Vector>
int test_copy_scaled_add(int test_sucess_count, const int dimension, const double epsilon)
{
    Vector range_y;
    Vector const_z;
    pnla::vector_init_range_elements(range_y, dimension);
    pnla::vector_init_constant_elements(const_z, dimension, 7.5);

    auto start_time = std::chrono::high_resolution_clock::now();
    pnla::vector_copy(range_y, const_z);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time in milliseconds for vector_copy is :  " << run_time.count() << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    pnla::vector_scale(const_z, 5.0);         // Scaling by a factor α
    end_time = std::chrono::high_resolution_clock::now();
    run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time in milliseconds for vector_scale is :  " << run_time.count() << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    pnla::vector_scaled_addition(range_y, const_z, -0.2);  // Rescaling by a factor 0f -(1/α) and adding with the vector before scaling.
    end_time = std::chrono::high_resolution_clock::now();
    run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time in milliseconds for vector_scaled_addition is :  " << run_time.count() << std::endl;

    double norm_function = pnla::vector_euclidean_norm(range_y);

    if(abs(norm_function) > epsilon)
    {
        std::cout << "Either the copy_vector or the vector scaling or the scaled vector addition does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;

}

/**
 * @brief Appends the value of the test_sucess initialized with 0 by a positive integer if any of the test fails.
 * 
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq"
 * @param dimension initialized with 0.
 * @param epsilon tolerable error in computation.
 * @return int Returns either 0 or a positive integer.
 * 
 */
template<typename Vector>
int test_vector_routines(const int dimension, const double epsilon)
{   
    int test_sucess = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    test_sucess += test_const_norm<Vector>(test_sucess, dimension, epsilon);

    test_sucess += test_range_norm<Vector>(test_sucess, dimension, epsilon);

    test_sucess += test_copy_scaled_add<Vector>(test_sucess, dimension, epsilon);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Total time in milliseconds is :  " << run_time.count() << std::endl;

    return test_sucess;
}

/**
 * @brief Test Programm for pnla's vector structs/classes.
 * 
 * @param argc Number of arguments.
 * @param argv Programm argument list.
 * @return int On sucess = 0, if return value != 0 some test failed.
 * 
 */
int main(int argc, char *argv[])
{
    int dim = 20;


 	if(argc == 2)
	{
        dim = std::stoi(argv[1]);
	}

    if(argc == 3)
	{
        dim = std::stoi(argv[1]);
        const int nr_of_threads = std::stoi(argv[2]);
        omp_set_num_threads(nr_of_threads);
	}
    //int nr_of_threads = 1;

    const double epsilon(std::numeric_limits<double>::epsilon()); 
    int test_result = 0;


    std::cout<<"Test Sequential Vector"<<std::endl;

    // instantiating the template function with the struct "vector_seq".
    test_result = test_vector_routines<pnla::vector_seq>(dim, epsilon);

    std::cout<<"Test OMP Vector"<<std::endl;

    // instantiating the template function with the struct "vector_seq".
    test_result = test_vector_routines<pnla::vector_omp>(dim, epsilon);

    if(test_result !=0 )
        return test_result;

    return test_result;
}



