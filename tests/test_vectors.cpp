/**
 * @file test_vectors.cpp
 * @author Thomas Rau (thomas1.rau@uni-bayreuth.de)
 * @brief Bare bone structure of a small unit test for your vector routines
 * @version 0.1
 * @date 2021-03-14
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <limits> //epsilon()
#include "vector_seq.h"
#include <cmath>


//
template<typename Vector>
int test_const_norm(int test_sucess_count, const int dimension, const double epsilon)
{
    Vector const_norm_x;   //May be a better name with vector
    pnla::vector_init_constant_elements(const_norm_x, dimension, 1.0);
    double norm_function = pnla::vector_euclidean_norm(const_norm_x);
    double norm_target = sqrt(dimension);

    if(abs(norm_function - norm_target) > epsilon)
    {
        std::cout << "Either the constant element initialization or the Eucledian norm does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;

}

//
template<typename Vector>
int test_range_norm(int test_sucess_count, const int dimension, const double epsilon)
{
    Vector range_norm_y;
    pnla::vector_init_range_elements(range_norm_y, dimension);
    double norm_function = pnla::vector_euclidean_norm(range_norm_y);
    double norm_target = 0.0;

    for(int i = 0; i < dimension; i++)
    {
        norm_target += static_cast<double>(i * i);
    }

    norm_target = sqrt(norm_target);

    if(abs(norm_function - norm_target) > epsilon)
    {
        std::cout << "The constant element initialization does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;
}

//
template<typename Vector>
int test_copy_scaled_add(int test_sucess_count, const int dimension, const double epsilon)
{
    Vector range_y;
    Vector const_z;
    pnla::vector_init_range_elements(range_y, dimension);
    pnla::vector_init_constant_elements(const_z, dimension, 7.5);
    pnla::vector_copy(range_y, const_z);
    pnla::vector_scale(const_z, 5.0);
    pnla::vector_scaled_addition(range_y, const_z, -0.2);

    double norm_function = pnla::vector_euclidean_norm(range_y);

    if(abs(norm_function) > epsilon)
    {
        std::cout << "Either the copy_vector or the vector scaling or the scaled vector addition does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;

}

//
template<typename Vector>     //in the template function definition "Vector" is an alias for the struct
int test_vector_routines(const int dimension, const double epsilon)
{   
    int test_sucess = 0;
 
    //
    test_sucess += test_const_norm<Vector>(test_sucess, dimension, epsilon);

    //
    test_sucess += test_range_norm<Vector>(test_sucess, dimension, epsilon);

    //
    test_sucess += test_copy_scaled_add<Vector>(test_sucess, dimension, epsilon);

    return test_sucess;
}

/**
 * @brief Test Programm for pnla's vector structs/classes
 * 
 * @param argc Number of arguments  
 * @param argv Programm argument list
 * @return int On sucess = 0, if return value != 0 some test failed
 */
int main(int argc, char *argv[])
{
    int dim = 20;


 	if(argc == 2)
	{
         dim = std::stoi(argv[1]);
	}

    const double epsilon(std::numeric_limits<double>::epsilon()); 
    int test_result = 0; 

    std::cout<<"Test Sequential Vector"<<std::endl;
    /// call of test_vector should look something like this
    test_result = test_vector_routines<pnla::vector_seq>(dim, epsilon);

    if(test_result !=0 )
        return test_result;

    return test_result;
}



