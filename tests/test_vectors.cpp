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
#include <string>
#include <limits> //epsilon()
#include "vector_seq.h"
#include <cmath>

using namespace std;


/**
 * @brief Bare bone function template for testing. The beauty of templates is, that you
 *        have to write this test only once and then it should work for all of your 
 *        vector versions 
 * 
 * @tparam Vector 
 * @param dimension 
 * @param epsilon 
 * @return int 
 */

//
int test_const_elem_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        const double const_value = 4.0;
        
        v1.vector_init_constant_elements(z, const_value);
        if(abs(z[2] - const_value) < epsilon)
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_range_elem_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        
        v1.vector_init_range_elements(z);
        if(abs(z[2] - 2) < epsilon)
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_std_double_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        std::vector<double> std_double(z.size());
        v1.vector_init_std_doubles(z);
        if(abs(z[2] - std_double[2]) < epsilon)
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_copy_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        std::vector<double> y = v1.vector_copy(z);
        if(abs(z[2] - y[2]) < epsilon)
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_scale_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        const double scale_factor = 6.5;
        v1.vector_scale(z, scale_factor);

        if(abs((z[2]/scale_factor) - 1.0) < epsilon)
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_dot_product_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        std::vector<double> y(dimension, 30.0);
        const double dot_product_calculated = dimension * y[0] * z[0];
        double dot_product_from_function = v1.vector_dot_product(z, y);

        if(abs((dot_product_calculated - dot_product_from_function) < epsilon))
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_eucledian_norm_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        const double norm_calculated = sqrt(dimension);
        double norm_from_function = v1.vector_euclidean_norm(z);

        if(abs((norm_calculated - norm_from_function) < epsilon))
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
int test_scaled_addition_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::Vector_n_dim v1(dimension);
        std::vector<double> z = v1.get_init_vector();
        std::vector<double> y (dimension, 1.0);
        const double scale_factor = 6.5;
        const double result_scaled_addition = y[2] + (scale_factor * z[2]);

        v1.vector_scaled_addition(z, y, scale_factor);

        if(abs((result_scaled_addition - y[2])) < epsilon)
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }

//
template<typename Vector>
int test_vector_routines(const int dimension, const double epsilon)
{   
    /// Here you have to implement testing routines for pnla's vector 
    /// structures/classes. If your test fails, set test success to a non zero value;

    int test_sucess = 0;
    const int test_count = 0;
    std::cout << dimension << std::endl;
    test_sucess += test_const_elem_vector(test_count, dimension, epsilon);
    test_sucess += test_range_elem_vector(test_count, dimension, epsilon);
    test_sucess += test_std_double_vector(test_count, dimension, epsilon);
    test_sucess += test_copy_vector(test_count, dimension, epsilon);
    test_sucess += test_scale_vector(test_count, dimension, epsilon);
    test_sucess += test_dot_product_vector(test_count, dimension, epsilon);
    test_sucess += test_eucledian_norm_vector(test_count, dimension, epsilon);
    test_sucess += test_scaled_addition_vector(test_count, dimension, epsilon);

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
    test_result = test_vector_routines<pnla::Vector_n_dim>(dim, epsilon);

    // Just for illustration of template function and to surpress warnings
    //test_result = test_vector_routines<double>(dim, epsilon);

    if(test_result !=0 )
        return test_result;

    return test_result;
}



