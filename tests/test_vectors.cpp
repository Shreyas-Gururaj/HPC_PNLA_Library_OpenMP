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
#include <vector_seq.h>
#include <cmath>


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
template<typename Vector>
int test_vector_routines(const int dimension, const double epsilon)
{   
    int test_sucess = 0;
    
    test_sucess += test_const_elem_vector(0, dimension, epsilon);
    test_sucess += test_range_elem_vector(0, dimension, epsilon);

    /// Here you have to implement testing routines for pnla's vector 
    /// structures/classes. If your test fails, set test success to a non zero value;

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



