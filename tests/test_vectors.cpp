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


const double initial_values = 5.0;


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
/*
//
int test_const_elem_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        const double const_value = 4.0;

        pnla::vector_seq obj_const_elem(dimension, initial_values);
        obj_const_elem.vector_init_constant_elements(obj_const_elem, const_value);

        for(unsigned int i =0; i < obj_const_elem.values.size(); i++)
        {
            if(abs(obj_const_elem.values[i] - const_value) < epsilon)
            {
                continue;                                                  //Avoid using continue
            }
            else
            {
                test_sucess_count += 1;
            }
        }
        
        return test_sucess_count;
    }

//

int test_range_elem_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::vector_seq obj_range_elem(dimension, initial_values);
        obj_range_elem.vector_init_range_elements(obj_range_elem);

        for(unsigned int i =0; i < obj_range_elem.values.size(); i++)
        {
            if(abs(obj_range_elem.values[i] - i) < epsilon)
            {
                continue;
            }
            else
            {
                test_sucess_count += 1;
            }
        }
        
        return test_sucess_count;
    }

//
int test_std_double_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::vector_seq obj_std_double(dimension, initial_values);
        obj_std_double.vector_init_std_doubles(obj_std_double);
        std::vector<double> std_vector(obj_std_double.values.size());

        for(unsigned int i =0; i < obj_std_double.values.size(); i++)
        {
            if(abs(obj_std_double.values[i] - std_vector[i]) < epsilon)
            {
                continue;
            }
            else
            {
                test_sucess_count += 1;
            }
        }
        
        return test_sucess_count;
    }

//
int test_copy_vector(int test_sucess_count, const int dimension, const double epsilon)
    {

        pnla::vector_seq obj_initial(dimension, initial_values);
        pnla::vector_seq obj_copy(dimension, 7.5);
        obj_copy.vector_copy(obj_initial, obj_copy);

        for(unsigned int i =0; i < obj_initial.values.size(); i++)
        {
            if(abs(obj_copy.values[i] - obj_initial.values[i]) < epsilon)
            {
                continue;
            }
            else
            {
                test_sucess_count += 1;
            }
        }
        
        return test_sucess_count;
    }

//
int test_scale_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        const double scaling_factor = 9.6;
        pnla::vector_seq obj_scale(dimension, initial_values);
        pnla::vector_seq obj_target(dimension, initial_values);
        obj_scale.vector_scale(obj_scale, scaling_factor);
        std::vector<double> std_vector(obj_scale.values.size());

        for(unsigned int i =0; i < obj_scale.values.size(); i++)
        {
            if(abs((obj_scale.values[i]/scaling_factor) - obj_target.values[i]) < epsilon)
            {
                continue;
            }
            else
            {
                test_sucess_count += 1;
            }
        }
        
        return test_sucess_count;
    }

//
int test_eucledian_norm_vector(int test_sucess_count, const int dimension, const double epsilon)
    {
        pnla::vector_seq obj_norm(dimension, 1.0);
        const double norm_target = sqrt(dimension);                                               
        double norm_from_function = obj_norm.vector_dot_product(obj_norm, obj_norm);

        if(abs((norm_target - norm_from_function) < epsilon))
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
        pnla::vector_seq obj_scaled_add_x(dimension, 0.0);
        pnla::vector_seq obj_scaled_add_y(dimension, 0.0);                                               
        double norm_from_function = obj_scaled_add_x.vector_dot_product(obj_scaled_add_x, obj_scaled_add_y);

        if(abs((norm_from_function - 0.0) < epsilon))
        {
            return test_sucess_count;
        }
        else
        {
            return test_sucess_count += 1;
        }
    }
*/

//
int test_const_norm(int test_sucess_count, const int dimension, const double epsilon)
{
    pnla::vector_seq const_norm_x;
    pnla::vector_init_constant_elements(const_norm_x, dimension, 1.0);
    double norm_function = pnla::vector_euclidean_norm(const_norm_x);
    double norm_target = sqrt(dimension);

    if(abs(norm_function - norm_target) > epsilon)
    {
        return test_sucess_count + 1;
        std::cout << "Either the constant element initialization or the Eucledian norm does not work as intended";
    }

    return test_sucess_count;

}

//
int test_range_norm(int test_sucess_count, const int dimension, const double epsilon)
{
    pnla::vector_seq const_norm_y;
    pnla::vector_init_range_elements(const_norm_y, dimension);
    double dot_function = pnla::vector_euclidean_norm(const_norm_y);
    double dot_target = 0.0;

    for(int i = 0; i < dimension; i++)
    {
        
    }

    if(abs(dot_function - dot_target) > epsilon)
    {
        return test_sucess_count + 1;
        std::cout << "Either the constant element initialization or the Eucledian norm does not work as intended";
    }

    return test_sucess_count;
}

//
template<typename Vector>     //in the template function definition "Vector" is an alias for the struct
int test_vector_routines(const int dimension, const double epsilon)
{   
    int test_sucess = 0;
 
    //
    test_sucess += test_const_norm(test_sucess, dimension, epsilon);



    /// Here you have to implement testing routines for pnla's vector 
    /// structures/classes. If your test fails, set test success to a non zero value;

    /*int test_sucess = 0;
    const int test_count = 0;
    std::cout << dimension << std::endl;
    test_sucess += test_const_elem_vector(test_count, dimension, epsilon);
    test_sucess += test_range_elem_vector(test_count, dimension, epsilon);
    test_sucess += test_std_double_vector(test_count, dimension, epsilon);
    test_sucess += test_copy_vector(test_count, dimension, epsilon);
    test_sucess += test_scale_vector(test_count, dimension, epsilon);
    test_sucess += test_eucledian_norm_vector(test_count, dimension, epsilon);
    test_sucess += test_scaled_addition_vector(test_count, dimension, epsilon);*/

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

    // Just for illustration of template function and to surpress warnings
    //test_result = test_vector_routines<double>(dim, epsilon);

    if(test_result !=0 )
        return test_result;

    return test_result;
}



