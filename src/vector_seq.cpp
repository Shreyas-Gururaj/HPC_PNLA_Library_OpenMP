/**
 * @file vector_seq.cpp
 * @author Thomas Rau (thomas1.rau@uni-bayreuth.de)
 * @brief Barebone structure of a pnla source file
 * @version 0.1
 * @date 2021-03-14
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <vector_seq.h>
#include<cmath>

namespace pnla{


    /// Here is the place for the actual implementation of your sequential
    /// vector routines

    //
    void Vector_n_dim::vector_init_constant_elements(std::vector<double> &init_vector, const double constant_value)

    {
        for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            init_vector[elem_index] = constant_value;
        }
    }

    //
    void Vector_n_dim::vector_init_range_elements(std::vector<double> &init_vector)
    {
        for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            init_vector[elem_index] = elem_index;
        }
    }

    //
    void Vector_n_dim::vector_init_std_doubles(std::vector<double> &init_vector)
    {
        std::vector<double> std_vector(init_vector.size());
        for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            init_vector[elem_index] = std_vector[elem_index];
        }

    }

    //
    void Vector_n_dim::vector_copy(std::vector<double> &init_vector, std::vector<double> &to_copy_vector)

    {
        if(init_vector.size() == to_copy_vector.size())
        {
            for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
            {
                to_copy_vector[elem_index] = init_vector[elem_index];
            }
        }

         else
        {
            std::cout << "Dimensions of the vectors do not match and cannot be copied" << std::endl;
        }
    }

    //
    void Vector_n_dim::vector_scale(std::vector<double> &init_vector, const double scaling_factor)
    {
        for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            init_vector[elem_index] = init_vector[elem_index] * scaling_factor;
        }
    }

    //
    void Vector_n_dim::vector_dot_product(std::vector<double> &init_vector, std::vector<double> &vector_for_dot_prod)
    {
        double result_dot_prod = 0;

        if(init_vector.size() == vector_for_dot_prod.size())
        {
            for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
            {
                result_dot_prod += init_vector[elem_index] * vector_for_dot_prod[elem_index];
            }

            std::cout << "The dot product of vector x and vector y is" << " : " << result_dot_prod;
        }

        else
        {
            std::cout << "Dimensions of the vectors do not match and dot product cannot be performed" << std::endl;
        }
    }

    //
    void Vector_n_dim::vector_euclidean_norm(std::vector<double> &init_vector)
    {
        double squared_sum = 0;
        double result_euclidean_norm;

        for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            squared_sum += init_vector[elem_index] * init_vector[elem_index];
            result_euclidean_norm = sqrt(squared_sum);
        }

        std::cout << "The Euclidean norm is" << " : " << result_euclidean_norm;

    }
   
   //
   void Vector_n_dim::vector_scaled_addition(std::vector<double> &init_vector, std::vector<double> &to_scale_vector, const double scaling_factor_add)
   {
       if(init_vector.size() == to_scale_vector.size())
       {
           for(unsigned int elem_index = 0; elem_index < init_vector.size(); elem_index++)
           {
               to_scale_vector[elem_index] = to_scale_vector[elem_index] + scaling_factor_add * init_vector[elem_index];
           }
       }
    
   }

} // end namespace pnla