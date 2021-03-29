/**
 * @file vector_seq.cpp
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief Implementation of all the routines for omp vectors.
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "vector_omp.h"
#include <cmath>
#include <omp.h>

namespace pnla{


    /// The actual implementation of your omp vector routines

    void vector_init_constant_elements(vector_omp &obj_const_elem, const int vector_dimension, const double constant_value)

    {
        ///resize vector
        obj_const_elem.values.resize(vector_dimension);
        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < obj_const_elem.values.size(); elem_index++)
        {
            obj_const_elem.values[elem_index] = constant_value;
        }
    }

    void vector_init_range_elements(vector_omp &obj_range_elem, const int vector_dimension)
    {
        ///resize vector
        obj_range_elem.values.resize(vector_dimension);
        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < obj_range_elem.values.size(); elem_index++)
        {
            obj_range_elem.values[elem_index] = static_cast<double>(elem_index);
        }
    }

    void vector_init_std_doubles(vector_omp &obj_std_doubles, std::vector<double> &std_vector, const int vector_dimension)
    {
        ///resize vector
        obj_std_doubles.values.resize(vector_dimension);
        std_vector.resize(obj_std_doubles.values.size());
        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < obj_std_doubles.values.size(); elem_index++)
        {
            obj_std_doubles.values[elem_index] = std_vector[elem_index];
        }

    }

    void vector_copy(const vector_omp &obj_initial, vector_omp &obj_to_be_copied)

    {
        ///resize vector
        obj_to_be_copied.values.resize(obj_initial.values.size());
        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < obj_initial.values.size(); elem_index++)
        {
            obj_to_be_copied.values[elem_index] = obj_initial.values[elem_index];
        }
    }

    void vector_scale(vector_omp &obj_scale, const double scaling_factor)
    {
        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < obj_scale.values.size(); elem_index++)
        {
            obj_scale.values[elem_index] = obj_scale.values[elem_index] * scaling_factor;
        }
    }

    double vector_dot_product(const vector_omp &obj_dot_x, const vector_omp &obj_dot_y)
    {
        double result_dot_prod = 0.0;

        #pragma omp parallel for shared(obj_dot_x, obj_dot_y) reduction(+:result_dot_prod)
        for(unsigned int elem_index = 0; elem_index < obj_dot_x.values.size(); elem_index++)
            {
                result_dot_prod += obj_dot_x.values[elem_index] * obj_dot_y.values[elem_index];
            }

        return result_dot_prod;
    }

    double vector_euclidean_norm(const vector_omp &obj_norm_x)
    {
        double result_euclidean_norm = sqrt(vector_dot_product(obj_norm_x, obj_norm_x));
        return result_euclidean_norm;

    }
   
   void vector_scaled_addition(vector_omp &obj_scaled_add_y, const vector_omp &obj_scaled_add_x, const double scaling_factor_add)
   {
      #pragma omp parallel for shared(obj_scaled_add_y, obj_scaled_add_x) schedule(static)
      for(unsigned int elem_index = 0; elem_index < obj_scaled_add_x.values.size(); elem_index++)
           {
               obj_scaled_add_y.values[elem_index] = obj_scaled_add_y.values[elem_index] + scaling_factor_add * obj_scaled_add_x.values[elem_index];
           }
   }

} // end namespace pnla