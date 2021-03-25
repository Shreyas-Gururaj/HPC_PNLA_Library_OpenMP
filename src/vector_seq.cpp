/**
 * @file vector_seq.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-03-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "vector_seq.h"
#include <cmath>

namespace pnla{


    /// Here is the place for the actual implementation of your sequential
    /// vector routines

    //
    void vector_init_constant_elements(vector_seq &obj_const_elem, const int vector_dimension, const double constant_value)

    {
        obj_const_elem.values.resize(vector_dimension);
        for(unsigned int elem_index = 0; elem_index < obj_const_elem.values.size(); elem_index++)
        {
            obj_const_elem.values[elem_index] = constant_value;
        }
    }

    //
    void vector_init_range_elements(vector_seq &obj_range_elem, const int vector_dimension)
    {
        obj_range_elem.values.resize(vector_dimension);
        for(unsigned int elem_index = 0; elem_index < obj_range_elem.values.size(); elem_index++)
        {
            obj_range_elem.values[elem_index] = static_cast<double>(elem_index);
        }
    }

    //
    void vector_init_std_doubles(vector_seq &obj_std_doubles, std::vector<double> &std_vector, const int vector_dimension)
    {
        obj_std_doubles.values.resize(vector_dimension);
        std_vector.resize(obj_std_doubles.values.size());
        for(unsigned int elem_index = 0; elem_index < obj_std_doubles.values.size(); elem_index++)
        {
            obj_std_doubles.values[elem_index] = std_vector[elem_index];
        }

    }

    //
    void vector_copy(const vector_seq &obj_initial, vector_seq &obj_to_be_copied)

    {
        obj_to_be_copied.values.resize(obj_initial.values.size());
        for(unsigned int elem_index = 0; elem_index < obj_initial.values.size(); elem_index++)
        {
            obj_to_be_copied.values[elem_index] = obj_initial.values[elem_index];
        }
    }

    //
    void vector_scale(vector_seq &obj_scale, const double scaling_factor)
    {
        for(unsigned int elem_index = 0; elem_index < obj_scale.values.size(); elem_index++)
        {
            obj_scale.values[elem_index] = obj_scale.values[elem_index] * scaling_factor;
        }
    }

    //
    double vector_dot_product(const vector_seq &obj_dot_x, const vector_seq &obj_dot_y)
    {
        double result_dot_prod = 0.0;

        for(unsigned int elem_index = 0; elem_index < obj_dot_x.values.size(); elem_index++)
            {
                result_dot_prod += obj_dot_x.values[elem_index] * obj_dot_y.values[elem_index];
            }

        return result_dot_prod;
    }

    //
    double vector_euclidean_norm(const vector_seq &obj_norm_x)
    {
        double result_euclidean_norm = sqrt(vector_dot_product(obj_norm_x, obj_norm_x));
        return result_euclidean_norm;

    }
   
   //
   void vector_scaled_addition(vector_seq &obj_scaled_add_y, const vector_seq &obj_scaled_add_x, const double scaling_factor_add)
   {
      for(unsigned int elem_index = 0; elem_index < obj_scaled_add_x.values.size(); elem_index++)
           {
               obj_scaled_add_y.values[elem_index] = obj_scaled_add_y.values[elem_index] + scaling_factor_add * obj_scaled_add_x.values[elem_index];
           }
   }

} // end namespace pnla