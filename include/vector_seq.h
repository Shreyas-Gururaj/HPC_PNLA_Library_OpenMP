/**
 * @file vector_seq.h
 * @author Thomas Rau (thomas1.rau@uni-bayreuth.de)
 * @brief Barebone structure of a pnla header file
 * @version 0.1
 * @date 2021-03-14
 * 
 * @copyright Copyright (c) 2021
 * 
 */

/// Remember, every header needs it's guard !
#ifndef __VECTOR_SEQ_H__
#define __VECTOR_SEQ_H__S

#include <iostream>
#include <vector>
#include <string>

using namespace std;

/// put everything of the PNLA-library into the suitable namespace
namespace pnla{

    /// Start here with the defenition of your struct/class for storing a vector 
    /// of doubles
    
    class vector_seq      //look for a better name for the class.
    {
        public:

            //
            unsigned int vector_dimension;

            //
            std::vector<double> values; 

            //
            vector_seq() = delete;

            
            vector_seq(const unsigned int dimension, const double initial_values)
            {
                vector_dimension= dimension;
                values.assign(vector_dimension, initial_values);
            }

            //
            void vector_init_constant_elements(vector_seq &obj_const_elem, const double constant_value);

            //
            void vector_init_range_elements(vector_seq &obj_range_elem);

            //
            void vector_init_std_doubles(vector_seq &obj_std_doubles);

            //
            void vector_copy(const vector_seq &obj_initial, vector_seq &obj_to_be_copied);

            //
            void vector_scale(vector_seq &obj_scale, const double scaling_factor);

            //
            double vector_dot_product(const vector_seq &obj_dot_x, const vector_seq &obj_dot_y);

            //
            double vector_euclidean_norm(const vector_seq &obj_norm_x);

            //
            void vector_scaled_addition(vector_seq &obj_scaled_add_y, const vector_seq &obj_scaled_add_x, const double scaling_factor_add);

           //look for a better name for the vector.


    };


}//end namespace pnla
 
#endif // __VECTOR_SEQ_H__