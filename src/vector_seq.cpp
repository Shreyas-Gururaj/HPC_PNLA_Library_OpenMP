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

namespace pnla{


    /// Here is the place for the actual implementation of your sequential 
    /// vector routines 

    void Vector_n_dim::inititalize_constant_elements(std::vector<double> &init_vector, const double constant_value)

    {
        for(int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            init_vector[elem_index] = constant_value;
        }
    }

    //
    void Vector_n_dim::inititalize_range_elements(std::vector<double> &init_vector)
    {
        for(int elem_index = 0; elem_index < init_vector.size(); elem_index++)
        {
            init_vector[elem_index] = elem_index;
        }
    }

    //
    // void Vector_n_dim::inititalize_std_doubles(std::vector<double> &init_vector)
    // {

    // }

} // end namespace pnla