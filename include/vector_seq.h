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
    
    class Vector_n_dim
    {
        public:

            //
            Vector_n_dim() = delete;

            
            Vector_n_dim(const unsigned int dimension)
            {
                vector_dimension= dimension;
                x.assign(vector_dimension, 1.0);
            }

            //
            unsigned int get_vector_dimension() 
            { 
             return vector_dimension; 
            }

            //
            std::vector <double> get_init_vector()
            {
             return x;
            }


            //
            void vector_init_constant_elements(std::vector<double> &init_vector, const double constant_value);

            //
            void vector_init_range_elements(std::vector<double> &init_vector);

            //
            void vector_init_std_doubles(std::vector<double> &init_vector);

            //
            std::vector<double> vector_copy(std::vector<double> &init_vector);

            //
            void vector_scale(std::vector<double> &init_vector, const double scaling_factor);

            //
            double vector_dot_product(std::vector<double> &init_vector, std::vector<double> &vector_for_dot_prod);

            //
            void vector_euclidean_norm(std::vector<double> &init_vector);

            //
            void vector_scaled_addition(std::vector<double> &init_vector, std::vector<double> &to_scale_vector, const double scaling_factor_add);

        private:

            //
            unsigned int vector_dimension;

            //
            std::vector <double> x;


    };


}//end namespace pnla

#endif // __VECTOR_SEQ_H__