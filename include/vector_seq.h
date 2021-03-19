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
#define __VECTOR_SEQ_H__

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
            int get_vector_dimension() 
            { 
             return vector_dimension; 
            }

            //
            std::vector <double> get_init_vector()
            {
             return x;
            }
            


            //
            void inititalize_constant_elements(std::vector<double> &init_vector, const double constant_value);

            //
            void inititalize_range_elements(std::vector<double> &init_vector);

            //
            //gitvoid inititalize_std_doubles(std::vector<double> &init_vector);

        private:

            //
            unsigned int vector_dimension;

            //
            std::vector <double> x;


    };


}//end namespace pnla

#endif // __VECTOR_SEQ_H__