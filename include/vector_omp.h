/**
 * @file vector_omp.h
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief Header files for vector_omp functions requird for the PCG algorithm.
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

/// Remember, every header needs it's guard !
#ifndef __VECTOR_OMP_H__
#define __VECTOR_OMP_H__

#include <vector>
#include <memory>
#include<utility>

/// put everything of the PNLA-library into the suitable namespace.

/**
 * @brief The struct contains "vector_dimension" to initialize a std vector of doubles of given dimension in 3 different ways i.e, 
 * 1. Vector containing constant elements.
 * 2. Vector containing range of elements.
 * 3. Vector containing elements from std vector of doubles given as argument.
 */
namespace pnla{

    /// Start here with the defenition of struct/class for storing a vector of doubles.
    
    struct vector_omp        
    {
        //Unsigned integer as dimension of the vector is always a positive integer.
        unsigned int vector_dimension;

        //Std vector of doubles to store the elemnts in the vector of given dimension.
        //std::unique_ptr<double[]> values;
        std::shared_ptr<double[]> values;

    };

        /**
         * @brief Takes the dimension and a double of "constsnt_value" as arguments to initialize
         *        our vector instance with the constant value as all the elements.
         * 
         * @param obj_const_elem Instance of the struct to initialize vector of a given dimension in place.
         * @param vector_dimension The required dimension of the initialized vector.
         * @param constant_value The value of the constant elements.
         */
        void vector_init_constant_elements(vector_omp &obj_const_elem, const int vector_dimension, const double constant_value);

        /**
         * @brief Takes the dimension as an argument to initialize our vector instance with the range
         *        of the given dimension as the elements of our vector.
         *
         * @param obj_range_elem Instance of the struct to initialize vector of a given dimension in place.
         * @param vector_dimension The required dimension of the initialized vector.
         */
        void vector_init_range_elements(vector_omp &obj_range_elem, const int vector_dimension);

        /**
         * @brief Takes the dimension, an instance of our vector stuct and a std vector of doubles with the same
         *        dimension having some values in their elements as arguments. The function copies the values 
         *        element-wise into our vector instance.
         * 
         * @param obj_std_doubles Instance of the struct to hold the values of given std vector of doubles in place.
         * @param std_vector The std vector od doubles whose values have to be copied to our vector.
         * @param vector_dimension The required dimension of the initialized vector.
         */
        void vector_init_std_doubles(vector_omp &obj_std_doubles, std::vector<double> &std_vector, const int vector_dimension);

        /**
         * @brief Takes one instance of the struct as a constant and copies its values element-wise to another instance
         *        of the stuct.
         * 
         * @param obj_initial Instance of the struct to hold the values of given constant struct instance in place.
         * @param obj_to_be_copied The constant struct whose values are to be copied.
         */
        void vector_copy(const vector_omp &obj_initial, vector_omp &obj_to_be_copied);

        /**
         * @brief Takes an instance of the struct and a double "scaling_factor" as arguments and scales the elments
         *        of our vector by a magnitude of scaling factor.
         * 
         * @param obj_scale Instance of the struct to hold the values of scaled elements in place.
         * @param scaling_factor The factor to which the elements are to be scaled.
         */
        void vector_scale(vector_omp &obj_scale, const double scaling_factor);

        /**
         * @brief Computes the dot product of 2 vectors contained in the 2 instances of the struct given as arguments.
         *        Vector x and y are constant throughout the computation.
         * 
         * @param obj_dot_x Arbitrary instance of the struct.
         * @param obj_dot_y Arbitrary instance of the struct.
         * @return double Returns the computed dot product as a double which can be used for further computations.
         */
        double vector_dot_product(const vector_omp &obj_dot_x, const vector_omp &obj_dot_y);

        /**
         * @brief Computes the eucledian norm by just taking the squre root of the dot product from the previous function.
         *        Vector x is constant throughout the computation.
         * 
         * @param obj_norm_x Instance of the struct containing vector required for computing its norm.
         * @return double Returns the computed Eucledian norm as a double which can be used for further computations.
         */
        double vector_euclidean_norm(const vector_omp &obj_norm_x);

        /**
         * @brief Computes scaled vector addition(y := y + αx) of given vectors x and y. Vector x is a constant
         *        throughout the computation.
         * 
         * @param obj_scaled_add_y Instance of the struct containig vector y.
         * @param obj_scaled_add_x Instance of the struct containig vector x.
         * @param scaling_factor_add Scaling factor α.
         */
        void vector_scaled_addition(vector_omp &obj_scaled_add_y, const vector_omp &obj_scaled_add_x, const double scaling_factor_add);

        /**
         * @brief Allocates the memory for the vector created using unique pointer.
         * 
         * @param unique Arbitrary instance of the struct.
         * @param dim The required dimension of the initialized vector.
         * 
         */
        void allocate(vector_omp &unique, const unsigned int dim);



}//end namespace pnla
 
#endif // __VECTOR_OMP_H__