/**
 * @file CRS_Matrix.h
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief Header files for sequential CRS_Matrix functions requird for the PCG algorithm.
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CRS_MATRIX_H__
#define __CRS_MATRIX_H__

#include <vector>
#include "vector_seq.h"

namespace pnla{

    /**
     * @brief Struct "CRS_Matrix contains 3 vectors to store the given matrix in the compressed row storage format
     *        required to initialize and compute matrix vector multiplication for reformulated for CRS format.
     * 
     */
    struct CRS_Matrix
    {
        // Stores the number of rows of the given matrix A.
        unsigned int total_num_of_rows;
        // Stores the number of non zero elements in the matrix A.
        unsigned int total_non_zero_elements;

        // Stores the values of the non zero entries of the matrix A. (Av)
        std::vector<double> Matrix_non_zero_elements;
        // Stores the column indices of the respective non zero entries. (Ja)
        std::vector<int> Col_indices_non_zero_elements;
        // Stores the row indices of the respective non zero entries. (Ia)
        std::vector<int> Row_indices_non_zero_elements;
    };

            /**
             * @brief Stores the values of a stiffness matrix A obtained from FD or FE discretization (Sparse, Positive definite 
             *        and symmetric) in the CRS format containing 3 vectors i.e, 
             *        1. Vector to store all the non zero entries (Av).
             *        2. Vector to store the column indices of the respective non zero entries (Ja).
             *        3. Vector to store the row indices of the respective non zero entries (Ia).
             * 
             * @param obj_CRS_init Instance of the struct to hold the 3 vectors in the CRS format.
             * @param num_of_rows The number of rows of the input matrix A.
             * @param num_non_zero_entries The number of non zero elements in the Matrix A.
             * @param values This vector stores the values of all the non zero entries row-wise. (std vector of doubles)
             * @param col_index This vector stores the column indices of all the non zero entries respectively. (std vector of int)
             * @param row_index This vector stores the row indices of all the non zero entries respectively. (std vector of int)
             */
            void CRS_Matrix_initialization(CRS_Matrix &obj_CRS_init, const unsigned int num_of_rows, 
                                           const unsigned int num_non_zero_entries, const std::vector<double> &values, 
                                           const std::vector<int> &col_index, const std::vector<int> &row_index);

            /**
             * @brief Computes the scaled matrix vector multiplication (y := ??Ax + ??y). Vector x remains constant
             *        throughout the computation. ?? and ?? are constant doubles given as arguments.
             * 
             * @param obj_CRS_Matrix Instance of the struct to compute the scaled matrix vector multiplication.
             * @param x arbitrary instance of a "vector_seq" struct used for compuation.
             * @param y arbitrary instance of a "vector_seq" struct used for compuation.
             * @param alpha Scaling factor for AX.
             * @param beta Scaling factor for y.
             */
            void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix &obj_CRS_Matrix, const vector_seq &x,
                                                         vector_seq &y, const double alpha, const double beta);

}//end namespace pnla

#endif // __CRS_MATRIX_H__