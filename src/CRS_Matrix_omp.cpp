/**
 * @file CRS_Matrix.cpp
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief Implementation of all the routines for CRS_Matrix_omp.
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "CRS_Matrix_omp.h"
#include <cmath>
#include <omp.h>


namespace pnla{

     void CRS_Matrix_initialization(CRS_Matrix_omp &obj_CRS_init, const unsigned int num_of_rows, const unsigned int num_non_zero_entries, 
                                    const std::vector<double> &values, const std::vector<int> &col_index, const std::vector<int> &row_index)
    {
        allocate_non_zero(obj_CRS_init, num_non_zero_entries);
        allocate_columns(obj_CRS_init, num_non_zero_entries);
        allocate_rows(obj_CRS_init, num_of_rows + 1);

        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < num_non_zero_entries; elem_index++)
        {
            obj_CRS_init.Matrix_non_zero_elements[elem_index] = values[elem_index];
        }

        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < num_non_zero_entries; elem_index++)
        {
            obj_CRS_init.Col_indices_non_zero_elements[elem_index] = col_index[elem_index];
        }

        #pragma omp parallel for schedule(static)
        for(unsigned int elem_index = 0; elem_index < num_of_rows; elem_index++)
        {
            obj_CRS_init.Row_indices_non_zero_elements[elem_index] = row_index[elem_index];
        }

    }

    void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix_omp &obj_CRS_Matrix, const vector_omp &x, 
                                                 vector_omp &y, const double alpha, const double beta)
    {

        allocate(y, static_cast<int>(y.vector_dimension));

        #pragma omp parallel for schedule(static)
        for(unsigned int i = 0; i < obj_CRS_Matrix.total_num_of_rows - 1; i++)
        {
            y.values[i] *= beta;

            for(int k = obj_CRS_Matrix.Row_indices_non_zero_elements[i]; k < obj_CRS_Matrix.Row_indices_non_zero_elements[i + 1]; k++)
                { 
                    y.values[i] += alpha * obj_CRS_Matrix.Matrix_non_zero_elements[k] * x.values[obj_CRS_Matrix.Col_indices_non_zero_elements[k]];
                }
        }

    }

    void allocate_non_zero(CRS_Matrix_omp &unique, const unsigned int non_zero)
   {
       unique.total_non_zero_elements = non_zero;
       unique.Matrix_non_zero_elements.reset(new double[non_zero]);
   }


    void allocate_columns(CRS_Matrix_omp &unique, const unsigned int non_zero)
   {
       unique.total_non_zero_elements = non_zero;
       unique.Col_indices_non_zero_elements.reset(new int[non_zero]);
   }

    void allocate_rows(CRS_Matrix_omp &unique, const unsigned int rows)
   {
       unique.total_num_of_rows = rows;
       unique.Row_indices_non_zero_elements.reset(new int[rows]);
   }
   
}// end namespace pnla