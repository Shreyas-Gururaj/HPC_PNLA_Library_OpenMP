/**
 * @file CRS_Matrix.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-03-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include "CRS_Matrix.h"
#include <cmath>

namespace pnla{

    //
     void CRS_Matrix_initialization(CRS_Matrix &obj_CRS_init, const unsigned int num_of_rows, const unsigned int num_non_zero_entries, 
                                    const std::vector<double> &values, const std::vector<int> &col_index, const std::vector<int> &row_index)
    {
        obj_CRS_init.total_non_zero_elements = num_non_zero_entries;
        obj_CRS_init.total_num_of_rows = num_of_rows;

        obj_CRS_init.Matrix_non_zero_elements.resize(num_non_zero_entries);
        obj_CRS_init.Col_indices_non_zero_elements.resize(num_non_zero_entries);
        obj_CRS_init.Row_indices_non_zero_elements.resize(num_of_rows);

        for(unsigned int elem_index = 0; elem_index < num_non_zero_entries; elem_index++)
        {
            obj_CRS_init.Matrix_non_zero_elements[elem_index] = values[elem_index];
        }

        for(unsigned int elem_index = 0; elem_index < num_non_zero_entries; elem_index++)
        {
            obj_CRS_init.Col_indices_non_zero_elements[elem_index] = col_index[elem_index];
        }

        for(unsigned int elem_index = 0; elem_index < num_of_rows; elem_index++)
        {
            obj_CRS_init.Row_indices_non_zero_elements[elem_index] = row_index[elem_index];
        }

    }

    //
    void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix &obj_CRS_Matrix, const vector_seq &x, 
                                                 vector_seq &y, const double alpha, const double beta)
    {
        for(unsigned int i = 0; i < obj_CRS_Matrix.Row_indices_non_zero_elements.size() - 1; i++)
        {
            y.values[i] *= beta;
            
            for(int k = obj_CRS_Matrix.Row_indices_non_zero_elements[i]; k < obj_CRS_Matrix.Row_indices_non_zero_elements[i + 1]; k++)
                { 
                    y.values[i] += alpha * obj_CRS_Matrix.Matrix_non_zero_elements[k] * x.values[obj_CRS_Matrix.Col_indices_non_zero_elements[k]];
                }
            
        }
       
    }
}// end namespace pnla