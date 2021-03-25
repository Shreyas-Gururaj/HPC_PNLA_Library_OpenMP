
#ifndef __CRS_MATRIX_H__
#define __CRS_MATRIX_H__

#include <vector>
#include "vector_seq.h"


namespace pnla{

    struct CRS_Matrix
    {
        unsigned int total_num_of_rows;
        unsigned int total_non_zero_elements;

        std::vector<double> Matrix_non_zero_elements;
        std::vector<int> Col_indices_non_zero_elements;
        std::vector<int> Row_indices_non_zero_elements;
    };


            void CRS_Matrix_initialization(CRS_Matrix &obj_CRS_init, const unsigned int num_of_rows, const unsigned int num_non_zero_entries, 
                                           const std::vector<double> &values, const std::vector<int> &col_index, 
                                           const std::vector<int> &row_index);

            //
            void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix &obj_CRS_Matrix, const vector_seq &x,
                                                         vector_seq &y, const double alpha, const double beta);

}//end namespace pnla

#endif // __CRS_MATRIX_H__