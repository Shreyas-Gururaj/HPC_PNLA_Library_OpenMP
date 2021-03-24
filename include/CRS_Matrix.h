
#ifndef __CRS_MATRIX_H__
#define __CRS_MATRIX_H__

#include <iostream>
#include <vector>
#include <string>


namespace pnla{

    struct CRS_Matrix
    {
            unsigned int total_num_of_rows;
            unsigned int total_non_zero_elements;

            std::vector<double> Matrix_non_zero_elements;
            std::vector<double> Col_indices_non_zero_elements;
            std::vector<double> Row_indices_non_zero_elements;
    };


            //
            void CRS_Matrix_initialization(CRS_Matrix &obj_CRS_init, const unsigned int num_of_rows,
                                           const unsigned int num_non_zero_entries);

            //
            void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix &obj_CRS_Matrix, const std::vector<double> x, 
                                                         std::vector<double> y, const double alpha, const double beta);

}//end namespace pnla

#endif // __CRS_MATRIX_H__