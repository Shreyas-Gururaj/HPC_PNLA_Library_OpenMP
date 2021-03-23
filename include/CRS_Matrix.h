
#ifndef __CRS_MATRIX_H__
#define __CRS_MATRIX_H__

#include <iostream>
#include <vector>
#include <string>


namespace pnla{

    class CRS_Matrix
    {
        public:

            unsigned int total_num_of_rows;
            unsigned int total_num_of_columns;
            unsigned int total_non_zero_elements;

            std::vector<double> Matrix_non_zero_elements;
            std::vector<double> Col_indices_non_zero_elements;
            std::vector<double> Row_indices_non_zero_elements;

            CRS_Matrix() = delete;


            //
            CRS_Matrix(const unsigned int num_of_rows, const unsigned int num_of_columns, const unsigned int num_non_zero_entries)
            {
                total_num_of_rows = num_of_rows;
                total_num_of_columns = num_of_columns;
                total_non_zero_elements = num_non_zero_entries;

                Matrix_non_zero_elements.assign(num_non_zero_entries, 0.0);
                Col_indices_non_zero_elements.assign(num_non_zero_entries, 0.0);
                Row_indices_non_zero_elements.assign(num_of_rows + 1, 0.0);
            }

            //
            void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix &obj_CRS_Matrix, const std::vector<double> x, 
                                                         std::vector<double> y, const double alpha, const double beta);
    };
}//end namespace pnla

#endif // __CRS_MATRIX_H__