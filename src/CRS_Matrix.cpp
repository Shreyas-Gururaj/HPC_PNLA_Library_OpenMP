#include "CRS_Matrix.h"
#include <cmath>

namespace pnla{

    //
    void CRS_Matrix_initialization(CRS_Matrix &obj_CRS_init, const unsigned int num_of_rows, 
                                   const unsigned int num_non_zero_entries)
    {
        obj_CRS_init.total_non_zero_elements = num_non_zero_entries;
        obj_CRS_init.total_num_of_rows = num_of_rows;

        obj_CRS_init.Matrix_non_zero_elements.assign(obj_CRS_init.total_non_zero_elements, 0.0);
        obj_CRS_init.Col_indices_non_zero_elements.assign(obj_CRS_init.total_non_zero_elements, 0.0);
        obj_CRS_init.Row_indices_non_zero_elements.assign(obj_CRS_init.total_num_of_rows, 0.0);
    }

    //
    void CRS_scaled_matrix_vector_multiplication(const CRS_Matrix &obj_CRS_Matrix, const std::vector<double> x, 
                                                 std::vector<double> y, const double alpha, const double beta)
    {
        for(unsigned int i = 0; i < obj_CRS_Matrix.total_num_of_rows; i++)
        {
            y[i] = 0.0;
            
            for(unsigned int k = obj_CRS_Matrix.Row_indices_non_zero_elements[i]; k < obj_CRS_Matrix.Row_indices_non_zero_elements[i + 1] - 1; k++)
                {
                    y[i] = (alpha * (obj_CRS_Matrix.Matrix_non_zero_elements[k] * x[obj_CRS_Matrix.Row_indices_non_zero_elements[k]])) + (beta * y[i]);
                }
        }

       
    }

    

}// end namespace pnla