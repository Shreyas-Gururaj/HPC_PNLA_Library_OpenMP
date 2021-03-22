#include "CRS_Matrix.h"
#include <cmath>

namespace pnla{

    //
    void CRS_scaled_matrix_vector_multiplication(CRS_Matrix_Representation CRS_Matrix, std::vector<double> y)
    {
        pnla::CRS_Matrix_Representation CRS_1(10, 10, 20);
        
        std::vector<double> A_v = CRS_1.get_non_zero_vector();
        std::vector<double> J_a = CRS_1.col_ind_vector();
        std::vector<double> I_a = CRS_1.row_ind_vector();

    }

    

}// end namespace pnla