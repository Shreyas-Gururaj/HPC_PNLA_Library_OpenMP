#include "PCG_Solver.h"
#include "CRS_Matrix.h"
#include "vector_seq.h"

namespace pnla{


template<typename Matrix, typename Vector>
void PCG_Result(const pnla::CRS_Matrix &CRS_Matrix_seq_A, const pnla::vector_seq &b_RHS_Vector, pnla::vector_seq &x_PCG_result, const double rel_accuracy,
                  const unsigned int max_iterations)
{

    Vector residue_vector;
    pnla::vector_copy(b_RHS_Vector, residue_vector);               // r_0 = b
    pnla::CRS_scaled_matrix_vector_multiplication(CRS_Matrix_seq_A, x_PCG_result, residue_vector, -1.0, 1.0);   // r_0 = r_0 - A * x.

    Vector vector_P , vector_V, vector_W;     //As required in the algorithm.

    pnla::vector_copy(residue_vector, vector_P);   //Assuming conditioner is an identity matrix.
    pnla::vector_copy(residue_vector, vector_V);    //Assuming conditioner is an identity matrix.
    pnla::vector_copy(residue_vector, vector_W);

    double rho = pnla::vector_dot_product(vector_V, residue_vector);
    double gamma = rho;     //required for further computations as rho(k) cannot be stored on the run.
    double alpha = 0.0;
    unsigned int k = 0;     //For the while loop

    double norm_b = pnla::vector_euclidean_norm(b_RHS_Vector);
    const double relative_accuracy = relative_accuracy * norm_b;

    do
    {
        pnla::CRS_scaled_matrix_vector_multiplication(CRS_Matrix_seq_A, vector_P, vector_W, 1.0, 0.0);  // W(k) = A * P(k)

        alpha = rho / pnla::vector_dot_product(vector_P, vector_W); //As defined

        pnla::vector_scaled_addition(x_PCG_result, vector_P, alpha); // X(k+1) = X(k) + alpha * P(K)
        pnla::vector_scaled_addition(residue_vector, vector_W, -alpha); // R(k+1) = R(k) - alpha * W(K)
        pnla::vector_copy(residue_vector, vector_V);                    // V(k+1) = C * R(k+1), Assuming C as identity.
        pnla::vector_dot_product(vector_V, residue_vector);

        gamma = rho / gamma;            // computes gamma(k) = rho(k+1) / rho(k)

        pnla::vector_scaled_addition(vector_P, vector_V, (1.0 / gamma));
        pnla::vector_scale(vector_P, gamma);
        k = k + 1;

    } while ((pnla::vector_euclidean_norm(residue_vector) > relative_accuracy) && (k <= max_iterations));
}

}// end namespace pnla