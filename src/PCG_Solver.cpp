/**
 * @file PCG_Solver.cpp
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief PCG solver implementation.
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "PCG_Solver.h"
#include <iostream>

namespace pnla{


template<typename Matrix, typename Vector>
int PCG_Result(const Matrix &CRS_Matrix_A, const Vector &b_RHS_Vector, Vector &x_PCG_result, const double rel_accuracy,
                  const unsigned int max_iterations)
{

    Vector residue_vector;
    //b_RHS_Vector.vector_dimension = b_RHS_Vector.values.size();
    vector_copy(b_RHS_Vector, residue_vector);               // initiates the residue_vector, r_0 with the RHS b_RHS_Vector.
    CRS_scaled_matrix_vector_multiplication(CRS_Matrix_A, x_PCG_result, residue_vector, -1.0, 1.0);   //


    Vector vector_P , vector_V, vector_W;     //As required in the algorithm.

    vector_copy(residue_vector, vector_P);   //Assuming conditioner is an identity matrix.
    vector_copy(residue_vector, vector_V);    //Assuming conditioner is an identity matrix.
    vector_init_constant_elements(vector_W, b_RHS_Vector.vector_dimension, 0.0);

    double rho = vector_dot_product(vector_V, residue_vector);
    double gamma = 0.0;     //required for further computations as rho(k) cannot be stored on the run otherwise.
    double alpha = 0.0;
    unsigned int k = 0;     //For the while loop.

    double norm_b = vector_euclidean_norm(b_RHS_Vector);
    const double relative_accuracy = rel_accuracy * norm_b;

    do
    {
        CRS_scaled_matrix_vector_multiplication(CRS_Matrix_A, vector_P, vector_W, 1.0, 0.0);  // computes W(k) = A * P(k).

        alpha = rho / vector_dot_product(vector_P, vector_W); // computes α(k) := ρ(k)/(P(k), W(k)).

        vector_scaled_addition(x_PCG_result, vector_P, alpha); // computes X(k+1) = X(k) + alpha * P(K).
        vector_scaled_addition(residue_vector, vector_W, -alpha); // computes R(k+1) = R(k) - alpha * W(K).
        vector_copy(residue_vector, vector_V);                    // computes V(k+1) = C * R(k+1), Assuming C as identity.

        
        gamma = rho;
        rho = vector_dot_product(vector_V, residue_vector);       //computes ρ(k+1) = (V(k+1), R(k+1)).
        gamma = rho / gamma;                                     // computes γk = ρ(k+1) / ρ(k).

        vector_scale(vector_P, gamma);  // Pre-computation for computing P(k+1) = V(k+1) + γ(k) * P(k).
        vector_scaled_addition(vector_P, vector_V, 1.0); // computes P(k+1) = V(k+1) + γ(k) * P(k).
    
        k = k + 1;

    } while ((vector_euclidean_norm(residue_vector) > relative_accuracy) && (k <= max_iterations));
      // While loop exits if either required relative accuracy is achieved or the number of iterations reaches max_iterations.
    //std::cout << vector_euclidean_norm(residue_vector) << std::endl;
    return k;  // Returns the number of iterations required to reach the specified relative accuracy.
}

template int PCG_Result(const CRS_Matrix &CRS_Matrix_A, const vector_seq &b_RHS_Vector, vector_seq &x_PCG_result, const double rel_accuracy,
                  const unsigned int max_iterations);

template int PCG_Result(const CRS_Matrix_omp &CRS_Matrix_A, const vector_omp &b_RHS_Vector, vector_omp &x_PCG_result, const double rel_accuracy,
                  const unsigned int max_iterations);
}// end namespace pnla