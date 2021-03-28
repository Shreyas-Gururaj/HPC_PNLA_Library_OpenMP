/**
 * @file PCG_Solver.h
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief 
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PCG_Solver_H__
#define __PCG_Solver_H__

#include <vector>
#include "vector_seq.h"
#include "CRS_Matrix.h"



namespace pnla{


/**
 * @brief Template function to compute the unknown x in the linear system of equations Ax = b using Pre-conditioned 
 *        conjugate gradient method iteratively from an arbitrary x_0. Same function can be used for parallel and
 *        distributed vector and matrix routines.
 * 
 * @tparam Matrix Template argument which can be instantiated with classes having the same charecteristics of struct "CRS_Matrix"
 * @tparam Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq"
 * @param CRS_Matrix_A Arbitrary CRS_Matrix to hold the values, column and row indices of non zero entries of Matrix A
 * @param b_RHS_Vector Arbitrary Vector to hold the values of the RHS load vector b.
 * @param x_PCG_result Arbitrary Vector containing the initial x_0
 * @param rel_accuracy The desired relative accuracy to indicate when the computation can stop.
 * @param max_iterations The maximum number of steps to indicate when the computation can stop.
 * @return int Returns the number of iterations taken to obtain desired relative accuracy.
 * 
 */
template<typename Matrix, typename Vector>
int PCG_Result(const Matrix &CRS_Matrix_A, const Vector &b_RHS_Vector, Vector &x_PCG_result, const double rel_accuracy,
                  const unsigned int max_iterations);


}//end namespace pnla

#endif // __PCG_Solver_H__