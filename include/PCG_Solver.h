#ifndef __PCG_Solver_H__
#define __PCG_Solver_H__

#include <vector>
#include "vector_seq.h"
#include "CRS_Matrix.h"



namespace pnla{


//
template<typename Matrix, typename Vector>
int PCG_Result(const Matrix &CRS_Matrix_A, const Vector &b_RHS_Vector, Vector &x_PCG_result, const double rel_accuracy,
                  const unsigned int max_iterations);


}//end namespace pnla

#endif // __PCG_Solver_H__