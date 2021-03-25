#ifndef __PCG_Solver_H__
#define __PCG_Solver_H__

#include <vector>
#include "vector_seq.h"
#include "CRS_Matrix.h"



namespace pnla{


//
template<typename Matrix, typename Vector>
void PCG_Result(const pnla::CRS_Matrix &CRS_Matrix_seq, const pnla::vector_seq &b, pnla::vector_seq &x_initial, const double rel_accuracy,
                  const unsigned int max_iterations);


}//end namespace pnla

#endif // __PCG_Solver_H__