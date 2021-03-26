#include <vector>
#include <iostream>
#include "FD_linear_system.h"
#include "vector_seq.h"
#include "CRS_Matrix.h"
#include "PCG_Solver.h"
#include <limits> //epsilon()
#include <cmath>

//
template<typename Matrix, typename Vector>
double PCG_convergence(Vector b_pcg, Vector x_PCG_result, Vector b_RHS, Matrix A, const double epsilon)
{
    pnla::vector_copy(b_RHS, b_pcg);
    pnla::CRS_scaled_matrix_vector_multiplication(A, x_PCG_result, b_pcg, 1.0, 0.0);
    pnla::vector_scaled_addition(b_pcg, b_RHS, -1.0);
    double norm_b_pcg = pnla::vector_euclidean_norm(b_pcg);
    
    return norm_b_pcg;
}

//
template<typename Matrix, typename Vector>
int poisson_2d(int test_sucess_count, const int inner_points, const double epsilon)
{
    //
    const double h = 1/static_cast<double>(inner_points + 1);
    FD_Linear_System obj_FD_LS(inner_points, h);

    //
    const unsigned int dimension = obj_FD_LS.get_dofs();
    double relative_accuracy = 1e-14;

    //
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> rows;
    obj_FD_LS.get_crs_matrix_vectors(values, columns, rows);
    const unsigned int num_of_rows = (rows.size());
    const unsigned int num_non_zero = (values.size());
    
    //
    std::vector<double> x, b;
    obj_FD_LS.get_x(x);
    obj_FD_LS.get_b(b);
    const int x_size = x.size();
    const int b_size = b.size();
    
    //
    Vector x_FD;
    pnla::vector_init_std_doubles(x_FD, x, x_size);
    //
    Vector b_FD;
    pnla::vector_init_std_doubles(b_FD, b, b_size);

    //
    Matrix A;
    pnla::CRS_Matrix_initialization(A, num_of_rows, num_non_zero, values, columns, rows);

    //
    pnla::PCG_Result<pnla::CRS_Matrix, pnla::vector_seq>(A, b_FD, x_FD, relative_accuracy, dimension);

    //
    Vector b_pcg;
    double norm_b = PCG_convergence<pnla::CRS_Matrix, pnla::vector_seq>(b_pcg, x_FD, b_FD, A, epsilon);
    //std::cout << "The norm of b is :  " << norm_b << std::endl;

    if(abs(norm_b) > epsilon)
    {
        std::cout << "PCG did not converge as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;

}


int main(int argc, char *argv[])
{
    int inner_points = 20;

    int test_result = 0;
    const double epsilon(std::numeric_limits<double>::epsilon());

    if(argc == 2)
	{
         inner_points = std::stoi(argv[1]);
	}

    test_result += poisson_2d<pnla::CRS_Matrix, pnla::vector_seq>(test_result, inner_points, epsilon);

    if(test_result !=0 )
        return test_result;

    return test_result;
    return 0;
}

