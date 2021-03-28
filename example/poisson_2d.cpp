#include <vector>
#include <iostream>
#include "FD_linear_system.h"
#include "vector_seq.h"
#include "CRS_Matrix.h"
#include "PCG_Solver.h"
#include <limits> //epsilon()
#include <cmath>

/**
 * @brief Checks the convergence of the PCG towards the actual X.
 * 
 * @tparam Matrix Matrix Template argument which can be instantiated with classes having the same charecteristics of struct "CRS_Matrix".
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq".
 * @param x_PCG_result The resultant X obtained after using PCG_Solver.
 * @param b_RHS Given RHS load vector.
 * @param A Given stifness matrix A in the CRS format.
 * @param epsilon tolerable error in computation.
 * @return double Returns the norm of the residue. (b - A * x)
 * 
 */
template<typename Matrix, typename Vector>
double PCG_convergence(Vector x_PCG_result, Vector b_RHS, Matrix A, const double epsilon)
{
    // Create a new vector b_pcg and initialize with b_RHS.
    Vector b_pcg;
    pnla::vector_copy(b_RHS, b_pcg);

    // Computes A * X and stores in b_pcg. X is the result obtained from the PCG_Solver.
    pnla::CRS_scaled_matrix_vector_multiplication(A, x_PCG_result, b_pcg, 1.0, 0.0);

    // Computes the residue and compare with the given b_RHS.
    pnla::vector_scaled_addition(b_pcg, b_RHS, -1.0);
    double norm_b_pcg = pnla::vector_euclidean_norm(b_pcg);

    // Check for absolute accuracy.
    double norm_b_RHS = pnla::vector_euclidean_norm(b_RHS);
    norm_b_pcg = norm_b_pcg / norm_b_RHS;
    
    return norm_b_pcg;
}

/**
 * @brief Application which solves the poisson equation on the unit square.
 * 
 * @tparam Matrix Matrix Template argument which can be instantiated with classes having the same charecteristics of struct "CRS_Matrix".
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq".
 * @param test_sucess_count initialized with 0.
 * @param inner_points Given as an argument to determine the FD_Linear system A, X and B.
 * @param epsilon tolerable error in computation.
 * @return int Returns either 0 or 1.
 * 
 */
template<typename Matrix, typename Vector>
int poisson_2d(int test_sucess_count, const int inner_points, const double epsilon)
{
    // Computes the step size required to initialize CRS_Matrix obtained from FD_linear_system.
    const double h = 1/static_cast<double>(inner_points + 1);
    FD_Linear_System obj_FD_LS(inner_points, h);

    // gets the dofs required to fix the max_iterations for the PCG solver.
    const unsigned int dimension = obj_FD_LS.get_dofs();

    // To get CRS Matrix from the FD_linear_system.
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> rows;
    obj_FD_LS.get_crs_matrix_vectors(values, columns, rows);
    const unsigned int num_of_rows = (rows.size());
    const unsigned int num_non_zero = (values.size());
    
    // To get the vector X and the RHS vector B from the FD_linear_system.
    std::vector<double> x, b;
    obj_FD_LS.get_x(x);
    obj_FD_LS.get_b(b);
    const int x_size = x.size();
    const int b_size = b.size();
    
    // Stores the Vector X in the sequential vector X_FD.
    Vector x_FD;
    pnla::vector_init_std_doubles(x_FD, x, x_size);
    // Stores the Vector B in the sequential vector B_FD.
    Vector b_FD;
    pnla::vector_init_std_doubles(b_FD, b, b_size);

    // Inititalizes the CRS_matrix.
    Matrix A_FD;
    pnla::CRS_Matrix_initialization(A_FD, num_of_rows, num_non_zero, values, columns, rows);

    // PCG solver called.
    double relative_accuracy = 1e-10;
    int iterations = pnla::PCG_Result<pnla::CRS_Matrix, pnla::vector_seq>(A_FD, b_FD, x_FD, relative_accuracy, dimension);

    // Check for convergence and the number of iterations taken to converge.
    double norm_b = PCG_convergence<pnla::CRS_Matrix, pnla::vector_seq>(x_FD, b_FD, A_FD, epsilon);
    
    std::cout << "The norm of b is :    " << norm_b << std::endl;
    std::cout << "Number of iterations taken to converge is :    " << iterations << std::endl;

    if(abs(norm_b) > epsilon)
    {
        std::cout << "PCG did not converge as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;

}

/**
 * @brief Test Programm for pnla's PCG_Solver.
 * 
 * @param argc Number of arguments.
 * @param argv Programm argument list.
 * @return int On sucess = 0, if return value != 0 some test failed.
 * 
 */
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

