/**
 * @file test_CRS_Matrix.cpp
 * @author Shreyas Gururaj (Shreyas.Gururaj@uni-bayreuth.de)
 * @brief 
 * @version 0.1
 * @date 2021-03-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <string>
#include <limits> //epsilon()
#include "CRS_Matrix.h"
#include "FD_linear_system.h"
#include <cmath>

/**
 * @brief Checks for the correctness of the functions "CRS_Matrix_initialization" and 
 *        "CRS_scaled_matrix_vector_multiplication" simultaneously.
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
int test_Matrix_init(int test_sucess_count, const int inner_points, const double epsilon)
{
    // Computes the step size required to initialize CRS_Matrix obtained from FD_linear_system.
    const double h = 1/static_cast<double>(inner_points + 1);
    FD_Linear_System obj_FD_LS(inner_points, h);

    // To get CRS Matrix from the FD_linear_system.
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> rows;
    obj_FD_LS.get_crs_matrix_vectors(values, columns, rows);

    // To get the test vector X containing the row sum of A * X(1).
    std::vector<double> test_x;
    obj_FD_LS.get_test_vector(test_x);
    const int test_x_size = test_x.size();

    // Storing the values obtained from the FD_linear system in the CRS_Matrix.
    const unsigned int num_of_rows = (rows.size());
    const unsigned int num_non_zero = (values.size());
    Matrix CRS_matrix_A;
    pnla::CRS_Matrix_initialization(CRS_matrix_A, num_of_rows, num_non_zero, values, columns, rows);

    // To create two vectors and carry out scaled matrix vector multiplication. Vector X has all 1's and Vector Y is a 0 vector.
    Vector x;
    Vector y;
    pnla::vector_init_constant_elements(x, num_of_rows, 1.0);
    pnla::vector_init_constant_elements(y, num_of_rows, 0.0);
    pnla::CRS_scaled_matrix_vector_multiplication(CRS_matrix_A, x, y, 1.0,  1.0);

    // Comparing the result with the test vector x.
    Vector text_x_seq;
    pnla::vector_init_std_doubles(text_x_seq, test_x, test_x_size);
    pnla::vector_scaled_addition(y, text_x_seq, -1.0);
    const double norm_y = pnla::vector_euclidean_norm(y);

    std::cout << "the norm is: " << norm_y << std::endl;

    if(abs(norm_y) > epsilon)
    {
        std::cout << "Either the std double initialization or the CRS initilization or the scaled matrix vector addition does not work as intended";
        return test_sucess_count + 1;
    }

    return test_sucess_count;
}

/**
 * @brief Appends the value of the test_sucess initialized with 0 by a positive integer if any of the test fails.
 * 
 * @tparam Matrix Matrix Template argument which can be instantiated with classes having the same charecteristics of struct "CRS_Matrix".
 * @tparam Vector Vector Template argument which can be instantiated with classes having the same charecteristics of struct "vector_seq".
 * @param inner_points Given as an argument to determine the FD_Linear system A, X and B.
 * @param epsilon tolerable error in computation.
 * @return int int Returns either 0 or a positive integer.
 */
template<typename Matrix, typename Vector>     //in the template function definition "Vector" is an alias for the struct
int test_vector_routines(const int inner_points, const double epsilon)
{   
    int test_sucess = 0;
 
    test_sucess += test_Matrix_init<Matrix, Vector>(test_sucess, inner_points, epsilon);

    return test_sucess;
}

/**
 * @brief Test Programm for pnla's CRS_Matrix structs/classes.
 * 
 * @param argc Number of arguments.
 * @param argv Programm argument list.
 * @return int On sucess = 0, if return value != 0 some test failed.
 */
int main(int argc, char *argv[])
{
    int total_inner_points = 10;


 	if(argc == 2)
	{
        total_inner_points = std::stoi(argv[1]);
	}

    const double epsilon(std::numeric_limits<double>::epsilon()); 
    int test_result = 0; 

    std::cout<<"Test sequential CRS_Matrix"<<std::endl;
    /// call of test_vector should look something like this
    test_result = test_vector_routines<pnla::CRS_Matrix, pnla::vector_seq>(total_inner_points, epsilon);

    std::cout << test_result << std::endl;
    if(test_result !=0 )
        return test_result;

    return test_result;
}



