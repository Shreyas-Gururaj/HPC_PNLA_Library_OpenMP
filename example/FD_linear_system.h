/**
 * @file FD_linear_system.h
 * @author Thomas Rau (thomas1.rau@uni-bayreuth.de)
 * @brief Header file for Finite Differences Linear System Class
 * @version 0.1
 * @date 2021-03-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __FD_LINEAR_SYSTEM_H__
#define __FD_LINEAR_SYSTEM_H__

#include <vector>


/**
 * @brief Class provides matrix resulting from finite difference discretization of
 *        2D Poisson's equation on a square. Resulting systems equals FE discretization
 *        if courants elements are used.
 */
class FD_Linear_System{

public:

    /// deletes default constructor 
    FD_Linerar_System() = delete;

    /**
     * @brief Construct a new fd linear system object
     * 
     * @param N Number of inner points per dimension
     * @param h Stepsize or mesh width
     */
    FD_Linear_System(const int N, const double h):
    m_N(N), m_h(h),pi(3.14159265358979323846){};

    /**
     * @return Number of Degrees of freedom of system 
     */
    inline int get_dofs()const 
    {
        return m_N*m_N;
    };

    /**
     * @brief Get the crs matrix vectors 
     * 
     * @param[in] values arbitrary double vector 
     * @param[out] values holds non-zero entries of FD-System matrix 
     * @param[in] col_indices arbitrary int vector 
     * @param[out] col_indices holds column indices of non-zero elements
     * @param[in] row_indices arbitrary int vector 
     * @param[out] row_indices holds row indices of non-zero elements
     */
    void get_crs_matrix_vectors(std::vector<double> &values,
                                std::vector<int> &col_indices,
                                std::vector<int> &row_indices) const;


    /**
     * @brief Get the test vector. Vector is the row sum of the matrix entries.
     *        Useful for testing your crs-matrix vector multiplication
     * 
     * @param[in] test_vector arbitrary double vector
     * @param[out] test_vector holds the row sum on matrix
     */
    void get_test_vector(std::vector<double> &test_vector)const;

    /**
     * @brief Get the x vector of the FD linear system Ax = b
     * 
     * @param[in] x arbitrary double vector
     * @param[out] x holds the the solution of the system 
     */
    void get_x(std::vector<double> &x)const;

    /**
     * @brief Gets the right hand side of FD system Ax = b
     * 
     * @param[in] b arbitrary double vector
     * @param[out] b holds the right hand side 
     */
    void get_b(std::vector<double> &b)const;

private:

    /// Stores number of inner points per dimension
    const int m_N;

    /// Stores mesh size of discretization
    const double m_h;


    /// stores \f$ \pi \f$
    const double pi;

};

#endif // __FD_LINEAR_SYSTEM_H__