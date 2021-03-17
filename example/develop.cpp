/**
 * @file develop.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-03-04
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <vector>
#include "my_utilities.h"
#include <vector_seq.h>

/**
 * @brief Small programm to show use of pnla, cmake and templates
 * 
 * @return int 
 */
int main()
{
    std::cout<<std::endl;
    std::cout<<"Start Test"<<std::endl<<std::endl;

    /// 
    const int dimension = 10;

    std::vector<double> test_vector(dimension, 1.0);
    print_array(test_vector, test_vector.size(),"test_vector");

    /// Shows the use of a vector struct/class of your PNLA-implementation
    //pnla::your_vector x;

    std::cout<< "The End " << std::endl<<std::endl;
    
    return 0;
}
