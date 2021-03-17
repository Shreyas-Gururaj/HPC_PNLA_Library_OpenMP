/**
 * @file my_utilities.h
 * @author Thomas Rau (thomas1.rau@uni-bayreuth.de)
 * @brief Shows usage of function templates. Provided function templates may be useful 
 * @version 0.1
 * @date 2021-03-14
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef __MY_UTILITIES_H__
#define __MY_UTILITIES_H__

#include <iostream>
#include <string>


/**
 * @brief Prints the content of array on console in a structured way.
 * 
 * @tparam array_type Some type which behaves like an array meaning it provides the []
 *         operator of elementwise access. Implemented for:
 *         - std::vector<double>
 * @param array Array holding some data to be printed
 * @param array_lenght length of array.
 * @param name Name of array to be printed on console
 * 
 * Since definition of fuction template is in cpp file, 
 */
template<typename array_type>
void print_array(const array_type &array, const int array_lenght, const std::string &name);


/**
 * @brief Prints vector on console in structured fashion.
 * 
 * @tparam Vector Class/structure with public member values and dimension.
 *         values is provides array like reading access to elements via [] operator
 *         dimension describes the number of elements stored in values.
 * @param v some vector,
 * @param name name to be printed on console
 * 
 * As long as you implement to the interface described in tparam you can use this 
 * function to print your vector version on the console. You may find them useful.
 */
template<typename Vector>
void print_vector(const Vector &v, const std::string &name)
{
    print_array(v.values, v.dimension, name);
}



#endif // __MY_UTILITIES_H__