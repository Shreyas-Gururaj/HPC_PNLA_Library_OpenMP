/**
 * @file my_utilites.cpp
 * @author Thomas Rau (thomas1.rau@uni-bayreuth.de)
 * @brief Shows usage of function templates. Provided function templates may be useful
 * @version 0.1
 * @date 2021-03-14
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "my_utilities.h"
#include <vector>


/**
 * Define fuctiontemplate. It means you provide a recipe to the compiler on how to 
 * build a fuction. 
 */
template<typename array_type>
void print_array(const array_type &array, const int array_lenght, const std::string &name) 
{
    std::cout << name << ":" <<std::endl; 

    for(int i=0; i<array_lenght;i++)
        std::cout<< array[i]<<std::endl;
    
    std::cout <<std::endl<<std::endl;  
}

/**
 * This tells the compiler to actually bild a function print_array using a function 
 * template. 
 * Here print_array is buildt for std::vector<double> 
 * This is needed if you put definition of function template into cpp file. In that way
 * you trade flexibility against security regarding the user because now the function 
 * template only can be used with actually pre build versions
 */
template void print_array(const std::vector<double> &v, const int array_lenght, 
                          const std::string &name); 