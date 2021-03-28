
# To add an executable named dev made out of develop.cpp and my_utilities.cpp
add_executable(dev develop.cpp my_utilities.cpp)

# Tells where headers of target are to be found
target_include_directories(dev PUBLIC {CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(dev PUBLIC ../include)

# Links first_app to library pnla
target_link_libraries(dev pnla)



##### For poisson_2d

# Tests need to be added as executables first
    
add_executable(test_poisson_2d poisson_2d.cpp FD_linear_system.cpp)

# Tells where headers of target are to be found
#target_include_directories(test_poisson_2d {CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(test_poisson_2d PUBLIC ../include)

# Links test_poisson_2d against library pnla
target_link_libraries(test_poisson_2d pnla)

### in what follows, we create a automated testing enviroment with cmake 

enable_testing()

# Add a testing call of test_poisson_2d with no arguments
add_test(NAME test_poisson_2d_default COMMAND test_poisson_2d)

# Add a testing call of test_poisson_2d with one argument (inner points = 100)
add_test(NAME test_poisson_2d_inner_100 COMMAND test_poisson_2d 100)

# # Add a testing call of test_poisson_2d with one argument (inner points = 200)
add_test(NAME test_poisson_2d_inner_200 COMMAND test_poisson_2d 200)

# # Add a testing call of test_poisson_2d with one argument (inner points = 400)
add_test(NAME test_poisson_2d_inner_400 COMMAND test_poisson_2d 400)

# # Add a testing call of test_poisson_2d with one argument (inner points = 800)
add_test(NAME test_poisson_2d_inner_800 COMMAND test_poisson_2d 800)

# # Add a testing call of test_poisson_2d with one argument (inner points = 1600)
add_test(NAME test_poisson_2d_inner_1600 COMMAND test_poisson_2d 1600)
