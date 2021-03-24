
# Tests need to be added as executables first
#add_executable(test_vectors test_vectors.cpp ../src/vector_seq.cpp)   
add_executable(test_vectors test_vectors.cpp)  

# Tells where headers of target are to be foundclea
#target_include_directories(test_vectors PUBLIC {CMAKE_CURRENT_SOURCE_DIR} )
target_include_directories(test_vectors PUBLIC ../include)

# Links test_vectors against library pnla
target_link_libraries(test_vectors pnla)

### in what follows, we create a automated testing enviroment with cmake 

enable_testing()

### Test vectors
# Add a testing call of test_vectors with no arguments
add_test(NAME Test_vectors_default COMMAND test_vectors)

# Add a testing call of test_vectors with one argument (dimension = 1000)
add_test(NAME Test_vectors_dim_1000 COMMAND test_vectors 1000)

# Add a testing call of test_vectors with one argument (dimension = 10000)
add_test(NAME Test_vectors_dim_10000 COMMAND test_vectors 10000)

# Add a testing call of test_vectors with one argument (dimension = 10000000)
#add_test(NAME Test_vectors_dim_10000000 COMMAND test_vectors 10000000)


