# First Thing, Check minimum  version of CMake needed
cmake_minimum_required(VERSION 3.5)

# Set project name and some stuff 
project(PNLA VERSION 1.0
		LANGUAGES C CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Set's the environment varialble CMAKE_MODULE_PATH where different 
# cmake modules in the format of .cmake files are located
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

# load compiler flags defined in CompilerFlags.cmake
include(CompilerFlags)


# Here your actual project follows
# Our project is the creation of a Parallel Numerical Linear Algebra library called 
# it will be called

# Adds a source file directory named PNLA_SOURCE
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src PNLA_SOURCE)

# Will create a library called pnla from the sources found in PNLA_SOURCE
add_library(pnla ${PNLA_SOURCE} )

# Tell where header files for pnla are located
target_include_directories(pnla PUBLIC include)


# add (optional) subdirectories with corresponding responsibility, each 
# need it's own CMakeLists.txt file!
add_subdirectory(example)
add_subdirectory(tests)
