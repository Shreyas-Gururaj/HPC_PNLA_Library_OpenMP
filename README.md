# Parallel Numerical Linear Algebra (Library)

***
## Build Project
This project uses **CMake** as building tool.  
In order to build the project, execute following commands in current directory:

mkdir build  
cd build  
cmake ..  
make  

***
Doxygen was used for documentation.

***
Example usages can be found in "example"(Resolves a Finite Element discretized 2d_poisson system with the given RHS load vector.)

***
CRS format was used for storage and matrix operation routines.

***
Tesing routines are implemented in "tests".

***
Please have a brief look at the <b> "HPC_Lab_WiSe 2020_21.pptx" </b> which is available in this repo to see the scaling in the openmp parallelized implementation of the preconditioned conjugate gradient (PCG) solver tested for liniear equation systems upto 2.5 million unknowns.

***
Please note that the degrees of freedom of the system is equal to num_inner_pts**2.
