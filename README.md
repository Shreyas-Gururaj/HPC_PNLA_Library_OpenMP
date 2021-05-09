# Parallel Numerical Linear Algebra (Library)

***
##Build Project
This project uses **CMake** as building tool.  
In order to build the project, execute following commands in current directory:

mkdir build  
cd build  
cmake ..  
make  

***
Example usages can be found in "example".

***
Tesing routines will implemented by you in the "tests" file.

***
Please have a brief look at the <b> "HPC_Lab_WiSe 2020_21.pptx" </b> which is available in this repo to see the scaling in the openmp parallelized implementation of the preconditioned conjugate gradient (PCG) solver tested for upto 2.5 million system of equations. It also contains the different.

***
Please note that the degrees of freedom of the system is equal to num_inner_pts**2
