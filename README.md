# diffev
Differential evolution algorithm

The program read in the snapshot configurations of solute-solvent system
and corresponding vibratinoal frequency shifts and optimize the polynomial
model of the vibrational frequency shifts[1] using the differential evolution
algorithm. 

## 'source' directory
The source directory contains seven Fortran code files. The six files excluding
'singular_par_module.f90' were taken from online and modified by me according to 
my need. The file 'singular_par_module.f90' implements the calculation of the 
polynomial terms from the input configurations and implements singular value 
decomposition using the LAPACK routine GESVD provided by the Intel Math Kernel 
Library. The following command makes the executable 'de':

>ifort -o de singular_par_module.f90 derived.f90 randperm.f90 left_win.f90 objfun.f90 deopt.f90 rundeopt.f90 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include  ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

Here $MKLROOT represents the location of the Intel Math Kernel Library.

## 'example' directory
Three input files (input.de, svd.in, even.dat) are used.
The output 'param.dat' contains the optimized coefficients of the polynomial terms.
In the last part of the stout output shows the relative contributions of each order of the terms, which are not yet normalized.

# Reference
[1] Kijeong Kwac and Minhaeng Cho, J. Chem. Phys. 151, 134112 (2019)
