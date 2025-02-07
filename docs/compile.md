# Compiling and running

## Prerequisites

The following must be installed in order to compile and run the code.

- A FORTRAN compiler (GNU/Intel)
- MPI (OpenMPI/MPICH)
- FFTW
- LAPACK (Intel MKL is also an option)
- HDF5

## Compiling 

A `Makefile` is provided along with the code to build the executable on Linux systems. The provided makefile uses the parallel HDF5 wrapper of the GNU FORTRAN compiler `h5pfc`. To use the provided `Makefile`, simply navigate to the directory containing the `Makefile` and type
```
make clean && make
```

Experienced users may be able to modify the provided `Makefile` to compile the code with Intel FORTRAN compilers and MKL.

## Running 

The built executable can be run by navigating to the directory where the executable is located and typing 
```
mpirun -np N ./afid N1 N2
```
or 
```
mpiexec -np ./afid N1 N2
```
where N is the total number of parallel processes divided into N1 times N2 processes along each of the wall-parallel directions. Note that `N` = `N1` $\times$ `N2`.