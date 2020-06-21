# Quaternions

A Fortran implementation of a quaternion class

## Getting started

```
export FC=gfortran # or ifort
mkdir build
cd build
cmake .. -DCMAKE_Fortran_COMPILER=$FC
make
./src/test_quaternions
```

## Prerequisites
- Fortran compiler
  - GNU (`gfortran`), > 8.1
  - Intel (`ifort`) > 18.1
  - PGI (`pgfortran`), not tested but might work
- cmake > 3.10
