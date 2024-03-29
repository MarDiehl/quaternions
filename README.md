# Quaternions

A Fortran implementation of a [quaternion](https://en.wikipedia.org/wiki/Quaternion) class

## Getting started

```
export FC=gfortran # or ifort or pgfortran
mkdir build
cd build
cmake .. -DCMAKE_Fortran_COMPILER=$FC
make
./src/test_quaternions
```

## Prerequisites
- Fortran compiler
  - GNU, version 8.0 or newer (`gfortran`)
  - Intel, version 18.0 or newer (`ifort`)
  - PGI, not tested (`pgfortran`)
- cmake, version 3.10 or newer

## License

[GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html)
