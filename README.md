# HAZ

[![Build Status](https://travis-ci.org/arkottke/HAZ.svg?branch=develop)](https://travis-ci.org/arkottke/HAZ)

Probabilistic Seismic Hazard Analysis written by Norm Abrahamson.

## Compiling

_HAZ_ is written in FORTRAN and may be compiled with Intel's _ifort_ compiler or
GNU's _gfortran_ compiler.

### Windows

On Windows, the easiest way to get a build system running is via [MSYS2](https://msys2.github.io/). After following the installation
instructions, open an _MinGW_ terminal and issue the following commands to
install the required dependencies. Note, that these commands install the 64-bit
versions. If you prefer to use 32-bit versions replace "x86_64" with "i686".
```
pacman -S git mingw-w64-x86_64-gcc-libgfortran mingw-w64-x86_64-cmake
```
After installing the required dependencies, _HAZ_ can be checked out of Git and
built using the following:
```
git clone https://github.com/abrahamson/HAZ.git HAZ
cd  HAZ
mkdir build
cd build
cmake .. -G "MSYS Makefiles"
make
```
This build will require DLLs provided by _MSYS2_. A static build that includes
the required DLLs, can be built with:
```
cmake .. -G "MSYS Makefiles" -DSTATIC=ON
make
```

The _HAZ_ executable can then be found under `HAZ/build/HAZ.exe`.
