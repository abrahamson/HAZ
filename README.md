# HAZ

[![Build Status](https://travis-ci.org/abrahamson/HAZ.svg?branch=master)](https://travis-ci.org/abrahamson/HAZ)

Probabilistic Seismic Hazard Analysis written by Norm Abrahamson.

## Compiling

_HAZ_ is written in FORTRAN and may be compiled with Intel's _ifort_ compiler or
GNU's _gfortran_ compiler.

### Windows

The following set of instructions use 64-bit versions of the installers and
packages (named *x86_64*), which should be appropriate for modern operating
systems. If you have an older, 32-bit, you should use the *i686* version. The
sytem type can be found under the "System type:" field of System Information.
On Windows 7, this can be accessed by: Start Menu -> Right Click "My Computer"
-> Select Properties. On Windows 10, this can be accessed by opeining the Start
Menu, then typing "About your PC.

1. Download and install MSYS2 following the directions at the
   [website](http://www.msys2.org/). Note, that both *x86_64* and *i686*
   installers are available. MSYS2 provides a bash shell, revision control
   systems and the like for building native Windows applications using
   MinGW-w64 toolchains.

2. The installation process starts a MSYS bash shell, by default. Using this
   shell, we will update the packages to the latest version with: `pacman
   -Syu`. This update typically requires that the MSYS shell be closed after
   the initial update. Close the shell, create a new shell (named "MSYS2
   MSYS"), and run `pacman -Syu` once more. Now the system is completely
   up-to-date.

3. Install the dependences for HAZ45 with

   ```
   pacman -S git mingw-w64-x86_64-gcc-fortran make mingw-w64-x86_64-cmake 
   ```

   This command will install:
     * `git` for working with the code repository, 
     * `gfortran` compiler, 
     * `cmake` and `make` for automating the build, and
     * the required supporting files.
  
4. Open a "MSYS MinGW 64-bit" shell. This is different from the MSYS shell and
   modifies the `PATH` environment to include 64-bit specific programs -- which
   were just installed. You can check that you have the appropriate system by
   checking the versions of the programs that are needed:

   ```
   akottke@ASBD1028 MINGW64 ~
   $ echo $MINGW_PREFIX
   /mingw64
   akottke@ASBD1028 MINGW64 ~
   $ make --version
   GNU Make 4.2.1
   Built for x86_64-pc-msys
   Copyright (C) 1988-2016 Free Software Foundation, Inc.
   License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
   This is free software: you are free to change and redistribute it.
   There is NO WARRANTY, to the extent permitted by law.

   akottke@ASBD1028 MINGW64 ~
   $ cmake --version
   cmake version 3.9.2

   CMake suite maintained and supported by Kitware (kitware.com/cmake).
   ```

   If you get errors at this step, stop and review the steps to ensure that
   `make` and `cmake` are properly installed.

5. Use `git` to checkout the `HAZ` source code and start the build.

   ```
   git clone https://github.com/abrahamson/HAZ.git HAZ
   cd  HAZ
   mkdir build
   cd build
   cmake .. -G "MSYS Makefiles"
   make
   ```

   Let's review this process. The repository is cloned into the "HAZ" directory
   with:

   ```
   akottke@ASBD1028 MINGW64 ~
   $ git clone https://github.com/abrahamson/HAZ.git HAZ
   Cloning into 'HAZ'...
   remote: Counting objects: 6156, done.
   remote: Total 6156 (delta 0), reused 0 (delta 0), pack-reused 6156
   Receiving objects: 100% (6156/6156), 7.21 MiB | 2.44 MiB/s, done.
   Resolving deltas: 100% (4351/4351), done.
   Checking out files: 100% (1446/1446), done.
   ```

   Next, we move into the "HAZ" directory, and create a "build" directory.
   Creating this directory separates the files built during the compliation
   from the human created source code.

   ```
   cd  HAZ
   mkdir build
   cd build
   ```

   Next, we run `cmake` to create `make` files that are using in the compiling
   and linking of the source code. This process takes a moment for `cmake` to
   determine the configuration and capabilities of the system.

   ```
   akottke@ASBD1028 MINGW64 ~/HAZ/build
   $  cmake .. -G "MSYS Makefiles"
   -- The C compiler identification is GNU 7.2.0
   -- The CXX compiler identification is GNU 7.2.0
   -- Check for working C compiler: C:/msys64/mingw64/bin/gcc.exe
   -- Check for working C compiler: C:/msys64/mingw64/bin/gcc.exe -- works
   -- Detecting C compiler ABI info
   -- Detecting C compiler ABI info - done
   -- Detecting C compile features
   -- Detecting C compile features - done
   -- Check for working CXX compiler: C:/msys64/mingw64/bin/g++.exe
   -- Check for working CXX compiler: C:/msys64/mingw64/bin/g++.exe -- works
   -- Detecting CXX compiler ABI info
   -- Detecting CXX compiler ABI info - done
   -- Detecting CXX compile features
   -- Detecting CXX compile features - done
   -- The Fortran compiler identification is GNU 7.2.0
   -- Check for working Fortran compiler: C:/msys64/mingw64/bin/gfortran.exe
   -- Check for working Fortran compiler: C:/msys64/mingw64/bin/gfortran.exe  -- works
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Checking whether C:/msys64/mingw64/bin/gfortran.exe supports Fortran 90
   -- Checking whether C:/msys64/mingw64/bin/gfortran.exe supports Fortran 90 -- yes
   --> Fortran compiler: gfortran.exe
   --> CMAKE_Fortran_COMPILER full path: C:/msys64/mingw64/bin/gfortran.exe
   --> CMAKE_Fortran_FLAGS: -ffixed-line-length-132 -fno-automatic
   -- Configuring done
   -- Generating done
   -- Build files have been written to: C:/msys64/home/akottke/HAZ/build
   ```

   Finally, the code is compiled with `make`:

   ```
   akottke@ASBD1028 MINGW64 ~/HAZ/build
   $ make
   Scanning dependencies of target HAZ
   [  2%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/ATTEN.F.obj
   C:/msys64/home/akottke/HAZ/HAZ files/ATTEN.F:9141:132:

           pause
                                                                                                                                       1
   Warning: Deleted feature: PAUSE statement at (1)
   [  5%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/CLDIST.F.obj
   [  8%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/Directivity.f.obj
   [ 11%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/Directivity_DPP.f.obj
   [ 14%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/Directivity_JWL.f.obj
   [ 17%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/FLT_SEG.F.obj
   [ 20%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/FRISK_IO.F.obj
   C:/msys64/home/akottke/HAZ/HAZ files/FRISK_IO.F:267:132:

            pause
                                                                                                                                       1
   Warning: Deleted feature: PAUSE statement at (1)
   C:/msys64/home/akottke/HAZ/HAZ files/FRISK_IO.F:103:39:

          call CheckDim ( nProb, MAX_Prob, 'MAX_PROB ' )
                                          1
   Warning: Character length of actual argument shorter than of dummy argument 'name' (9/10) at (1) [-Wargument-mismatch]
   [ 23%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/GC2.f.obj
   [ 26%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/HAZ_INT.F.obj
   [ 29%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/HAZ_MAIN2.f.obj
   [ 32%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/INTERP.F.obj
   [ 35%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/MAGPROB.F.obj
   [ 38%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/MATRIX.F.obj
   [ 41%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/MED_ATT.F.obj
   [ 44%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/NGA2008.F.obj
   [ 47%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/NGA2008_TRAdjusted_V3.F.obj
   [ 50%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/NGAWest2.F.obj
   [ 52%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/NGAWest2_Vert.F.obj
   [ 55%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/NPROB.F.obj
   [ 58%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/SIGMA.f.obj
   [ 61%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/atten_eq_tw.f.obj
   [ 64%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/atten_set2.f.obj
   [ 67%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/attenper.f.obj
   [ 70%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/bcHydro_subduction.f.obj
   [ 73%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/calcDepthProb.f.obj
   [ 76%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/determ.f.obj
   [ 79%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/directivity_bayless.f.obj
   [ 82%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/eus_atten2.f.obj
   [ 85%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/nRuptures.f.obj
   [ 88%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/readFlt.f.obj
   [ 91%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/readGrid.f.obj
   [ 94%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/setRates.f.obj
   [ 97%] Building Fortran object CMakeFiles/HAZ.dir/HAZ_files/soilAMp.f.obj
   [100%] Linking Fortran executable HAZ.exe
   [100%] Built target HAZ
   ```
6. `HAZ` can then be run with `HAZ.exe`:

   ```
   akottke@ASBD1028 MINGW64 ~/HAZ/build
   $ ./HAZ.exe
    *********************************
    *          Hazard Code          *
    *        Release HAZ45.2        *
    *      Tagged Jan 2, 2017       *
    *********************************

    Enter number of cases to run in batch mode.
             (For a single run enter 0)
   ```

   Note that you have to run `HAZ` from the "MSYS MinGW 64-bit" shell.
   Alternatively, you can figure the build to compile a static version with:

   ```
   cmake .. -G "MSYS Makefiles" -DSTATIC=ON
   make
   ```

7. If you modify the source code, the program can recompiled using with `make`
   -- none of the previous steps are required. However, if you add additional
   files, then the `cmake` step must be repeated. All files in the "HAZ files"
   directory ending in ".f" or ".F" are compiled and linked.

