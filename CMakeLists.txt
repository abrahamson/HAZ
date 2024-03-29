# Adapted from: https://cmake.org/Wiki/CMakeForFortranExample
# CMake project file for HAZ

cmake_minimum_required (VERSION 2.8.12)
project (HAZ)
enable_language (Fortran)

option(STATIC "Enable static linking for binary distribution")

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the HAZ directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran

  # Flags required for source formatting
  set (CMAKE_Fortran_FLAGS "-ffixed-line-length-132 -fno-automatic")
  # Static compliation
  if (STATIC)
      message("--> Static build.")
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static")
  endif (STATIC)
  # Optimization flags
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3") # -ffpe-trap=invalid,zero,overflow,underflow,denormal -ffpe-summary=none")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

message ("--> Fortran compiler: " ${Fortran_COMPILER_NAME})
message ("--> CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message ("--> CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})

# build executables
file (GLOB SOURCES src/*.f)
set (TARGET "HAZ")

add_executable (${TARGET} ${SOURCES})
set_target_properties(${TARGET} PROPERTIES LINKER_LANGUAGE Fortran)
