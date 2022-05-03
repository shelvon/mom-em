cmake_minimum_required (VERSION 2.6)
project (mom-em)
enable_language(Fortran)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/../bin/")

#---- check system environment variables ----
message("---- Checking system environment variables")
message("---- -- $ENV{FC} --version: ${CMAKE_Fortran_COMPILER_VERSION}")# only detected by cmake3 on Linux
message("---- -- Host system: ${CMAKE_HOST_SYSTEM_NAME}")
cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)
message("---- -- Host name: ${HOST_NAME}")

#---- run libraries check ----
message("---- Checking installed libraries")
# check MKL root folder
if(DEFINED ENV{MKLROOT})
  # set bash environment variable to CMAKE variable
  set(MKLROOT "$ENV{MKLROOT}")
  message("---- -- The intel mkl root folder: ${MKLROOT}")
else()
  message("---- -- The MKLROOT is not set, use default folder!")
  set(MKLROOT "/opt/intel/mkl")
endif()

# check HDF5 root folder
if(DEFINED ENV{HDF5_ROOT})
  set(HDF5_ROOT "$ENV{HDF5_ROOT}")
  message("---- -- The hdf5 root folder: ${HDF5_ROOT}")
else()
  message("---- -- The HDF5_ROOT is not set, use default folder!")
  set(HDF5_ROOT "/opt/local")
  message("---- -- The hdf5 root folder: ${HDF5_ROOT}")
endif()

#---- set flags for building and linking
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(FCFLAGS " -cpp -m64")# "-pipe" speed up compiling time in gcc
  set(FCFLAGS_Release " -fopenmp -O3 -ffree-line-length-0")
  
  #set(FCFLAGS_Debug " -g -Wall -Wextra -Warray-temporaries -Wconversion \
  #	-ffpe-trap=zero,overflow,underflow -finit-real=nan -fbacktrace \
  #	-fcheck=all -fopenmp -ffree-line-length-0")
  
  # Be careful with -fcheck=all, it may emit run-time errors as the option ‘pointer’
  # enable generation of run-time checks for pointers and allocatables and make
  # the debugging impossible.
  set(FCFLAGS_Debug " -g -fopenmp -O0 -Wall -Wno-unused-dummy-argument -fbounds-check -fbacktrace -ffree-line-length-0")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(FCFLAGS " -fpp -heap-arrays")
  set(FCFLAGS_Release " -fpp -qopenmp -parallel -O2")
  
  #set(FCFLAGS_Debug " -g -debug full -warn all -traceback -check all")
  set(FCFLAGS_Debug " -g -debug full -traceback -check all")
else()
  message("---- Unsupported fortran compiler!")
endif()

#message("${CMAKE_HOST_SYSTEM_VERSION}")
#return()
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(HDF5_include "lib64/gfortran/modules")
    #set(LDFLAGS_MKL " -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread")
    #set(LDFLAGS_MKL " -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread")
    set(LDFLAGS_MKL " -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -liomp5 -lpthread")
  endif()
  
  set(libdir_MKL " -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64")
  
  set(libdir_HDF5 " -L${HDF5_ROOT}/lib")
  set(LDFLAGS_HDF5 " -Wl,-z,relro -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz")

  # modification in deepin or ubuntu 
  # if(CMAKE_HOST_SYSTEM_VERSION MATCHES "deepin")
  set(HDF5_include "include/hdf5/serial")
  set(libdir_HDF5 " -L${HDF5_ROOT}/lib/x86_64-linux-gnu/hdf5/serial -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial")
  # endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(HDF5_include "include")

  # "-L${MKLROOT}/lib": set path to mkl libraries, followed by
  # "-Wl,-rpath,${MKLROOT}/lib" which is a special linker option on macOS.
  # "-L${MKLROOT}/../compiler/lib" set path to other intel libraries, e.g. libiomp5 (intel openmp)
  # "-Wl,-rpath,${MKLROOT}/../compiler/lib" is again the extra linker command on macOS.
  set(libdir_MKL " -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib")
  
  #set(LDFLAGS_MKL " -lmkl")
  set(LDFLAGS_MKL " -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread")
  #If libiomp5 doesn't work, use libgomp instead and then uncomment the following line
  #set(LDFLAGS_MKL " -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lgomp -lpthread")
  
  set(libdir_HDF5 " -L${HDF5_ROOT}/lib -Wl,-rpath,${HDF5_ROOT}/lib")
  set(LDFLAGS_HDF5 " -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz")
endif()

set(FCFLAGS_MKL " -I${MKLROOT}/include")
set(FCFLAGS_HDF5 " -I${HDF5_ROOT}/${HDF5_include} -D_USE_HDF5")
message("--------------------------------")

if (NOT (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# set to debug mode, if not on cluster 
if((${HOST_NAME} MATCHES "narvi") )
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
  message("---- Building release")

  # Support OpenMP.
  # May want to add also "-vec-report1" to get vectorization reports.
  set(CMAKE_Fortran_FLAGS "${FCFLAGS_Release} ${FCFLAGS} ${FCFLAGS_MKL} ${FCFLAGS_HDF5}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
  message("---- Building debug")

  set(CMAKE_Fortran_FLAGS "${FCFLAGS_Debug} ${FCFLAGS} ${FCFLAGS_MKL} ${FCFLAGS_HDF5}")
  # set(CMAKE_Fortran_FLAGS " -g -cpp -m64 -fopenmp -O0 -Wall")# for test
elseif(CMAKE_BUILD_TYPE MATCHES "profile")
  message("---- Building profile")

  set(CMAKE_Fortran_FLAGS "-pg ${FCFLAGS} ${FCFLAGS_MKL} ${FCFLAGS_HDF5}")
else()
  message(FATAL_ERROR "Invalid build type!")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${libdir_MKL} ${LDFLAGS_MKL} ${libdir_HDF5} ${LDFLAGS_HDF5} -lm -ldl")
#set(CMAKE_EXE_LINKER_FLAGS " -ldl -lm")# for test

message("compile flags: ${CMAKE_Fortran_FLAGS}")
message("linker flags: ${CMAKE_EXE_LINKER_FLAGS}")
# return()
# file(GLOB_RECURSE sources  *.f90 *.f)
# add_executable(mom ${sources})

add_executable(mom
  toms644/zbsubs.f
  toms644/machcon.f90
  fson/src/fson.f90
  fson/src/fson_path_m.f90
  fson/src/fson_string_m.f90
  fson/src/fson_value_m.f90
  aux.f90
  bc.f90
  bessel.f90
  common.f90
  constants.f90
  cs.f90
  data.f90
  diffr.f90
  dipole.f90
  ffields.f90
  focal.f90
  green.f90
  greenprd.f90
  h5_wrapper.f90
  intbase.f90
  interface.f90
  io.f90
  linalg.f90
  material.f90
  mesh.f90
  modes.f90
  nfields.f90
  nfpost.f90
  nlbulk.f90
  nlsurf.f90
  #pec.f90
  pupil.f90
  quad.f90
  rwgf.f90
  solver.f90
  solver_vie.f90
  source.f90
  srcint.f90
  strat.f90
  symmetry.f90
  sysmat.f90
  time.f90
  vie.f90
  main.f90
  )

#	#----doxygen generation
#	#message("${CMAKE_CURRENT_SOURCE_DIR}")
#	#message("${PROJECT_SOURCE_DIR}")
#	#message("${CMAKE_CURRENT_BINARY_DIR}")
#	#set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.in)
#	set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/../doxygen.conf)
#	#configure_file(${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf @ONLY)
#	add_custom_target(doc-doxygen
#	COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
#	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#	COMMENT "Generating html documentation with Doxygen" VERBATIM
#	)
#	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ../doc-doxygen)
