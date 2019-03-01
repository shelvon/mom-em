cmake_minimum_required (VERSION 2.6)
project (mom-em)
enable_language (Fortran)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/../bin/")

######## run system check ########
message("################################")
# check if the environment variable ${MKLROOT} is defined or not
	#message("$ENV{MKLROOT}")
if(DEFINED ENV{MKLROOT})
	# set bash environment variable to CMAKE variable
	set(MKLROOT "$ENV{MKLROOT}")
	message(">> The intel mkl root folder: \n ${MKLROOT}")
else()
	message(">> Please set MKLROOT first!")
endif()

# check HDF5 environment variables
if(DEFINED ENV{HDF5_DIR})
	set(HDF5_DIR "$ENV{HDF5_DIR}")
	message(">> The hdf5 root folder: \n ${HDF5_DIR}")
else()
    message(">> Please set HDF5_DIR first!")
endif()

######## check system variables ########
message("################################")
message(">> Fortran compiler: $ENV{FC}")
message(">> $ENV{FC} version: ${CMAKE_Fortran_COMPILER_VERSION}")
message(">> Host system: ${CMAKE_HOST_SYSTEM_NAME}")
cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)
message(">> Host name: ${HOST_NAME}")
message("################################")

if( (CMAKE_Fortran_COMPILER MATCHES "ifort") OR (CMAKE_Fortran_COMPILER MATCHES "h5fc") )
	set(fortranFlagsRelease "-fpp -qopenmp -parallel -O2 -heap-arrays -I${MKLROOT}/include -I${HDF5_DIR}/include")
	#	set(fortranFlagsDebug "-g -debug full -heap-arrays -warn all -traceback -check all -warn all -I${MKLROOT}/include")
	# set(fortranFlagsDebug "-g -fpp -debug full -heap-arrays -check all -assume none -traceback -I${MKLROOT}/include -I${HDF5_DIR}/include")
	set(fortranFlagsDebug "-g -fpp -debug full -heap-arrays -traceback -I${MKLROOT}/include -I${HDF5_DIR}/include")# assume none causes error for hdf5 compilation
	if(CMAKE_SYSTEM_NAME MATCHES "Linux")
		# linker on Linux
		message("## set the linker for Linux platform")
		set(exeLinkerFlags "-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${HDF5_DIR}/lib -Wl,-rpath -Wl,${HDF5_DIR}/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl -lm")
		# set(exeLinkerFlags "-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${HDF5_DIR}/lib ${HDF5_DIR}/lib/libhdf5hl_fortran.a -lhdf5hl_fortran ${HDF5_DIR}/lib/libhdf5_hl.a -lhdf5_hl ${HDF5_DIR}/lib/libhdf5_fortran.a -lhdf5_fortran ${HDF5_DIR}/lib/libhdf5.a -lhdf5 -lz -ldl -lm")
	elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")	
		# linker on macOS
		message("## set the linker for macOS platform")
		set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${HDF5_DIR}/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl -lm")
		# set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${HDF5_DIR}/lib ${HDF5_DIR}/lib/libhdf5hl_fortran.a -lhdf5hl_fortran ${HDF5_DIR}/lib/libhdf5_hl.a -lhdf5_hl ${HDF5_DIR}/lib/libhdf5_fortran.a -lhdf5_fortran ${HDF5_DIR}/lib/libhdf5.a -lhdf5 -lz -ldl -lm")
		# set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${HDF5_DIR}/lib ${HDF5_DIR}/lib/libhdf5hl_fortran.a ${HDF5_DIR}/lib/libhdf5_hl.a ${HDF5_DIR}/lib/libhdf5_fortran.a ${HDF5_DIR}/lib/libhdf5.a -lhdf5 -lz -ldl -lm")
	endif()
	
	message("################################")
	if(${HOST_NAME} MATCHES "paimenharakka")
	 	message(">>>> $ENV{FC} version is 17")
		set(BUILD_TYPE "debug")
		# set(BUILD_TYPE "release")
	elseif(${HOST_NAME} MATCHES "shelvon.local")
	# elseif(${HOST_NAME} MATCHES "bogon") # hostname@tianjin
	 	message(">>>> $ENV{FC} version is 17")
		set(BUILD_TYPE "debug")
	else()
	 	message(">> $ENV{FC} version is lower than 17")
	 	message(">>>> on cluster, set BUILD_TYPE to release")
	 	message(">>>> modify release and linker flags for ifort 13")
		set(BUILD_TYPE "debug")	
		#set(fortranFlagsRelease "-fpp -openmp -mkl=parallel -O2")
		#set(exeLinkerFlags "-liomp5 -lpthread -lm -ldl")
	endif()
	
elseif(CMAKE_Fortran_COMPILER MATCHES "gfortran")
	message("## gnu fortran compiler")
	set(fortranFlagsRelease "-cpp -fopenmp  -std=gnu -ffree-line-length-0 -DMKL_ILP64 -m64 -I${MKLROOT}/include -O3")
	set(fortranFlagsDebug "-g -cpp -ggdb -std=gnu -DMKL_LP64 -m64 -I${MKLROOT}/include -Wall -Wextra -Warray-temporaries -Wconversion -ffpe-trap=zero,overflow,underflow -finit-real=nan -fbacktrace -fcheck=all -fopenmp -ffree-line-length-0")
	set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
	# set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")

	set(BUILD_TYPE "debug")
else()
	message("## Unsupported fortran compiler!")
endif()

if(BUILD_TYPE MATCHES "release")
 message("Building release")

 # Support OpenMP.
 # May want to add also "-vec-report1" to get vectorization reports.
 set(CMAKE_Fortran_FLAGS "${fortranFlagsRelease}")
elseif(BUILD_TYPE MATCHES "debug")
 message("Building debug")

 set(CMAKE_Fortran_FLAGS "${fortranFlagsDebug}")
elseif(BUILD_TYPE MATCHES "profile")
 message("Building profile")

 set(CMAKE_Fortran_FLAGS "-pg")
else()
	message(FATAL_ERROR "Invalid build type!")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${exeLinkerFlags}")

add_executable(mom
toms644/zbsubs.f
toms644/machcon.f90
fson/fson.f90
fson/fson_path_m.f90
fson/fson_string_m.f90
fson/fson_value_m.f90
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
main.f90
material.f90
mesh.f90
modes.f90
nfields.f90
nfpost.f90
nlbulk.f90
nlsurf.f90
pec.f90
pupil.f90
quad.f90
rwgf.f90
solver_vie.f90
solver.f90
source.f90
srcint.f90
strat.f90
symmetry.f90
sysmat.f90
time.f90
vie.f90
)
