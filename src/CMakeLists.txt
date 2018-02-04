cmake_minimum_required (VERSION 2.6)
project (for90-mom2)
enable_language (Fortran)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/../bin/")

######## run system check ########
message("################################")
# check if the environment variable ${MKLROOT} is defined or not
	#message("$ENV{MKLROOT}")
if(DEFINED ENV{MKLROOT})
	# set bash environment variable to CMAKE variable
	set(MKLROOT "$ENV{MKLROOT}")
	message("## intel mkl root folder: \n ${MKLROOT}")
else()
	message("## please set MKLROOT first!")
endif()

#message($ENV{FC})
if(CMAKE_Fortran_COMPILER MATCHES "ifort")
	set(fortranFlagsRelease "-fpp -openmp -parallel -O2")
	#	set(fortranFlagsDebug "-g -debug full -heap-arrays -warn all -traceback -check all -I${MKLROOT}/include")
	set(fortranFlagsDebug "-g -debug full -heap-arrays -assume none -traceback -I${MKLROOT}/include")
	if(CMAKE_SYSTEM_NAME MATCHES "Linux")
		# linker on Linux
		message("## set the linker for Linux platform")
		set(exeLinkerFlags "-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
	elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")	
		# linker on macOS
		message("## set the linker for macOS platform")
		set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
	endif()
	
elseif(CMAKE_Fortran_COMPILER MATCHES "gfortran")
#elseif(CMAKE_Fortran_COMPILER MATCHES "gfortran-mp-4.8")
	message("## gnu fortran compiler")
	set(fortranFlagsRelease "-fopenmp  -std=f2008 -ffree-line-length-0 -DMKL_ILP64 -m64 -I${MKLROOT}/include -O3")
	set(fortranFlagsDebug "-g -ggdb -DMKL_LP64 -m64 -I${MKLROOT}/include -Wall -Wextra -Warray-temporaries -Wconversion -ffpe-trap=zero,overflow,underflow -finit-real=nan -fbacktrace -fcheck=all -fopenmp -std=f2008 -ffree-line-length-0")
	set(exeLinkerFlags "-L${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
	# set(exeLinkerFlags "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")

else()
	message("## Unsupported fortran compiler!")
endif()
message("################################")


set(BUILD_TYPE "debug")

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
aux.f90
constants.f90
ffields.f90
focal.f90
interface.f90
mesh.f90
nlsurf.f90
solver.f90
modes.f90
strat.f90
vie.f90
bc.f90
cs.f90
intbase.f90
nfields.f90
pec.f90
solver_vie.f90
symmetry.f90
bessel.f90
diffr.f90
green.f90
linalg.f90
nfpost.f90
quad.f90
source.f90
pupil.f90
sysmat.f90
common.f90
dipole.f90
greenprd.f90
main.f90
nlbulk.f90
rwgf.f90
srcint.f90
time.f90
)
