CC = icpc #g++-6 #icpc #

F = gfortran-6

HOME=/usr/local

# CFLAGS = -O2 -c -mkl -Wall -Wno-unused-variable -Wno-sign-compare -Wunused-but-set-variable  -ffast-math -funroll-loops -march=native -finline-functions  -std=c++14 -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp

# CFLAGS = -O3 -ffast-math -c -mkl -Wall -Wno-sign-compare -Wunused-but-set-variable -xCORE-AVX2 -msse4 -funroll-loops -finline-functions -mtune=native  -std=c++14 -I/home/hamish/Research/Eigen -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp

#-parallel -qopenmp
CFLAGS = -fast -Ofast -parallel -qopenmp -ffast-math -c -mkl -Wall -Wno-sign-compare -Wunused-but-set-variable -xCORE-AVX2 -msse4 -std=c++14 -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp

# -ipo
CLINK =  -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp -lblas -mkl -ipo

# FFLAGS= -c -I/home/hamish/Research/SHTOOLS-4.0/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/Research/SHTOOLS-4.0/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/Research/SHTOOLS-4.0/lib -lSHTOOLS -Llib -lfftw3 -llapack

FFLAGS= -O3 -c -I/home/hamish/Research/SHTOOLS/modules -m64 -fPIC  -ffast-math -L/home/hamish/Research/SHTOOLS/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas
FLINK =  -lgfortran -L/home/hamish/Research/SHTOOLS/lib -lSHTOOLS -Llib -lfftw3 -llapack

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: ODIS

# ODIS: legendre.o legendreDeriv.o extractSHCoeff.o main.o mathRoutines.o tidalPotentials.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o viscosity.o interpolation.o advection.o solver.o
# 	$(CC) legendre.o legendreDeriv.o extractSHCoeff.o $(FLINK) main.o mathRoutines.o tidalPotentials.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o viscosity.o interpolation.o advection.o solver.o -o ODIS -lgfortran -fopenmp $(CLINK)


ODIS: legendre.o legendreDeriv.o extractSHCoeffGG.o extractSHCoeffLL.o main.o globals.o outFiles.o  mesh.o tidalPotentials.o membraneConstants.o boundaryConditions.o initialConditions.o solver.o timeIntegrator.o drag.o spatialOperators.o temporalOperators.o energy.o sphericalHarmonics.o interpolation.o pressure.o tests.o
	$(CC) $(FLINK) legendre.o legendreDeriv.o extractSHCoeffGG.o extractSHCoeffLL.o $(FLINK) main.o globals.o outFiles.o  mesh.o tidalPotentials.o membraneConstants.o boundaryConditions.o initialConditions.o solver.o timeIntegrator.o drag.o spatialOperators.o temporalOperators.o energy.o sphericalHarmonics.o interpolation.o pressure.o tests.o -o ODIS -lgfortran  $(CLINK)

legendre.o: legendre.f95
	$(F) $(FFLAGS) legendre.f95
#
legendreDeriv.o: legendreDeriv.f95
	$(F) $(FFLAGS) legendreDeriv.f95
#
extractSHCoeffGG.o: extractSHCoeffGG.f95
	$(F) $(FFLAGS) extractSHCoeffGG.f95

extractSHCoeffLL.o: extractSHCoeffLL.f95
	$(F) $(FFLAGS) extractSHCoeffLL.f95

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

tidalPotentials.o: tidalPotentials.cpp
	$(CC) $(CFLAGS) tidalPotentials.cpp

outFiles.o: outFiles.cpp
	$(CC) $(CFLAGS) outFiles.cpp

globals.o: globals.cpp
	$(CC) $(CFLAGS) globals.cpp

mesh.o: mesh.cpp
	$(CC) $(CFLAGS) mesh.cpp

membraneConstants.o: membraneConstants.cpp
	$(CC) $(CFLAGS) membraneConstants.cpp

boundaryConditions.o: boundaryConditions.cpp
	$(CC) $(CFLAGS) boundaryConditions.cpp

initialConditions.o: initialConditions.cpp
	$(CC) $(CFLAGS) initialConditions.cpp

solver.o: solver.cpp
	$(CC)  $(CFLAGS) $(CLINK) solver.cpp

timeIntegrator.o: timeIntegrator.cpp
	$(CC)  $(CFLAGS) $(CLINK) timeIntegrator.cpp

drag.o: drag.cpp
	$(CC)  $(CFLAGS) $(CLINK) drag.cpp

# coriolis.o: coriolis.cpp
# 	$(CC)  $(CFLAGS) $(CLINK) coriolis.cpp

spatialOperators.o: spatialOperators.cpp
	$(CC)  $(CFLAGS) $(CLINK) spatialOperators.cpp

temporalOperators.o: temporalOperators.cpp
	$(CC)  $(CFLAGS) $(CLINK) temporalOperators.cpp

energy.o: energy.cpp
	$(CC)  $(CFLAGS) $(CLINK) energy.cpp

sphericalHarmonics.o: sphericalHarmonics.cpp
	$(CC) $(CFLAGS) sphericalHarmonics.cpp

interpolation.o: interpolation.cpp
	$(CC) $(CFLAGS) interpolation.cpp

pressure.o: pressure.cpp
	$(CC)  $(CFLAGS) $(CLINK) pressure.cpp

tests.o: tests.cpp
	$(CC)  $(CFLAGS) $(CLINK) tests.cpp

clean:
	rm -r *o ODIS NorthVelocity EastVelocity Displacement Grid Energy
