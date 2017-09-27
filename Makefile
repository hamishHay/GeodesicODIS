CC = g++-6

F = gfortran-6

HOME=/usr/local

#-mavx2 -mfma

CFLAGS = -O3 -c -Wall -Wno-unused-variable -Wno-sign-compare -Wunused-but-set-variable  -ffast-math -funroll-loops -finline-functions  -std=c++14 -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp

CLINK = -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp

FFLAGS= -c -I/home/hamish/Research/SHTOOLS-4.0/modules -m64 -fPIC -O3 -ffast-math -L/home/hamish/Research/SHTOOLS-4.0/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas

FLINK =  -lgfortran -L/home/hamish/Research/SHTOOLS-4.0/lib -lSHTOOLS -Llib -lfftw3 -llapack

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: ODIS

# ODIS: legendre.o legendreDeriv.o extractSHCoeff.o main.o mathRoutines.o tidalPotentials.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o viscosity.o interpolation.o advection.o solver.o
# 	$(CC) legendre.o legendreDeriv.o extractSHCoeff.o $(FLINK) main.o mathRoutines.o tidalPotentials.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o viscosity.o interpolation.o advection.o solver.o -o ODIS -lgfortran -fopenmp $(CLINK)


ODIS: legendre.o legendreDeriv.o extractSHCoeffGG.o extractSHCoeffLL.o main.o tidalPotentials.o outFiles.o globals.o mesh.o solver.o timeIntegrator.o drag.o coriolis.o spatialOperators.o temporalOperators.o energy.o sphericalHarmonics.o interpolation.o pressure.o
	$(CC) $(FLINK) legendre.o legendreDeriv.o extractSHCoeffGG.o extractSHCoeffLL.o $(FLINK) main.o tidalPotentials.o outFiles.o globals.o mesh.o solver.o timeIntegrator.o drag.o coriolis.o spatialOperators.o temporalOperators.o energy.o sphericalHarmonics.o interpolation.o pressure.o -o ODIS -lgfortran $(CLINK)

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

solver.o: solver.cpp
	$(CC)  $(CFLAGS) $(CLINK) solver.cpp

timeIntegrator.o: timeIntegrator.cpp
	$(CC)  $(CFLAGS) $(CLINK) timeIntegrator.cpp

drag.o: drag.cpp
	$(CC)  $(CFLAGS) $(CLINK) drag.cpp

coriolis.o: coriolis.cpp
	$(CC)  $(CFLAGS) $(CLINK) coriolis.cpp

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

clean:
	rm -r *o ODIS NorthVelocity EastVelocity Displacement Grid Energy
