CC = g++-6

F = gfortran-6

HOME=/usr/local

#-mavx2 -mfma

CFLAGS= -O3 -g -c -Wall -Wno-unused-variable -Wunused-but-set-variable -std=c++11 -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp -fopenmp

CLINK = -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp

FFLAGS= -c -I/home/hamish/Research/SHTOOLS-4.0/modules  -m64 -fPIC -O2 -ffast-math -L/home/hamish/Research/SHTOOLS-4.0/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas

FLINK = -lgfortran -L/home/hamish/Research/SHTOOLS-4.0/lib -lSHTOOLS -Llib -lfftw3

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: ODIS

ODIS: legendre.o legendreDeriv.o extractSHCoeff.o main.o mathRoutines.o tidalPotentials.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o solver.o
	$(CC) legendre.o legendreDeriv.o extractSHCoeff.o $(FLINK) main.o mathRoutines.o tidalPotentials.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o solver.o -o ODIS -lgfortran -fopenmp $(CLINK)

legendre.o: legendre.f95
	$(F) $(FFLAGS) legendre.f95

legendreDeriv.o: legendreDeriv.f95
	$(F) $(FFLAGS) legendreDeriv.f95

extractSHCoeff.o: extractSHCoeff.f95
	$(F) $(FFLAGS) extractSHCoeff.f95

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

mathRoutines.o: mathRoutines.cpp
	$(CC) $(CFLAGS) mathRoutines.cpp

tidalPotentials.o: tidalPotentials.cpp
	$(CC) $(CFLAGS) tidalPotentials.cpp

outFiles.o: outFiles.cpp
	$(CC) $(CFLAGS) outFiles.cpp

globals.o: globals.cpp
	$(CC) $(CFLAGS) globals.cpp

mesh.o: mesh.cpp
	$(CC) $(CFLAGS) mesh.cpp

field.o: field.cpp
	$(CC) $(CFLAGS) field.cpp

depth.o: depth.cpp
	$(CC) $(CFLAGS) depth.cpp

mass.o: mass.cpp
	$(CC) $(CFLAGS) mass.cpp

energy.o: energy.cpp
	$(CC) $(CFLAGS) energy.cpp

solver.o: solver.cpp
	$(CC)  $(CFLAGS) $(CLINK) solver.cpp

clean:
	rm -r *o ODIS NorthVelocity EastVelocity Displacement Grid Energy
