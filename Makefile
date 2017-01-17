CC = g++-6

F = gfortran-6

HOME=/usr/local

#-I/usr/local/hdf5/include -lhdf5 -L/usr/local/hdf5/lib

CFLAGS= -O2 -c -Wall  -std=c++11 -L/usr/local/hdf5/lib -I/usr/local/hdf5/include -lhdf5 -lhdf5_cpp

CLINK = -L/usr/local/hdf5/lib -I/usr/local/hdf5/include -lhdf5 -lhdf5_cpp

FFLAGS= -c -I/source/SHTOOLS/modules  -m64 -fPIC -O2 -ffast-math -L/source/SHTOOLS/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas

FLINK = -lgfortran -L/source/SHTOOLS/lib -lSHTOOLS -L/usr/local/lib -lfftw3

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: ODIS_TEST

ODIS_TEST: extractSHCoeff.o main.o mathRoutines.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o solver.o
	$(CC) extractSHCoeff.o $(FLINK) main.o mathRoutines.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o solver.o -o ODIS_TEST -lgfortran $(CLINK)

extractSHCoeff.o: extractSHCoeff.f95
	$(F) $(FFLAGS) extractSHCoeff.f95

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

mathRoutines.o: mathRoutines.cpp
	$(CC) $(CFLAGS) mathRoutines.cpp

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
