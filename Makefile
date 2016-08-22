CC=g++

HOME=/usr/local

#-I/usr/local/hdf5/include -lhdf5

CFLAGS= -Ofast -c -Wall  -std=c++11

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: ODIS

ODIS: test.o main.o mathRoutines.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o solver.o
	$(CC) test.o -lgfortran main.o mathRoutines.o outFiles.o globals.o mesh.o field.o depth.o mass.o energy.o solver.o -o ODIS

test.o: test.f90
	gfortran -c test.f90

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
	$(CC) $(CFLAGS) solver.cpp

clean:
	rm -r *o ODIS NorthVelocity EastVelocity Displacement Grid Energy
