CC=g++

CFLAGS= -O3 -c -Wall  -std=c++11

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: ODIS

ODIS: main.o mathRoutines.o outFiles.o globals.o mesh.o field.o mass.o energy.o solver.o
	$(CC) main.o mathRoutines.o outFiles.o globals.o mesh.o field.o mass.o energy.o solver.o -o ODIS

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

mass.o: mass.cpp
	$(CC) $(CFLAGS) mass.cpp

energy.o: energy.cpp
	$(CC) $(CFLAGS) energy.cpp

solver.o: solver.cpp
	$(CC) $(CFLAGS) solver.cpp

clean:
	rm -r *o ODIS *.txt NorthVelocity EastVelocity Displacement Grid
