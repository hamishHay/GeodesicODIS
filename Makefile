CC = g++ #icpc #

F = gfortran

SRC_DIR = src
OBJ_DIR = obj

print-%  : ; @echo $* = $($*)

CSRC = $(wildcard $(SRC_DIR)/*.cpp)
COBJ = $(CSRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# $(COBJ): src/%.o : src/%.c

FSRC = #$(wildcard $(SRC_DIR)/*.f95)
FOBJ = #$(FSRC:$(SRC_DIR)/%.f95=$(OBJ_DIR)/%.o)

EXE = ODIS


# HOME MACHINE
CFLAGS = -Ofast -ffast-math -c -lzlib -Wall -Iinclude -Wno-sign-compare -Wunused-but-set-variable  -msse4 -std=c++14 -I/home/hay/Research/eigen -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp
CLINK =  -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial  -Iinclude -I/home/hay/Research/eigen -lhdf5_cpp -lhdf5_serial #-lblas #-mkl #-ipo

FFLAGS= #-c -I/home/hamish/SHTOOLS/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/SHTOOLS/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
FLINK = # -lgfortran -L/home/hamish/SHTOOLS/lib -lSHTOOLS -Llib -lfftw3 -llapack -lgfortran

# WORK MACHINE
# CFLAGS = -Ofast -ffast-math -c -Wall -Iinclude -Wno-sign-compare -Wunused-but-set-variable  -msse4 -std=c++14 -I/home/hamish/Research/Eigen -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp
# CLINK =  -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial  -Iinclude -I/home/hamish/Eigen -lhdf5 -lhdf5_cpp -lblas #-mkl #-ipo

# FFLAGS= -c -I/home/hamish/Research/SHTOOLS/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/Research/SHTOOLS/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/Research/SHTOOLS/lib -lSHTOOLS -Llib -lfftw3 -llapack -lgfortran

# gfortran -I/home/hamish/SHTOOLS/modules -m64 -fPIC -O3 -std=f2003 -ffast-math -L/home/hamish/SHTOOLS/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas

# CFLAGS = -fast -Ofast -parallel -qopenmp -ffast-math -c -mkl -Wall -Iinclude -Wno-sign-compare -Wunused-but-set-variable -xCORE-AVX2 -msse4 -std=c++14 -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp
# CLINK =  -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial  -Iinclude -lhdf5 -lhdf5_cpp -lblas -mkl #-ipo

# FFLAGS= -c -I/home/hamish/Research/Research/SHTOOLS-4.0/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/Research/Research/SHTOOLS-4.0/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/Research/Research/SHTOOLS-4.0/lib -lSHTOOLS -Llib -lfftw3 -llapack


all: $(EXE)

$(EXE): $(FOBJ) $(COBJ)
	$(CC) $(FLINK) $^ -o $@  $(FLINK) $(CLINK) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< $ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f95
	$(F) $(FFLAGS) $< $ -o $@


clean:
	rm -r $(FOBJ) $(COBJ) ODIS OUTPUT.txt ERROR.txt
