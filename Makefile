CC = g++ #icpc #

F = gfortran

SRC_DIR = src
OBJ_DIR = obj
CONSTS_DIR = constants
OUTPUT_DIR = $(shell pwd)

print-%  : ; @echo $* = $($*)

CSRC = $(wildcard $(SRC_DIR)/*.cpp)
COBJ = $(CSRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
#  -DNDEBUG
# $(COBJ): src/%.o : src/%.c

FSRC = #$(wildcard $(SRC_DIR)/*.f95)
FOBJ = #$(FSRC:$(SRC_DIR)/%.f95=$(OBJ_DIR)/%.o)

EXE = $(OUTPUT_DIR)/ODIS


# HOME MACHINE
# CFLAGS = -Ofast -ffast-math -c -Wall -Iinclude -Wno-sign-compare -Wunused-but-set-variable  -msse4 -std=c++14 -I/home/hamish/Eigen -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp
# CLINK =  -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial  -Iinclude -I/home/hamish/Eigen -lhdf5 -lhdf5_cpp -lblas #-mkl #-ipo

# FFLAGS= -c -I/home/hamish/SHTOOLS/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/SHTOOLS/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/SHTOOLS/lib -lSHTOOLS -Llib -lfftw3 -llapack -lgfortran

# WORK MACHINE

# CFLAGS = -Ofast -I/usr/local/opt/libomp/include -c -Wall -Iinclude -I$(CONSTS_DIR)  -DNDEBUG -Wno-sign-compare -Wunused-but-set-variable -Xclang -fopenmp -msse4 -lstdc++ -std=c++17 -I/usr/local/ -I/Users/hamishhay/Desktop/Research/ -I/usr/local/Cellar/hdf5/1.12.2_2/include -march=native
# # CFLAGS = -O3 -c -Wall -Iinclude -DNDEBUG -diag-disable=10441 -Wno-sign-compare -Wunused-but-set-variable -fopenmp -msse4 -std=c++17 -I/usr/local/ -I/Users/hamishhay/Desktop/Research/ -I/usr/local/Cellar/hdf5/1.12.2_2/include -march=native #-L/usr/include/hdf5/serial -I/usr/include/hdf5/serial #-lhdf5 -lhdf5_cpp
# CLINK = -L/usr/local/opt/libomp/lib -L/usr/local/Cellar/hdf5/1.12.2_2/lib/ -L/usr/local/Cellar/hdf5/1.12.2_2/include -lhdf5 -lomp -lhdf5_cpp -lblas #-mkl #-ipo

# FFLAGS= -c -I/home/hamish/Research/SHTOOLS/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/Research/SHTOOLS/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/Research/SHTOOLS/lib -lSHTOOLS -Llib -lfftw3 -llapack -lgfortran

# FRAMEWORK 
CFLAGS = -c -O3 -Wall -Iinclude -Wno-sign-compare -I$(CONSTS_DIR) -Wunused-but-set-variable -msse4 -std=c++17 -I/usr/local/ -I/usr/local/include/hdf5/ -I/usr/include/hdf5/serial -march=native #-lhdf5 -lhdf5_cpp -Ofast -ffast-math 
CLINK =  -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -I/usr/local/include/hdf5 -lhdf5 -lhdf5_cpp  -lblas #-mkl #-ipo 

# FFLAGS= -c -I/home/hamish/Research/SHTOOLS/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/Research/SHTOOLS/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/Research/SHTOOLS/lib -lSHTOOLS -Llib -lfftw3 -llapack -lgfortran

# gfortran -I/home/hamish/SHTOOLS/modules -m64 -fPIC -O3 -std=f2003 -ffast-math -L/home/hamish/SHTOOLS/lib -lSHTOOLS -L/usr/local/lib -lfftw3 -lm -llapack -lblas

# CFLAGS = -fast -Ofast -parallel -qopenmp -ffast-math -c -mkl -Wall -Iinclude -Wno-sign-compare -Wunused-but-set-variable -xCORE-AVX2 -msse4 -std=c++14 -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp
# CLINK =  -L/usr/include/hdf5/serial -I/usr/include/hdf5/serial  -Iinclude -lhdf5 -lhdf5_cpp -lblas -mkl #-ipo

# FFLAGS= -c -I/home/hamish/Research/Research/SHTOOLS-4.0/modules -m64 -fPIC -Ofast -ffast-math -L/home/hamish/Research/Research/SHTOOLS-4.0/lib  -L/usr/local/lib -lfftw3 -lm -llapack -lblas
# FLINK =  -lgfortran -L/home/hamish/Research/Research/SHTOOLS-4.0/lib -lSHTOOLS -Llib -lfftw3 -llapack


all: $(EXE)


$(EXE): $(COBJ)
	$(CC) $^ -o $@  $(CLINK) 

# $(EXE): $(FOBJ) $(COBJ)
# $(CC) $(FLINK) $^ -o $@  $(FLINK) $(CLINK) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c -MMD $< $ -o $@

# $(OBJ_DIR)/%.o: $(SRC_DIR)/%.f95
# 	$(F) $(FFLAGS) $< $ -o $@

-include $(OBJ_DIR)/*.d

clean:
	rm -r $(FOBJ) $(COBJ) $(OBJ_DIR)/*.d ODIS OUTPUT.txt ERROR.txt
