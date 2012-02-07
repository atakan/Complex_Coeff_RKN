CC = gcc
FC = gfortran
CFLAGS = -march=core2
ALL_CFLAGS = -Wall -O2 -fpeel-loops -fomit-frame-pointer -ffast-math -DNDEBUG -fopenmp $(CFLAGS)
#ALL_CFLAGS = -Wall -O0 -g -fopenmp $(CFLAGS)
#ALL_CFLAGS = -Wall -O2 $(CFLAGS)

LIBS = -largtable2 -lgsl -lgslcblas -lm

## icc
#CC = /opt/intel/Compiler/11.1/064/bin/ia32/icc
#CFLAGS = -fno-alias -xSSE4.2 -fomit-frame-pointer -O2 -funroll-loops -openmp
#ALL_CFLAGS = -Wall -wd810,869,981,1418,1572 -DNDEBUG $(CFLAGS)
#LIBS = libnpy.a -L/opt/intel/Compiler/11.1/064/lib/ia32 -lgsl -lgslcblas -lm 


#BNAME = err_check_mock
BNAME = err_check
EXE = $(BNAME)

OBJS = $(BNAME).o \
       err_check_initcond.o \
       err_check_dynamics_integrator.o \
       err_check_dynamics_integrator_extras.o \
       err_check_utils.o

#pattern rule to compile object files from C and Fortran files
%.o : %.c Makefile err_check.h
	$(CC) $(ALL_CFLAGS) -c $< -o $@
%.o : %.f Makefile err_check.h
	$(FC) $(ALL_CFLAGS) -c $< -o $@

all: $(EXE)

$(EXE): $(OBJS) Makefile err_check.h
	$(CC) $(ALL_CFLAGS) $(OBJS) -o $(EXE) $(LIBS)

clean:
	rm $(OBJS) $(EXE)
