## Makefile
## -----------------------------------------------------
## MPI Version of the Spherical Acoustic Sun Simulator.
## Copyright 2006, Shravan Hanasoge
##
## Hansen Experimental Physics Laboratory
## 455 Via Palou way, Stanford
## CA 94305, USA
## Email: shravan@stanford.edu
## -----------------------------------------------------
##

OBJS=   driver.o        all_modules.o    RHS.o\
        step.o	initialize.o\
	RHS2D.o	derivatives.o

HDF5_FC = mpif90
HDF5_LD = mpif90

FC = h5pfc
LD = h5pfc
 
FFLAGS = -warn all -O2 ## -ffree-line-length-none  -fbounds-check  -fbacktrace -Wall -W
LDFLAGS = 

INCLUDE= fftw3.f

LIBS = -L/opt/local/cfitsio/cfitsio-3.350/lib/ -lcfitsio -lfftw3

COMMAND= sparc.x	

$(COMMAND): $(OBJS)
	$(LD) -o $@ $(FFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm *.o *.mod $(COMMAND)

all_modules.o:
derivatives.o:	all_modules.o
initialize.o:	all_modules.o derivatives.o
RHS.o:	all_modules.o		derivatives.o
driver.o:	all_modules.o	initialize.o	derivatives.o RHS.o	step.o
RHS2D.o:	all_modules.o	derivatives.o	RHS.o
step.o:	all_modules.o	RHS.o	RHS2D.o
