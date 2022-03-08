#
#
PROG =	emepctm
###################################################

include Makefile.SRCS

###################################################

# prefered netCDF 4.2.1.1 or later
LIBS = -lnetcdf -lnetcdff
LLIB = -L${NETCDF}/lib
INCL = -I${NETCDF}/include

F90 = mpiifort

# GNU gfortran compiler (version 4.4.3 or later)
#F90FLAGS = -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -O2

# Intel ifort compiler (comment out if gfortran used)
F90FLAGS = -shared-intel -r8 -recursive -O2

###################################################


LDFLAGS = $(F90FLAGS) $(LLIB) -o $(PROG) $(FOBJ) $(INCL) $(LIBS)


.SUFFIXES: $(SUFFIXES)  .f90

.f90.o:
	$(F90) $(F90FLAGS) $(INCL) -c $<


all:  $(PROG)


# Include the dependency-list (created by makedepf90)
include dependencies

$(PROG): $(FOBJ)
	 $(F90) $(LDFLAGS)
#

clean: diskclean

diskclean:
	rm -f $(PROG) *.o *.mod

##########################################################
