#-----------------------------------------------------
# Makefile to compile the Alpha_driver system on DROP
#-----------------------------------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
## FIXED Hardware parameters
#IDENTITY="Cray XC30 (edison)"
#CORES_PER_NODE=24
#NUMAS_PER_NODE=2
#
## Compilers and flags
#FC     = ftn -e m -J ${GACODE_ROOT}/modules
#F77    = ${FC}
#FOMP   =
#FMATH  = -s real64
#FOPT   = -O3 -hfp3
#FDEBUG = -eD -Ktrap=fp -m 1 -R bcp
#
## System math libraries
#LMATH = ${FFTW_POST_LINK_OPTS} ${FFTW_INCLUDE_OPTS}
##LMATH = -L${FFTW_DIR}/lib -Wl,-rpath=${FFTW_DIR}/lib \
##        -lfftw3_threads -lfftw3f_threads -lfftw3 -lfftw3f
#
## NetCDF
#NETCDF = ${NETCDF_DIR}/lib/libnetcdff.a ${NETCDF_DIR}/lib/libnetcdf.a
#NETCDF_INC = ${NETCDF_DIR}/include
#
## Archive 
#ARCH = ar cr

export EXTRA_LIBS = \
        ${GACODE_ROOT}/f2py/expro/expro_lib.a \
        ${GACODE_ROOT}/f2py/geo/geo_lib.a \
        ${GACODE_ROOT}/shared/math/math_lib.a \

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
#   FFLAGS=${FDEBUG}
endif

EXEC=Alpha_driver

LLIB=Alpha_lib

OBJECTS = Alpha_use_input.o \
          Alpha_use_output.o \
          Alpha_use_other.o \
          Alpha_use_sav_diffEP.o \
          Alpha_use_sav_QLdiffEP.o \
          Alpha_allocate.o \
          Alpha_read_input.o \
          Alpha_comp_eq_plasma.o \
          Alpha_comp_alpha_slowing.o \
          Alpha_transport.o \
          Alpha_time_dep_transport.o \
          Alpha_diffusivity.o \
          Alpha_QLdiffusivity.o \
          Alpha_write_output.o \
          Alpha_mainsub.o 

.SUFFIXES : .o .f90 .f .F

all: $(LLIB).a $(EXEC)
	rm ${GACODE_ROOT}/modules/ALPHA*.mod

$(EXEC): $(LLIB).a $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o $(LLIB).a $(EXTRA_LIBS)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o  $(EXEC) $(LLIB).a
