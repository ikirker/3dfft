####################################################
#  Makefile for parallel 3D FFT benchmark.         #
#                                                  #
#  By Ian Kirker.                                  #
####################################################

# Last revision 22/8/08

# Syntax: make LIB=[fftw3|fftw2|essl|mkl|acml] 
#			   CC=[gcc|pgcc|xlc|xlc_bg] 
#              MPICC=any 
#              SYSTEM=[generic|Antimony|ness|hector|hpcx|eddie|bluegene|marenostrum]
#              fft


# File variables.
SRC=A2A3D.c  \
	comms.c  \
	dataOps.c \
	decomposition.c \
	libDefs.c \
	main.c \
	options.c \
	performLocalTranspose.c \
	validateParameters.c 
	
OBJ=$(SRC:.c=.o)
HEADERS=$(SRC:.c=.h)

# Compilers - CC is underlying compiler, only used to determine correct optimisation flags.
MPICC=mpicc
CC=gcc

# Compiler-specific optimisation/general use flags.
gccOPTFLAGS= -qfast
pgccOPTFLAGS= -fast -Mipa=fast
xlc_bgOPTFLAGS= -O3 -qhot -qarch=440d -qtune=440 -qlanglvl=stdc99
xlcOPTFLAGS= -q64 -O5 -qlanglvl=stdc99

CFLAGS=$($(CC)OPTFLAGS)

############################################
# Flags to find and link fft libraries.    #
############################################

# Defaults
LIB=fftw3
SYSTEM=generic

#If you have FFTW3.2a3 or better, uncomment extra lines
fftw3_on_generic_flags= \
	-lfftw3 #\
	#-lfftw3_mpi \
	#-DHAS_AUTO
# If you don't have fftw2 mpi, remove -DHAS_AUTO 
#  and the fftw_mpi link.
fftw2_on_generic_flags= \
	-lfftw \
	-lfftw_mpi \
	-DHAS_AUTO
# Use this option instead if you've compiled fftw2 with type
#  prefixes.
#fftw2_on_generic_flags= \
#	-ldfftw \
#	-ldfftw_mpi \
#	-DTYPE_SPECIFIED \
#	-DHAS_AUTO 
acml_on_generic_flags= \
	-lacml
# If your version of essl doesn't include PESSL,
#  delete -lpessl, -lblacs and -DHAS_AUTO.
essl_on_generic_flags= \
	-lessl \
	-lpessl \
	-lblacs \
	-DHAS_AUTO
# I'm not sure if all versions of MKL include the parallel
#  functionality, but if yours doesn't, remove -DHAS_AUTO
#  and whatever -ls no longer apply.
mkl_on_generic_flags= \
	-lmkl_cdft \
	-lmkl_cdft_core \
	-lm \
	-lmkl \
	-lmkl_core \
	-lmkl_def \
	-lguide \
	-DHAS_AUTO

# Antimony (only actually has fftw3, fftw2) 
fftw3_on_antimony_flags= \
	-L/Users/ik/packages/fftw-3.2alpha3/lib \
	-I/Users/ik/packages/fftw-3.2alpha3/include \
	-lfftw3_mpi \
	/Users/ik/packages/fftw-3.2alpha3/lib/libfftw3.a \
	-DHAS_AUTO
# Non-parallel version
#fftw3_on_antimony_flags= \
#	-L/usr/local/include \
#	-I/usr/local/lib \
#	-lfftw3
fftw2_on_antimony_flags= \
	-lfftw \
	-lfftw_mpi \
	-DHAS_AUTO \
	-L/Users/ik/packages/fftw-2.1.5/lib \
	-I/Users/ik/packages/fftw-2.1.5/include
acml_on_antimony_flags= \
	-lacml
essl_on_antimony_flags= -lessl \
	-lpessl \
	-lblacs \
	-DHAS_AUTO
mkl_on_antimony_flags= \
	-lmkl_cdft \
	-lmkl_cdft_core \
	-lm \
	-lmkl \
	-lmkl_core \
	-lmkl_def \
	-lguide \
	-DHAS_AUTO
	
# Ness - fftw3, fftw2, acml and mkl. (Oddly, since mkl is an intel lib)
fftw3_on_ness_flags= \
	-L/home/s07/fftw/fftw-3.1.2/lib \
	-I/home/s07/fftw/fftw-3.1.2/include \
	-lfftw3
fftw2_on_ness_flags= \
	-lfftw \
	-lfftw_mpi \
	-DHAS_AUTO \
	-L/home/s07/fftw/fftw-2.1.5/lib \
	-I/home/s07/fftw/fftw-2.1.5/include 
acml_on_ness_flags= \
	-lacml \
	-lpgftnrtl
mkl_on_ness_flags= \
	-L/export/intel/mkl/10.0.1.014/lib/em64t \
	-I/export/intel/mkl/10.0.1.014/include \
	-lmkl_cdft \
	-lmkl_cdft_core \
	-lm \
	-lmkl \
	-lmkl_core \
	-lmkl_def \
	-lguide 

# HECToR - fftw3, fftw2, acml
fftw3_on_hector_flags= \
	-lfftw3
fftw2_on_hector_flags= \
	-ldfftw \
	-DHAS_AUTO \
	-ldfftw_mpi \
	-DFFTW_TYPE_SPECIFIED
acml_on_hector_flags= \
	-lacml \
	-lpgftnrtl

# HPCx - fftw3, fftw2, (p)essl
fftw3_on_hpcx_flags= \
	-L/hpcx/usrlocal/packages/fftw/fftw3_64/lib \
	-I/hpcx/usrlocal/packages/fftw/fftw3_64/include \
	-lm \
	-lfftw3
fftw2_on_hpcx_flags= \
	-L/hpcx/usrlocal/packages/fftw/lib \
	-I/hpcx/usrlocal/packages/fftw/include \
	-lm \
	-ldfftw \
	-DHAS_AUTO \
	-ldfftw_mpi \
	-DFFTW_TYPE_SPECIFIED
essl_on_hpcx_flags= \
	-lessl \
	-lpesslsmp \
	-lblacssmp \
	-lm \
	-DHAS_AUTO

#BlueGene - fftw3, fftw2, essl
fftw3_on_bluegene_flags= \
	-lm \
	-lfftw3 \
	-L/home/b00/ikirker/packages/440d-O3-qhot/fftw-3.2alpha3/lib \
	-I/home/b00/ikirker/packages/440d-O3-qhot/fftw-3.2alpha3/include
fftw2_on_bluegene_flags= \
	-lm \
	-lfftw \
	-lfftw_mpi \
	-DHAS_AUTO \
	-L/home/b00/ikirker/packages/440d-O3-qhot/fftw-2.1.5/lib \
	-I/home/b00/ikirker/packages/440d-O3-qhot/fftw-2.1.5/include
essl_on_bluegene_flags= \
	-lm \
	-lesslbg  \
	-L/opt/ibmmath/lib \
	-I/opt/ibmmath/include \
	-L/opt/ibmcmp/xlf/bg/10.1/blrts_lib \
	-lxlf90 \
	-lxlfmath

#Eddie - MKL, FFTW3,2
fftw3_on_eddie_flags= \
        -L/usr/local.local/Cluster-Apps/fftw/intel/64/3.1.2/lib \
        -I/usr/local.local/Cluster-Apps/fftw/intel/64/3.1.2/include \
        -lfftw3 \
        -lm
fftw2_on_eddie_flags= \
        -I/usr/local/Cluster-Apps/fftw/intel/64/2.1.5/double/include \
		-L/exports/applications/apps/fftw/intel/64/2.1.5/double/lib \
        -lfftw \
        -lm
mkl_on_eddie_flags= \
        -L/exports/applications/apps/intel/mkl/10.0.1.014/lib/em64t \
        -I/exports/applications/apps/intel/mkl/10.0.1.014/include \
        -lmkl_cdft \
        -lmkl_cdft_core \
        -lm \
        -lmkl \
        -lmkl_core \
        -lmkl_def \
        -lguide

# MareNostrum in Barcelona - FFTW3,2, ESSL
fftw3_on_marenostrum= \
        -I/gpfs/apps/FFTW/3.1.1/64/include \
		-L/gpfs/apps/FFTW/3.1.1/64/lib/ \
		-lfftw3 \
		-lm
fftw2_on_marenostrum= \
		-I/gpfs/apps/FFTW/2.1.5/64/include \
		-L/gpfs/apps/FFTW/2.1.5/64/lib \
		-lfftw \
		-DHAS_AUTO
essl_on_marenostrum= \
        -lessl \
        -lpessl \
        -lblacs \
        -DHAS_AUTO

LIBFLAGS=$($(LIB)_on_$(SYSTEM)_flags) -DFFT_$(LIB)

# This is empty by default, but allows the specification of 
#  extra command-line arguments (e.g. library locations) at
#  the make command line.
EXTRAFLAGS=  

####################################################
#   TARGETS
####################################################

all: fft

fft: $(OBJ) Makefile
	$(MPICC) $(CFLAGS)  -o $@-$(LIB)  $(OBJ) $(LIBFLAGS) $(EXTRAFLAGS)

clean:
	-rm -f fft-* $(OBJ) *.oo
	
# Removes object files but not executables
sweep:
	-rm -f $(OBJ) *.oo

.c.o : $(HEADERS) $(SRC) Makefile
	$(MPICC) $(CFLAGS) -c $(@:.o=.c) $(LIBFLAGS) $(EXTRAFLAGS)
	
