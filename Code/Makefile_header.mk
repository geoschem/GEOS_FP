#EOC
# $Id: Makefile_header.mk,v 1.1.1.1 2009/08/04 14:52:05 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_header.mk (in Code subdirectory)
#
# !DESCRIPTION: This sub-makefile defines the variables which specify
# compilation options for the different compiler/platform combinations.  
# Also, the default makefile compilation rules are specified here.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# The following variables are exported to the main-level Makefile:
#
# Variable   Description
# --------   -----------
# F90        Contains the Fortran compilation commands
# FREEFORM   Contains the command to force F90 "free format" compilation
# LD         Contains the command to link to libraries & make executable
# LINK_NC    Specifies the command to link to the HDF libraries on this system
#
# FFLAGS, DIR_HDF, LINK_NC are local variables that are not returned 
# to the "outside world".
#
# !REVISION HISTORY: 
#  04 Aug 2009 - R. Yantosca - Initial version
#  01 Apr 2010 - R. Yantosca - Modified for netCDF-4.0.1 compiled for MPI
#                              as built with the NASA baselibs.
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Make ifort the default compiler
ifndef COMPILER
#COMPILER = mpif90
COMPILER = ifort
endif

###############################################################################
# Include directory for NETCDF library
# Modify this accordingly for your system!
INC_NC  = -I$(BL_INC_NETCDF) -I$(BL_INC_HDF5)
###############################################################################
# Library link commands for NETCDF library
# Modify this accordingly for your system!
LINK_NC = \
-L$(BL_LIB_NETCDF) -lnetcdf \
-L$(BL_LIB_HDF5) -lhdf5_hl \
-L$(BL_LIB_HDF5) -lhdf5 \
-L$(BL_LIB_ZLIB) -lz
###############################################################################

#==============================================================================
# MPIF90 compilation options (default)
#==============================================================================
ifeq ($(COMPILER),mpif90) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS = -cpp -w -noalign -convert big_endian -g -traceback 
else
FFLAGS = -cpp -w -O2 -auto -noalign -convert big_endian -openmp
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS += -traceback
endif

F90      = mpif90 $(FFLAGS) $(INC_NC)
LD       = mpif90 $(FFLAGS)
FREEFORM = -free

endif

#==============================================================================
# IFORT compilation options
#==============================================================================
ifeq ($(COMPILER),ifort) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS = -cpp -w -noalign -convert big_endian -g -traceback 
else
FFLAGS = -cpp -w -O2 -auto -noalign -convert big_endian -openmp
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS += -traceback
endif

F90      = ifort $(FFLAGS) $(INC_NC)
LD       = ifort $(FFLAGS)
FREEFORM = -free

endif

#==============================================================================
# Portland Group (PGI) compilation options
#==============================================================================
ifeq ($(COMPILER),pgi) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS = -byteswapio -Mpreprocess -fast -Bstatic
else
FFLAGS = -byteswapio -Mpreprocess -fast -mp -Mnosgimp -DHE4 -Bstatic
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

F90      = pgf90 $(FFLAGS) $(INC_NC)
LD       = pgf90 $(FFLAGS)
FREEFORM = -Mfree

endif

#==============================================================================
# SunStudio compilation options
#==============================================================================
ifeq ($(COMPILER),sun) 

# Default compilation options
# NOTE: -native builds in proper options for whichever chipset you have!
FFLAGS = -fpp -fast -stackvar -xfilebyteorder=big16:%all -native -DHE4

# Additional flags for parallel run
ifndef DEBUG
FFLAGS += -openmp=parallel
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

#---------------------------------------------------------------
# If your compiler is under the name "f90", use these lines!
#F90      = f90 $(FFLAGS) $(INC_NC)
#LD       = f90 $(FFLAGS)
#---------------------------------------------------------------
# If your compiler is under the name "sunf90", use these lines!
F90      = sunf90 $(FFLAGS) $(INC_NC)
LD       = sunf90 $(FFLAGS)
##---------------------------------------------------------------
FREEFORM = -free

endif

#==============================================================================
# IBM/XLF compilation options
# NOTE: someone who runs on IBM compiler should check this !!!
#==============================================================================
ifeq ($(COMPILER),xlf) 

# Default compilation options
FFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x80000000 -qfixed -qsuffix=cpp=f -q64

# Add optimization options
FFLAGS += -O3 -qarch=auto -qtune=auto -qcache=auto -qmaxmem=-1 -qstrict -DHE4

# Add more options for parallel run
ifndef DEBUG
FFLAGS += -qsmp=omp:opt -WF,-Dmultitask -qthreaded
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

F90      = xlf90_r $(FFLAGS) $(INC_NC)
LD       = xlf90_r $(FFLAGS)
FREEFORM = -qrealsize=8

endif

#==============================================================================
# Default compilation rules for *.f, *.f90, *.F, *.F90 and *.c files
#==============================================================================
.SUFFIXES: .f .F .f90 .F90 .c
.f.o:                   ; $(F90) -c $*.f
.F.o:                   ; $(F90) -c $*.F
.f90.o:                 ; $(F90) -c $(FREEFORM) $*.f90 
.F90.o:                 ; $(F90) -c $(FREEFORM) $*.F90 

#==============================================================================
# Export global variables so that the main Makefile will see these
#==============================================================================
export F90
export FREEFORM
export LD
export LINK_NC
#EOC

#==============================================================================
# Print variables for testing/debugging purposes (uncomment if necessary)
#==============================================================================
#headerinfo:
#	@echo '####### in Makefile_header.mk ########' 
#	@echo "compiler: $(COMPILER)"
#	@echo "debug   : $(DEBUG)"
#	@echo "bounds  : $(BOUNDS)"
#	@echo "f90     : $(F90)"
#	@echo "ld      : $(LD)"
#	@echo "link_nc : $(LINK_NC)"
#	@echo "cc      : $(CC)"

