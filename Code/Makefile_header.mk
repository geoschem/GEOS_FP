#EOC
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
#  To build the programs, call "make" with the following syntax:
#
#    make TARGET [ OPTIONAL-FLAGS ]
#
#  To display a complete list of options, type "make help".
#
#  The following variables are accepted either as command-line options,
#  or may be defined in your ~/.cshrc or ~/.bashrc file:
#                                                                             .
#  Variable     Description
#  ----------   -----------
#  BIN_NETCDF   Specifies the path for netCDF etc. executables
#  INC_NETCDF   Specifies the path for netCDF etc. include files & modules
#  LIB_NETCDF   Specifies the path for netCDF etc. libraries
#                                                                             .
#  The following variables are exported to the main-level Makefile:
#
#  Variable   Description
#  --------   -----------
#  F90        Contains the Fortran compilation commands
#  FREEFORM   Contains the command to force F90 "free format" compilation
#  LD         Contains the command to link to libraries & make executable
#  LINK_NC    Specifies the command to link to the HDF libraries on this system
#
#  FFLAGS, DIR_HDF, LINK_NC are local variables that are not returned 
#  to the "outside world".
#
# !REVISION HISTORY: 
#  24 Oct 2011 - R. Yantosca - Initial version, based on GEOS-5
#  15 Feb 2012 - R. Yantosca - Now compile IFORT w/ -mcmodel=medium -i-dynamic
#  11 May 2012 - R. Yantosca - Now attempt to use nf-config, nc-config to
#                              obtain the library linking sequence.  This will
#                              make the Makefile much more portable.
#  11 May 2012 - R. Yantosca - Now use INC_NETCDF, BIN_NETCDF, LIB_NETCDF
#                              env variables to specify directory paths 
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Make ifort the default compiler
ifndef COMPILER
COMPILER := ifort
endif

###############################################################################
###                                                                         ###
###  Set linker commands for local and external libraries (incl. netCDF)    ###
###                                                                         ###
###############################################################################

ifdef NETCDF_INCLUDE
 NC_INC_CMD           := -I$(NETCDF_INCLUDE)
else
 NC_INC_CMD           := -I$(INC_NETCDF)
endif

ifdef NETCDF_FORTRAN_INCLUDE
 NC_INC_CMD           += -I$(NETCDF_FORTRAN_INCLUDE)
endif

# Get the version number (e.g. "4130"=netCDF 4.1.3; "4200"=netCDF 4.2, etc.)
NC_VERSION           :=$(shell nc-config --version)
NC_VERSION           :=$(shell echo "$(NC_VERSION)" | sed 's|netCDF ||g')
NC_VERSION           :=$(shell echo "$(NC_VERSION)" | sed 's|\.||g')
NC_VERSION_LEN       :=$(shell perl -e "print length $(NC_VERSION)")
ifeq ($(NC_VERSION_LEN),3)
 NC_VERSION          :=$(NC_VERSION)0
endif
ifeq ($(NC_VERSION_LEN),2) 
 NC_VERSION          :=$(NC_VERSION)00
endif

# Test if we have at least netCDF 4.2.0.0
AT_LEAST_NC_4200     :=$(shell perl -e "print ($(NC_VERSION) ge 4200)")

ifeq ($(AT_LEAST_NC_4200),1) 

  #-------------------------------------------------------------------------
  # netCDF 4.2 and higher:
  # Use "nf-config --flibs" and "nc-config --libs"
  # Test if a separate netcdf-fortran path is specified
  #-------------------------------------------------------------------------
  NC_LINK_CMD        := $(shell nf-config --flibs)
  NC_LINK_CMD        += $(shell nc-config --libs)

else

  #-----------------------------------------------------------------------
  # Prior to netCDF 4.2:
  # Use "nc-config --flibs"
  #-----------------------------------------------------------------------
  NC_LINK_CMD        := $(shell nc-config --flibs)

endif

#=============================================================================
#%%%%% FIX FOR USE WITH THE GEOS-Chem-Libraries (bmy, 1/13/15)
#%%%%% 
#%%%%% If your GEOS-Chem-Libraries netCDF/HDF5 package was built in one 
#%%%%% directory and then moved somewhere else, then nf-config and nc-config 
#%%%%% may not return the proper link directory path.  
#%%%%% 
#%%%%% To avoid this error, we shall test if the $GC_LIB environment variable 
#%%%%% contains the text "GEOS-Chem-Libraries".  (Recall that $GC_LIB is 
#%%%%% defined in either your .bashrc or .cshrc file depending on which Unix 
#%%%%% shell you use.)  If we find the text "GEOS-Chem-Libraries" in $GC_LIB, 
#%%%%% then we shall override the library path returned by nf-config and 
#%%%%% nc-config with the path specified by $GC_LIB.  This will ensure that 
#%%%%% we point to the location where the GEOS-Chem-Libraries are installed.
#%%%%%
#%%%%% NOTE: This fix should work for most users.  If it does not work, then
#%%%%% contact the GEOS-Chem Support Team (geos-chem-support@as.harvard.edu).
#%%%%%
REGEXP               :="GEOS-Chem-Libraries"
ifeq ($(shell [[ "$(LIB_NETCDF)" =~ $(REGEXP) ]] && echo true),true)
  NC_LINK_CMD        := $(filter -l%,$(NC_LINK_CMD))
  NC_LINK_CMD        :=-L$(LIB_NETCDF) $(NC_LINK_CMD)
endif
#=============================================================================

# For backwards compatibility
LINK_NC :=$(NC_LINK_CMD)
INC_NC  :=$(NC_INC_CMD)

###############################################################################
###                                                                         ###
###  Test if the netCDF library was built with compression enabled          ###
###                                                                         ###
###############################################################################

# Test if the "nf_def_var_deflate" function is defined in netcdf.inc
# Look for netcdf.inc where the netCDF-Fortran library is located
ifdef NETCDF_FORTRAN_INCLUDE
  GREP :=$(strip $(shell grep nf_def_var_deflate $(NETCDF_FORTRAN_INCLUDE)/netcdf.inc))
else
  GREP :=$(strip $(shell grep nf_def_var_deflate $(NETCDF_INCLUDE)/netcdf.inc))
endif

# Look for the second word of the combined search results
WORD                 :=$(word 2,"$(GREP)")

# If it matches "nf_def_var_deflate", then define Cpp flag NC_HAS_COMPRESSION 
ifeq ($(WORD),nf_def_var_deflate)
  USER_DEFS          += -DNC_HAS_COMPRESSION
endif

#==============================================================================
# MPIF90 compilation options 
#==============================================================================
ifeq ($(COMPILER),mpif90) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS   := -cpp -w -O0 -auto -noalign -mcmodel=medium -shared-intel -g -traceback
else
FFLAGS   := -cpp -w -O2 -auto -noalign -mcmodel=medium -shared-intel -openmp
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS   += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS   += -traceback
endif

# Add any extra definitions
FFLAGS   += $(USER_DEFS)

INCLUDE  := -module $(MOD) -I$(MOD) $(INC_NC)
F90      := mpif90 $(FFLAGS) $(INCLUDE)
LD       := mpif90 $(FFLAGS) $(INCLUDE)
FREEFORM := -free

endif

#==============================================================================
# IFORT compilation options (default)
#==============================================================================
ifeq ($(COMPILER),ifort) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS   := -cpp -w -O0 -auto -noalign -mcmodel=medium -shared-intel -g -traceback
else
FFLAGS   := -cpp -w -O2 -auto -noalign -mcmodel=medium -shared-intel -openmp
endif

# Add flag to denote if we are using the sample data (wh
ifdef USE_SAMPLE_DATA
FFLAGS   += -DUSE_SAMPLE_DATA
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS   += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS   += -traceback
endif

# Add any extra definitions
FFLAGS   += $(USER_DEFS)

INCLUDE  := -module $(MOD) -I$(MOD) $(INC_NC)
F90      := ifort $(FFLAGS) $(INCLUDE)
LD       := ifort $(FFLAGS) $(INCLUDE)
FREEFORM := -free

endif

#==============================================================================
# Portland Group (PGI) compilation options
#==============================================================================
ifeq ($(COMPILER),pgi) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS   := -byteswapio -Mpreprocess -fast -Bstatic
else
FFLAGS   := -byteswapio -Mpreprocess -fast -mp -Mnosgimp -DHE4 -Bstatic
endif

# Add flag to denote if we are using the sample data (wh
ifdef USE_SAMPLE_DATA
FFLAGS   += -DUSE_SAMPLE_DATA
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS   += -C
endif

# Add any extra definitions
FFLAGS   += $(USER_DEFS)

INCLUDE  := -module $(MOD) -I$(MOD) $(INC_NC)
F90      := pgf90 $(FFLAGS) $(INCLUDE)
LD       := pgf90 $(FFLAGS) $(INCLUDE)
FREEFORM := -Mfree

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

# Add flag to denote if we are using the sample data (wh
ifdef USE_SAMPLE_DATA
FFLAGS   += -DUSE_SAMPLE_DATA
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

# Add any extra definitions
FFLAGS   += $(USER_DEFS)

# Include options
INCLUDE  := -module $(MOD) -I$(MOD) $(INC_NC)

#---------------------------------------------------------------
# If your compiler is under the name "f90", use these lines!
#F90      = f90 $(FFLAGS) $(INCLUDE)
#LD       = f90 $(FFLAGS) $(INCLUDE)
#---------------------------------------------------------------
# If your compiler is under the name "sunf90", use these lines!
F90      = sunf90 $(FFLAGS) $(INCLUDE)
LD       = sunf90 $(FFLAGS) $(INCLUDE)
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

# Add any extra definitions
FFLAGS   += $(USER_DEFS)

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
#	@echo "NETCDF_INCLUDE : $(NETCDF_INCLUDE)"

