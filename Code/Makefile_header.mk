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

###############################################################################
###                                                                         ###
###  Set the default compiler, based on environment variable FC             ###
###                                                                         ###
###############################################################################

# %%%%% Test if Intel Fortran Compiler is selected %%%%%
REGEXP               :=(^[Ii][Ff][Oo][Rr][Tt])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)

  # If we are building GCHP, then set the compile command to "mpifort",
  # which invokes the MPI magic.  Otherwise set it to $(FC). (bmy, 10/17/16)
  COMPILER_FAMILY    :=Intel
  USER_DEFS          += -DLINUX_IFORT
endif

# %%%%% Test if GNU Fortran Compiler is selected %%%%%
REGEXP               :=(^[Gg][Ff][Oo][Rr][Tt][Rr][Aa][Nn])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
  COMPILER_FAMILY    :=GNU
  USER_DEFS          += -DLINUX_GFORTRAN
endif

# Compiler command
COMPILE_CMD        :=$(FC)

# %%%%% ERROR CHECK!  Make sure our compiler selection is valid! %%%%%
REGEXP               :=((-DLINUX_)?IFORT|GFORTRAN)
ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
  $(error $(ERR_CMPLR))
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

###############################################################################
###                                                                         ###
###  Define settings for the GNU FORTRAN COMPILER (aka gfortran)            ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER_FAMILY),GNU) 

  # Get the GNU Fortran version
  GNU_VERSIONTEXT    :=$(shell $(FC) --version)
  GNU_VERSION        :=$(word 4, $(GNU_VERSIONTEXT))
  GNU_VERSION        :=$(subst .,,$(GNU_VERSION))
  NEWER_THAN_447     :=$(shell perl -e "print ($(GNU_VERSION) gt 447)")

  # Base set of compiler flags
  FFLAGS             :=-cpp -w -std=legacy -fautomatic -fno-align-commons
  ifeq ($(IS_HPC),1)
    FFLAGS             += -fconvert=native
  else
    FFLAGS             += -fconvert=big-endian
  endif
  FFLAGS             += -fno-range-check

  # OPTIONAL: Add the GNU Fortran -march option, which compiles for a
  # specific computer architecture.  This may cause issues on some types
  # of CPUs (e.g. Intel), so we have left this as an optional argument.
  ifdef M_ARCH
    FFLAGS           += -march=$(M_ARCH)
  endif

  # Default optimization level for all routines (-O3)
  ifndef OPT
    # Options of interest
    #  -limf                Intel math libraries - machine must have them
    #  -O3                  Highest safe optimization level
    #  -march=native        Make the binary machine-specific. If in doubt, 
    #                        use a specific architecture, eg...
    #  -march=corei7-avx    Binary uses optimizations for 
    #                        Intel Sandy-Bridge Xeon (e.g. E5-2680)
    #  -mfpmath=sse         Use SSE extensions
    #  -funroll-loops       Enable loop unrolling
    OPT              := -O3 -funroll-loops
    #OPT              := -O3 -march=corei7-avx -mfpmath=sse -funroll-loops
  endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    #-fcheck=all would be more comprehensive but would force bounds checking
    FFLAGS           += -g -gdwarf-2 -gstrict-dwarf -O0
    FFLAGS           += -Wall -Wextra -Wconversion
    FFLAGS           += -Warray-temporaries -fcheck-array-temporaries
    TRACEBACK        := yes
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           += $(OPT)
  endif

  # Prevent any optimizations that would change numerical results
  #GFORTRAN_BAD#FFLAGS             += -fp-model source

  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fopenmp
  endif

  # Get Operating System (Linux = Linux; Darwin = MacOSX)
  ifndef UNAME
    UNAME            :=$(shell uname)
  endif

  # OSX compilation options
  ifeq ($(UNAME),Darwin)
    # This has not yet been tested
    $(error $(ERR_OSCOMP))
  #  FFLAGS           += -Wl,-stack_size,0x2cb410000  # 12 GB of stack space
  #  ifdef DEBUG
  #    FFLAGS         += -g0 -debug -save-temps -fpic -Wl,-no_pie
  #  endif
  endif

  # Add options for medium memory model.  This is to prevent G-C from 
  # running out of memory at hi-res, especially when using netCDF I/O.
  ifneq ($(UNAME),Darwin)
    #GFORTRAN_BAD#FFLAGS           += -mcmodel=medium -shared-intel
    FFLAGS           += -mcmodel=medium
  endif

  # Turn on checking for floating-point exceptions
  # These are approximately equivalent to -fpe0 -ftrapuv in IFORT
  # NOTE: GNU Fortran 4.4.7 does not allow for -finit-real-snan, so
  # we will only add this flag for versions newer than 4.4.7
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(FPE)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -ffpe-trap=invalid,zero,overflow
    ifeq ($(NEWER_THAN_447),1)
      FFLAGS           += -finit-real=snan
    endif
  endif
  ifeq ($(shell [[ "$(FPEX)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -ffpe-trap=invalid,zero,overflow
    ifeq ($(NEWER_THAN_447),1)
      FFLAGS           += -finit-real=snan
    endif
  endif

  # Add option for "array out of bounds" checking
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(BOUNDS)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fbounds-check
  endif

  # Also add traceback option
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TRACEBACK)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fbacktrace
    ifndef DEBUG
       FFLAGS += -g
    endif
  endif

  # Compile for use with the GNU profiler (gprof), if necessary
  ifeq ($(IS_GPROF),1) 
    FFLAGS           += -pg
  endif

  # Append the user options in USER_DEFS to FFLAGS
  FFLAGS             += $(USER_DEFS)

  # Include options (i.e. for finding *.h, *.mod files)
  INCLUDE :=-J$(MOD) $(NC_INC_CMD)

  # Set the standard compiler variables
  CC                 :=
  F90                :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE)
  F90ISO             :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE_ISO)
  LD                 :=$(COMPILE_CMD) $(FFLAGS)
  FREEFORM           := -ffree-form -ffree-line-length-none
  R8                 := -fdefault-real-8 -fdefault-double-8

endif

###############################################################################
###                                                                         ###
###  Define settings for the INTEL FORTRAN COMPILER (aka ifort)             ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER_FAMILY),Intel) 

  # Base set of compiler flags
  FFLAGS             :=-cpp -w -auto -noalign -convert big_endian

  # Default optimization level for all routines (-O2)
  ifndef OPT
    OPT              := -O2
  endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -g -O0 -check arg_temp_created -debug all
    TRACEBACK        := yes
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           += $(OPT) -vec-report0
  endif

  # Prevent any optimizations that would change numerical results
  FFLAGS             += -fp-model source

  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -openmp
  endif

  # Get Operating System (Linux = Linux; Darwin = MacOSX)
  ifndef UNAME
    UNAME            :=$(shell uname)
  endif

  # OSX compilation options
  ifeq ($(UNAME),Darwin)
    FFLAGS           += -Wl,-stack_size,0x2cb410000  # 12 GB of stack space
    ifdef DEBUG
      FFLAGS         += -g0 -debug -save-temps -fpic -Wl,-no_pie
    endif
  endif

  # Add options for medium memory model.  This is to prevent G-C from 
  # running out of memory at hi-res, especially when using netCDF I/O.
  ifneq ($(UNAME),Darwin)
    FFLAGS           += -mcmodel=medium -shared-intel
  endif

  # Turn on checking for floating-point exceptions
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(FPE)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fpe0 -ftrapuv
  endif
  ifeq ($(shell [[ "$(FPEX)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fpe0 -ftrapuv
  endif

  # Add special IFORT optimization commands
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(IPO)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -ipo -static
  endif

  # Add option for "array out of bounds" checking
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(BOUNDS)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -check bounds
  endif

  # Also add traceback option
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TRACEBACK)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -traceback
  endif

  # Compile for use with the GNU profiler (gprof), if necessary
  ifeq ($(IS_GPROF),1) 
    FFLAGS           += -p
  endif

  # Append the user options in USER_DEFS to FFLAGS
  FFLAGS             += $(USER_DEFS)

  # Include options (i.e. for finding *.h, *.mod files)
  INCLUDE            :=-module $(MOD) $(NC_INC_CMD)

  # Set the standard compiler variables
  CC                 :=
  F90                :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE)
  F90ISO             :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE_ISO)
  LD                 :=$(COMPILE_CMD) $(FFLAGS)
  FREEFORM           := -free
  #ifneq ($(shell [[ "$(HPC)" =~ $(REGEXP) ]] && echo true),true)
  #ifneq ($(HPC),yes)
    R8               := -r8
  #endif

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

