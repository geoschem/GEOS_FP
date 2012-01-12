#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in Code subdirectory)
#
# !DESCRIPTION: Makefile for the GEOS-5.7.x data processing code.
#\\
#\\
# !REMARKS:
# To build the program, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# BIN        Specifies the directory where the executable will be sent
# DOC        Specifies the directory where documentation will be built
# EXE        Specifies the name of the executable file
# F90        Contains the default Fortran compilation commands
# LD         Contains the command to link to libraries & make executable
# LIB        Specifies the directory where the libNcUtils.a files will be sent
# LINK_NC    Contains the comand to link to the netCDF library
# MOD        Specifies the directory where module files (*.mod) will be sent
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
#                                                                             .
# F90, LD, LINK_NC are included from "Makefile_header.mk"
#
# !REVISION HISTORY: 
#  24 Oct 2011 - R. Yantosca - Initial version, based on GEOS-5
#  11 Jan 2012 - R. Yantosca - Now compile Geos57A3*Module.F90 files
#  12 Jan 2012 - R. Yantosca - Split Geos57A3MstModule.F90 into two modules
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL   = /bin/sh
BIN     = ../bin
DOC     = ../doc
EXE     = Geos57Driver.x
LIB     = ../lib
MOD     = ../mod

# Include header file.  This returns variables F90, LD, LINK_NC
# as well as the default Makefile compilation rules for source code files.
include ./Makefile_header.mk

#==============================================================================
# Filename section
#
# NOTES:
# (1) The function $(wildcard m_*.F90) returns all files matching m_*.F90, etc.
# (2) The syntax $(SRC:.F90=.o) is a shorthand for the "patsubst" function. 
#==============================================================================

#----------------------
# Objects
#----------------------

# Source code & object files
SRC  := $(wildcard *.F) $(wildcard *.F90)            # Get these files
SRC  := $(filter-out CreateTemplateFile.F90,$(SRC))  # but skip these files
TMP  := $(SRC:.F=.o)
OBJ  := $(TMP:.F90=.o)

# Object file for executable
OBJx := $(EXE:.x=.o)

#-----------------------
# Dependency listing
#-----------------------
m_do_err_out.o           : m_do_err_out.F90

m_netcdf_io_close.o      : m_netcdf_io_close.F90        \
                           m_do_err_out.o

m_netcdf_io_create.o     : m_netcdf_io_create.F90       \
                           m_do_err_out.o

m_netcdf_io_define.o     : m_netcdf_io_define.F90       \
                           m_do_err_out.o

m_netcdf_io_get_dimlen.o : m_netcdf_io_get_dimlen.F90   \
                           m_do_err_out.o

m_netcdf_io_handle_err.o : m_netcdf_io_handle_err.F90   \
                           m_do_err_out.o

m_netcdf_io_open.o       : m_netcdf_io_open.F90         \
                           m_do_err_out.o

m_netcdf_io_read.o       : m_netcdf_io_read.F90         \
                           m_do_err_out.o

m_netcdf_io_write.o      : m_netcdf_io_write.F90        \
                           m_do_err_out.o

CharpakModule.o          : CharpakModule.F

Geos57InputsModule.o     : Geos57InputsModule.F90       \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57RegridModule.o     : Geos57RegridModule.F90       \
                           Geos57InputsModule.o

Geos57UtilityModule.o    : Geos57UtilityModule.F90      \
                           CharpakModule.o

Geos57A1Module.o         : Geos57A1Module.F90           \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57A3CldModule.o      : Geos57A3CldModule.F90        \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57A3DynModule.o      : Geos57A3DynModule.F90        \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57A3MstCModule.o     : Geos57A3MstCModule.F90       \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57A3MstEModule.o     : Geos57A3MstEModule.F90       \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57CnModule.o         : Geos57CnModule.F90           \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57I3Module.o         : Geos57I3Module.F90           \
                           CharpakModule.o              \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57UtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

Geos57Driver.o           : Geos57Driver.F90             \
                           Geos57InputsModule.o         \
                           Geos57RegridModule.o         \
                           Geos57A1Module.o             \
                           Geos57A3CldModule.o          \
                           Geos57A3DynModule.o          \
                           Geos57A3MstCModule.o         \
                           Geos57A3MstEModule.o         \
                           Geos57CnModule.o             \
                           Geos57I3Module.o

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: all doc docclean clean realclean help test driver

all: lib exe

lib: $(OBJ)
	$(AR) crs libGeos57.a $(OBJ)
	mv libGeos57.a $(LIB)

exe: $(OBJx)
	$(LD) $(OBJx) -L$(LIB) -lGeos57 $(LINK_NC) -o $(EXE)
	cp -f $(EXE) $(BIN)

doc:
	$(MAKE) -C $(DOC) doc

docclean:
	$(MAKE) -C $(DOC) clean

clean:
	rm -f *.o *.nc *.x *.a *.mod *.nc

realclean:
	$(MAKE) clean
	$(MAKE) docclean
	rm -f $(BIN)/*.x
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod

driver: all

#==============================================================================
# Help Screen
#==============================================================================

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'lib        Compiles code and creates the *.mod, *.o, and *.a files'
	@echo 'exe        Builds the executable file'
	@echo 'all        Synonym for "lib exe" (default target)'

	@echo 'check      Compiles code and runs the executable using GMAO sample data (for debugging)'
	@echo 'clean      Removes *.o *.mod and executable files in Code subdirectory'
	@echo 'realclean  Makes clean, then removes contents of bin, lib, doc subdirectories'
	@echo 'doc        Builds documentation files (*.ps, *.pdf) in doc subdirectory'
	@echo 'docclean   Removes documentation files from doc subdirectory'
	@echo 'help       Displays this help screen'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'TRACEBACK=yes  Compiles with -traceback flag (ifort only)'
#
#EOC