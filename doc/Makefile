#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in doc subdirectory)
#
# !DESCRIPTION: Makefile for building the documentation (in PDF and 
#  PostScript formats) for the GEOS-5.7.x data processing code, Makefiles,
#  and Perl scripts.
#\\
#\\
# !REMARKS:
# To build the documentation, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY: 
#  12 Jan 2012 - R. Yantosca - Initial version
#  12 Jan 2012 - R. Yantosca - Split source code doc from netCDF code doc
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL   = /bin/sh
BINDIR  = ../bin
CODEDIR = ../Code
MAINDIR = ../
PERLDIR = ../perl

#==============================================================================
# Default targets
#==============================================================================

.PHONY: doc clean

doc: srcdoc ncdoc makedoc perldoc

clean:
	rm -f *.tex *.ps *.pdf

#==============================================================================
# Build documentation for source code
#==============================================================================

# Source files
SRC1 =                                 \
./src.intro                            \
$(CODEDIR)/CreateTemplateFile.F90      \
$(CODEDIR)/CharpakModule.F	       \
$(CODEDIR)/Geos57A1Module.F90	       \
$(CODEDIR)/Geos57A3CldModule.F90       \
$(CODEDIR)/Geos57A3DynModule.F90       \
$(CODEDIR)/Geos57A3MstCModule.F90      \
$(CODEDIR)/Geos57A3MstEModule.F90      \
$(CODEDIR)/Geos57CnModule.F90	       \
$(CODEDIR)/Geos57Driver.F90	       \
$(CODEDIR)/Geos57I3Module.F90	       \
$(CODEDIR)/Geos57InputsModule.F90      \
$(CODEDIR)/Geos57RegridModule.F90      \
$(CODEDIR)/Geos57UtilityModule.F90     \
$(CODEDIR)/Geos57Driver.F90

# LaTeX etc. files
TEX1 = Geos57Driver.tex
DVI1 = Geos57Driver.dvi
PDF1 = Geos57Driver.pdf
PS1  = Geos57Driver.ps

srcdoc:
	rm -f $(TEX1)
	protex -sf $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	dvips $(DVI1) -o $(PS1)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for netCDF library code
#==============================================================================

# Source files
SRC4 =                                 \
./ncdf.intro                           \
$(CODEDIR)/m_do_err_out.F90	       \
$(CODEDIR)/m_netcdf_io_checks.F90      \
$(CODEDIR)/m_netcdf_io_close.F90       \
$(CODEDIR)/m_netcdf_io_create.F90      \
$(CODEDIR)/m_netcdf_io_define.F90      \
$(CODEDIR)/m_netcdf_io_get_dimlen.F90  \
$(CODEDIR)/m_netcdf_io_handle_err.F90  \
$(CODEDIR)/m_netcdf_io_open.F90	       \
$(CODEDIR)/m_netcdf_io_read.F90	       \
$(CODEDIR)/m_netcdf_io_write.F90

# LaTeX etc. files
TEX4 = NcdfUtilities.tex
DVI4 = NcdfUtilities.dvi
PDF4 = NcdfUtilities.pdf
PS4  = NcdfUtilities.ps

ncdoc:
	rm -f $(TEX4)
	protex -sf $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	dvips $(DVI4) -o $(PS4)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for Makefiles
#==============================================================================

# Source files 
SRC2 =                        \
./make.intro                  \
$(MAINDIR)/Makefile           \
$(CODEDIR)/Makefile           \
$(CODEDIR)/Makefile_header.mk \
./Makefile

# LaTeX etc. files
TEX2 = Geos57Driver_Makefiles.tex
DVI2 = Geos57Driver_Makefiles.dvi
PDF2 = Geos57Driver_Makefiles.pdf
PS2  = Geos57Driver_Makefiles.ps

makedoc:
	rm -f $(TEX2)
	protex -fS $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	dvips $(DVI2) -o $(PS2)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for Perl scripts
#==============================================================================

# Source files 
SRC3 =                    \
./perl.intro              \
$(PERLDIR)/doGeos57       \
$(PERLDIR)/doGeos57.input \
$(PERLDIR)/getGeos57      \
$(PERLDIR)/checkGeos57    \
$(PERLDIR)/sleepGeos57    \
$(PERLDIR)/schedGeos57    \
$(PERLDIR)/cleanJobs      \
$(PERLDIR)/cleanLogs      \
$(PERLDIR)/delGeos57      \
$(PERLDIR)/Dates.pm

# LaTeX etc. files
TEX3 = Geos57Driver_Perl.tex
DVI3 = Geos57Driver_Perl.dvi
PDF3 = Geos57Driver_Perl.pdf
PS3  = Geos57Driver_Perl.ps

perldoc:
	rm -f $(TEX3)
	protex -sfS $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	dvips $(DVI3) -o $(PS3)
	rm -f *.aux *.dvi *.log *.toc


#==============================================================================
#  Help Screen
#==============================================================================

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'doc      Builds documentation (*.ps, *.pdf) for both source code and makefiles (default)'
	@echo 'srcdoc   Builds documentation (*.ps, *.pdf) only for source code files'
	@echo 'makedoc  Builds documentation (*.ps, *.pdf) only for makefiles'
	@echo 'clean    Removes all files (*.tex, *.ps, *.pdf)'
	@echo 'help     Displays this help screen'

#EOC