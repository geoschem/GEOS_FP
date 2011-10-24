#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in doc subdirectory)
#
# !DESCRIPTION: Makefile for building the documentation (in PDF and 
#  PostScript formats) for the GEOS-5 data processing code, Makefiles,
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
#  28 Jul 2009 - R. Yantosca - Initial version
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

doc: srcdoc makedoc perldoc

clean:
	rm -f *.tex *.ps *.pdf

#==============================================================================
# Build documentation for source code
#==============================================================================

# Source files
SRC1 =                              \
./intro.txt                         \
$(CODEDIR)/MerraDriver.F90          \
$(CODEDIR)/MerraInputsModule.F90    \
$(CODEDIR)/MerraBinaryModule.F90    \
$(CODEDIR)/MerraRegridModule.F90    \
$(CODEDIR)/MerraA3Module.F90        \
$(CODEDIR)/MerraA6Module.F90        \
$(CODEDIR)/MerraCnModule.F90        \
$(CODEDIR)/MerraI6Module.F90        \
$(CODEDIR)/He4ErrorModule.F90       \
$(CODEDIR)/He4GridModule.F90        \
$(CODEDIR)/He4IncludeModule.F90     \
$(CODEDIR)/CharpakModule.F

# LaTeX etc. files
TEX1 = MerraDriver.tex
DVI1 = MerraDriver.dvi
PDF1 = MerraDriver.pdf
PS1  = MerraDriver.ps

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