#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (run-directory)
#
# !DESCRIPTION: Makefile for the GEOS-5.7.x data processing code.  This is a
#  "router" makefile, which calls to the Makefile in the Code directory.
#\\
#\\
# !REMARKS:
# To build the program, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY: 
# 12 Jan 2012 - R. Yantosca - Add "check" target to build & run the data
#                             processing code using the GMAO sample data
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL = /bin/sh
DIR   = ../Code
EXE0   = Geos57Driver0.x
EXE1   = Geos57Driver1.x
EXE2   = Geos57Driver2.x

#==============================================================================
# Makefile targets (for now, build w/ the sample data)
#==============================================================================

.PHONY: all norun debug driver lib check clean realclean doc docclean help

all: 
	$(MAKE) -C $(DIR) all

lib:
	$(MAKE) -C $(DIR) lib

check: 
	$(MAKE) -C $(DIR) all USE_SAMPLE_DATA=yes
	$(EXE) < 20110815

debug:
	$(MAKE) -C $(DIR) all
	$(EXE0) < 20110901 > log.0 &
	$(EXE1) < 20110901 > log.1 &
	$(EXE2) < 20110901 > log.2 &
	$(EXE3) < 20110901 > log.3 &
	$(EXE4) < 20110901 > log.4 &
	$(EXE5) < 20110901 > log.5 &

clean:
	$(MAKE) -C $(DIR) clean

realclean:
	$(MAKE) -C $(DIR) realclean
	rm -f GEOS572*.nc         

doc:
	$(MAKE) -C $(DIR) doc

docclean:
	$(MAKE) -C $(DIR) docclean

help:
	$(MAKE) -C $(DIR) help

driver: all

template:
	$(MAKE) -C $(DIR) -f Makefile_Template all

#EOC
