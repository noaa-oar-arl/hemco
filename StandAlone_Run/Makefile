#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/run subdirectory)
#
# !DESCRIPTION: Makefile for building an executable file for the standalone
#  HEMCO simulation.  It also creates a customized version of the HEMCO
#  Configuration file (HEMCO_Config.rc) and runs a short test simulation.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#
#   make -jN TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY: 
#  16 Jul 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT      :=../..
MOD       :=$(ROOT)/mod
LIB       :=$(ROOT)/lib
BIN       :=$(ROOT)/bin

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# Executable file
OBJx      :=hemco_standalone.o
RUN       :=./run
EXE       :=$(OBJx:.o=.x)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete list!                ###
###                                                                         ###
###############################################################################

.PHONY: all check clean debug

all: check

check:
	@$(MAKE) exe

clean:
	rm -f *.o *.mod *.a *.x *.nc *.log 
	rm -f $(RUN)/*.x

exe: $(OBJx)
	@echo "LINK: $(LINK_HCO)"
	$(LD) $(OBJx) $(LINK_HCO) -o $(EXE)
	cp -f $(EXE) $(RUN)

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "RUN     : $(RUN)"
	@echo "F90     : $(F90)"
	@echo "LINK    : $(LINK_HCO)"
#EOC
