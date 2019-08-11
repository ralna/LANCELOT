#
#  Build LANCELOT A and its SIF file decoder.
#
#  version for Unix systems using gfortran

#  Nick Gould, August 8th, 2019, for CGT Productions.

#  Basic system commands

RM=rm -f
CP=cp
MV=mv

#  Temporary directory for compilation
#  ( do NOT use TMP = . which would cause file overwriting)

TMP=/tmp

#  Default sizing program

SIZE=large

#  Default precision

PRECISION=double

#  Compiler flags (machine dependent)

FORTRAN=gfortran
FFLAGS= -O3 -ffixed-form

#  The compilation

install:
	$(FORTRAN) -o $(SIZE) $(SIZE).f
	$(FORTRAN) -o $(PRECISION) $(PRECISION).f
	cd ./sifdec ; $(MAKE) sifdec_$(PRECISION) \
                        PRECISION=$(PRECISION) SIZE=$(SIZE) \
                        FORTRAN="$(FORTRAN)" FFLAGS="$(FFLAGS)" \
                        RM="$(RM)" CP="$(CP)" MV="$(MV)" TMP="$(TMP)"
	cd ./optimizers ; $(MAKE) lancelot_$(PRECISION) \
                        PRECISION=$(PRECISION) SIZE=$(SIZE) \
                        FORTRAN="$(FORTRAN)" FFLAGS="$(FFLAGS)" \
                        RM="$(RM)" CP="$(CP)" TMP="$(TMP)"

# HELP!

help:
	@ echo "Possible choices are:"
	@ echo " - help: this information"
	@ echo " - lancelot: install LANCELOT A and its SIF file decoder"
	@ echo " -  options (default)"
	@ echo " -   FORTRAN=your fortran 77 compiler (gfortran)"
	@ echo " -   FFLAGS=your fortran compiler flags (-O3 -ffixed-form)"
	@ echo " -   SIZE=small medium large huge (large)"
	@ echo " -   PRECISION=single double (double)"
	@ echo " -   RM=rm command with force (/bin/rm -f)"
	@ echo " -   CP=cp command (/bin/cp)"
	@ echo " -   MV=mv command (/bin/mv)"
	@ echo " -   TMP=temporary workspace directory, not ./ (/tmp)"

#  Clean up afterwards

clean:
	cd ./sifdec ; $(MAKE) clean PRECISION=$(PRECISION)
	cd ./optimizers ; $(MAKE) clean PRECISION=$(PRECISION)
	$(RM) $(PRECISION) $(SIZE)
