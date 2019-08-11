#
#  Build everything for LANCELOT A
#
#  version for Unix systems using gfortran

#  Nick Gould, August 8th, 2019, for CGT Productions.

#  Basic system commands

RM=rm -f
CP=cp
MV=mv
SED=sed
CHMOD=chmod

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

#  The installation

lancelot:
	echo "FORTRAN=" $(FORTRAN)
	echo " FFLAGS=" $(FFLAGS)
	cd ./src ; $(MAKE) install \
                        PRECISION=$(PRECISION) SIZE=$(SIZE) \
                        FORTRAN="$(FORTRAN)" FFLAGS="$(FFLAGS)" \
                        RM="$(RM)" CP="$(CP)" MV="$(MV)" TMP="$(TMP)"
	$(SED) "s/^set FORTRAN=.*/set FORTRAN='$(FORTRAN)'/g" \
           ./bin/lan_default > $(TMP)/temporary_file
	$(SED) "s/^set FFLAGS=.*/set FFLAGS='$(FFLAGS)'/g" \
           $(TMP)/temporary_file > ./bin/lan

	$(SED) "s/^FORTRAN=.*/FORTRAN=$(FORTRAN)/g" \
           ./bin/compil_default > $(TMP)/temporary_file
	$(SED) "s/^FFLAGS=.*/FFLAGS=$(FFLAGS)/g" \
           $(TMP)/temporary_file > ./bin/compil
	$(CHMOD) a+x ./bin/lan
	$(RM) $(TMP)/temporary_file

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
	@ echo " -   SED=sed command (/bin/sed)"
	@ echo " -   CHMOD=chmod command (/bin/chmod)"
	@ echo " -   TMP=temporary workspace directory, not ./ (/tmp)"

#  Clean up afterwards

clean:
	cd ./src ; $(MAKE) clean PRECISION=$(PRECISION)
	cd ./bin ; $(RM) lan compil
