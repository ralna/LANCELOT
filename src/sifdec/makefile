#
#  Compiling the SIF file decoder.
#
#  version for Unix systems using gfortran

#  Ph. Toint, November 21rst, 1991, for CGT Productions.
#  This version August 8th, 2019


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
FFLAGSLD=

#  Filtering rule

.SUFFIXES: .o .f
.f.o:
	../$(PRECISION) < $< | ../$(SIZE) > $(TMP)/$<
	$(FORTRAN) -c $(FFLAGS) -o $(TMP)/$*.o $(TMP)/$<
	$(RM) $(TMP)/$<

#  list of object files created

OBJ=gps.o   decode.o  rest.o  inlanc.o  utils.o  makefn.o  makegr.o  printp.o \
    local.o runsd.o

#  The compilations

sifdec: sifdec_$(PRECISION)

sifdec_single:
	$(MAKE) $(OBJ) PRECISION=single SIZE=$(SIZE) TMP=$(TMP)
	cd $(TMP); $(FORTRAN) $(FFLAGSLD) -o sifdec_single $(OBJ); $(RM) $(OBJ)
	$(MV) $(TMP)/sifdec_single ../../bin/sifdec_single

sifdec_double:
	$(MAKE) $(OBJ) PRECISION=double SIZE=$(SIZE) TMP=$(TMP)
	cd $(TMP); $(FORTRAN) $(FFLAGSLD) -o sifdec_double $(OBJ); $(RM) $(OBJ)
	$(MV) $(TMP)/sifdec_double ../../bin/sifdec_double

#  HELP!

help:
	@ echo "Possible choices are:"
	@ echo " - help: this information"
	@ echo " - sifdec_single: the single precision problem decoder"
	@ echo " - sifdec_double: the double precision problem decoder"
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
	$(RM) ../../bin/sifdec_$(PRECISION)
