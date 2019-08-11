#
#  Produces object files for the LANCELOT optimizers
#
#  version for Unix systems using gfortran

#  Philippe Toint, November 25th, 1991, for CGT Productions.
#  This version August 8th, 2019

#  Basic system commands

RM=rm -f
CP=cp

#  Set directory for temporary object files
#  (do NOT use TMP = . which would cause file overwriting)

TMP=/tmp

#  Default precision

PRECISION=double

# Set directories for single and double precision object files

OD=../../objects/$(PRECISION)

#  Default sizing program

SIZE=large

#  Compiler flags

FORTRAN=gfortran
FFLAGS= -O3 -ffixed-form

#  Filtering rule

.SUFFIXES: .o .f
.f.o:
	../$(PRECISION) < $< | ../$(SIZE) > $(TMP)/$<
	$(FORTRAN) -c $(FFLAGS) -o $(OD)/$*.o $(TMP)/$<
	$(RM) $(TMP)/$<

#  List of object files to create

OB = asmbl.o   auglg.o   dumbr.o  cauch.o   cg.o     drche.o  lance.o  \
     drchg.o   elgrd.o   frntl.o  gtpgr.o   local.o  hsprd.o  inxac.o  \
     sort.o    linpac.o  assl.o   misc.o    mdchl.o  precn.o  others.o \
     runlan.o  sbmin.o   datal.o  duma27.o  duma31.o dumc29.o speci.o \
     scaln.o   initw.o   bndsl.o

#  The compilations

lancelot: lancelot_$(PRECISION)

lancelot_single:
	$(MAKE) $(OB) PRECISION=single SIZE=$(SIZE) TMP=$(TMP) \
                  OD=../../objects/single

lancelot_double:
	$(MAKE) $(OB) PRECISION=double SIZE=$(SIZE) TMP=$(TMP) \
                  OD=../../objects/double

# HELP!

help:
	@ echo "Possible choices are:"
	@ echo " - help: this information"
	@ echo " - lancelot_single: the single precision optimizer"
	@ echo " - lancelot_double: the double precision optimizer"
	@ echo " -  options (default)"
	@ echo " -   FORTRAN=your fortran 77 compiler (gfortran)"
	@ echo " -   FFLAGS=your fortran compiler flags (-O3 -ffixed-form)"
	@ echo " -   SIZE=small medium large huge (large)"
	@ echo " -   PRECISION=single double (double)"
	@ echo " -   RM=rm command with force (/bin/rm -f)"
	@ echo " -   CP=cp command (/bin/cp)"
	@ echo " -   TMP=temporary workspace directory, not ./ (/tmp)"

#  Clean up afterwards

clean:
	$(RM) ../../objects/double/*
