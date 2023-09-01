# 
# makefile for GMCALC v1.5.3          Aug 2, 2022
# http://people.physics.carleton.ca/~logan/gmcalc/
#

# put the path to LoopTools here (tested with LoopTools 2.15 & 2.16)
#LT = $(HOME)/Documents/work/looptools/LoopTools-2.15/x86_64-Darwin
LT = $(HOME)/Documents/work/looptools/LoopTools-2.16/arm64-Darwin
LIB = lib   # use this if libooptools.a is in $LT/lib
# LIB = lib64 # use this if libooptools.a is in $LT/lib64

# put the path to HiggsSignals 2 and HiggsBounds 5 here
HS2 = $(HOME)/Documents/work/HiggsSignals-2.2.1beta
# HB5 = $(HOME)/Documents/work/HiggsBounds-5.3.0beta
HB5 = $(HOME)/Documents/work/HiggsBounds-5.3.0beta

FF = gfortran
# Note that for our interface to HiggBounds/HiggsSignals to work properly,
# gfortran version 5 or later must be used.

gmpoint: gmpoint.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f src/gmconstr.f
	$(FF) -I$(LT)/include/ \
	gmpoint.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f src/gmconstr.f \
	-L$(LT)/$(LIB) -looptools \
	-o gmpoint.x

gmpoint-nolt: gmpoint.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops-dummy.F  \
src/gamvvof_final.f src/vegas.f src/gmconstr.f
	$(FF) gmpoint.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f src/gmconstr.f \
	-o gmpoint-nolt.x

gmscan: gmscan.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f src/gmconstr.f
	$(FF) -I$(LT)/include/ \
	gmscan.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f src/gmconstr.f \
	-L$(LT)/$(LIB) -looptools \
	-o gmscan.x

gmscan-nolt: gmscan.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops-dummy.F \
src/gamvvof_final.f src/vegas.f src/gmconstr.f
	$(FF) gmscan.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f src/gmconstr.f \
	-o gmscan-nolt.x

gmmg5: gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
src/gmspectrum.f src/gmqcd.f src/gmdecays.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f src/lininterp.f
	$(FF) -I$(LT)/include/ \
	gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
	src/gmspectrum.f src/gmqcd.f src/gmdecays.f \
	src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f src/lininterp.f \
	-L$(LT)/$(LIB) -looptools \
	-o gmmg5.x

gmmg5-nolt: gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
src/gmspectrum.f src/gmqcd.f src/gmdecays.f src/heteroloops-dummy.F \
src/gamvvof_final.f src/vegas.f 
	$(FF) gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
	src/gmspectrum.f src/gmqcd.f src/gmdecays.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f \
	-o gmmg5-nolt.x

gmhb5: gmhb5.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f src/gmconstr.f src/gmhbhs.f
	$(FF) -I$(LT)/include/ -I$(HB5) -I$(HS2) \
	gmhb5.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f src/gmconstr.f src/gmhbhs.f \
	-L$(LT)/$(LIB) -L$(HB5) -L$(HS2) -looptools -lHB -lHS \
	-o gmhb5.x

gmhb5-nolt: gmhb5.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops-dummy.F \
src/gamvvof_final.f src/vegas.f src/gmconstr.f src/gmhbhs.f
	$(FF) -I$(HB5) -I$(HS2) \
	gmhb5.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f src/gmconstr.f src/gmhbhs.f \
	 -L$(HB5) -L$(HS2) -lHB -lHS \
	-o gmhb5-nolt.x
