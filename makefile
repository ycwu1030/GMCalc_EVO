# 
# makefile for GMCALC v1.4.1         Sept 30, 2018
# http://people.physics.carleton.ca/~logan/gmcalc/
#

# tested with LoopTools 2.15
LT = $(HOME)/Documents/work/looptools/LoopTools/x86_64-Darwin

FF = gfortran

gmpoint: gmpoint.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f
	$(FF) -I$(LT)/include/ \
	gmpoint.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f \
	-L$(LT)/lib -looptools \
	-o gmpoint.x

gmpoint-nolt: gmpoint.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops-dummy.F  \
src/gamvvof_final.f src/vegas.f
	$(FF) gmpoint.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f \
	-o gmpoint-nolt.x

gmscan: gmscan.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f
	$(FF) -I$(LT)/include/ \
	gmscan.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f \
	-L$(LT)/lib -looptools \
	-o gmscan.x

gmscan-nolt: gmscan.f src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f src/heteroloops-dummy.F \
src/gamvvof_final.f src/vegas.f
	$(FF) gmscan.f src/gmdecays.f src/gmprint.f \
	src/gmindir.f src/gminit.f src/gmread.f src/gmspectrum.f \
	src/lininterp.f src/gmutils.f src/gmqcd.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f \
	-o gmscan-nolt.x

gmmg5: gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
src/gmspectrum.f src/gmqcd.f src/gmdecays.f src/heteroloops.F \
src/gamvvof_final.f src/vegas.f
	$(FF) -I$(LT)/include/ \
	gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
	src/gmspectrum.f src/gmqcd.f src/gmdecays.f \
	src/heteroloops.F \
	src/gamvvof_final.f src/vegas.f \
	-L$(LT)/lib -looptools \
	-o gmmg5.x

gmmg5-nolt: gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
src/gmspectrum.f src/gmqcd.f src/gmdecays.f src/heteroloops-dummy.F \
src/gamvvof_final.f src/vegas.f
	$(FF) gmmg5.f src/gmprint.f src/gminit.f src/gmread.f \
	src/gmspectrum.f src/gmqcd.f src/gmdecays.f \
	src/heteroloops-dummy.F \
	src/gamvvof_final.f src/vegas.f \
	-o gmmg5-nolt.x


