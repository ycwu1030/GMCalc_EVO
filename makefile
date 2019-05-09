#
# makefile for GMCALC v1.5.3          Aug 2, 2022
# http://people.physics.carleton.ca/~logan/gmcalc/
#

# put the path to LoopTools here (tested with LoopTools 2.15 & 2.16)
LTLIB = $(HOME)/Library/Mathematica/Applications/LoopTools/arm64-Darwin/lib

# put the path to HiggsSignals 2 and HiggsBounds 5 here
HS2 = $(HOME)/Documents/work/HiggsSignals-2.2.1beta
# HB5 = $(HOME)/Documents/work/HiggsBounds-5.3.0beta
HB5 = $(HOME)/Documents/work/HiggsBounds-5.3.0beta

SRCDIR := src
FF = gfortran
FFLAG = -I$(LTLIB)/../include
FLIBS = -L$(LTLIB) -looptools
FFLAGHBHS = -I$(HB5) -I$(HS2)
FLIBSHBHS = -L$(HS2) -L$(HB5) -lHB -lHS

SRC = src/gmdecays.f src/gmprint.f src/gmindir.f \
src/gminit.f src/gmread.f src/gmspectrum.f src/lininterp.f \
src/gmutils.f src/gmqcd.f \
src/gamvvof_final.f src/vegas.f src/gmconstr.f
SRCHBHS = $(SRCDIR)/gmhbhs.f
SRCLT = $(SRC) $(SRCDIR/heteroloops.F)
SRCNOLT = $(SRC) $(SRCDIR/heteroloops-dummy.F)
OBJ = $(patsubst $(SRCDIR)/%.f, $(SRCDIR)/%.o, $(SRC))
OBJHBHS = $(SRCDIR)/gmhbhs.o
OBJLT = $(OBJ) $(SRCDIR)/heteroloops.o
OBJNOLT = $(OBJ) $(SRCDIR)/heteroloops-dummy.o

all: gmpoint.x gmscan.x gmmg5.x gmhb5.x

all-nolt: gmpoint-nolt.x gmscan-nolt.x gmmg5-nolt.x gmhb5-nolt.x

gmpoint.x: gmpoint.f $(OBJLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJLT) $(FLIBS)

gmpoint-nolt.x: gmpoint.f $(OBJNOLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJNOLT) $(FLIBS)

gmscan.x: gmscan.f $(OBJLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJLT) $(FLIBS)

gmscan-nolt.x: gmscan.f $(OBJNOLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJNOLT) $(FLIBS)

gmmg5.x: gmmg5.f $(OBJLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJLT) $(FLIBS)

gmmg5-nolt.x: gmmg5.f $(OBJNOLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJNOLT) $(FLIBS)

gmhb5.x: gmhb5.f $(OBJLT) $(OBJHBHS)
	$(FF) $(FFLAG) $(FFLAGHBHS) -o $@ $< $(OBJLT) $(OBJHBHS) $(FLIBS)

gmhb5-nolt.x: gmhb5.f $(OBJNOLT) $(OBJHBHS)
	$(FF) $(FFLAG) $(FFLAGHBHS) -o $@ $< $(OBJNOLT) $(OBJHBHS) $(FLIBS)

%.x: %.f $(OBJLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJLT) $(FLIBS)

$(SRCDIR)/%.o: $(SRCDIR)/%.f
	$(FF) $(FFLAG) -c $< -o $@

$(SRCDIR)/gmhbhs.o: $(SRCDIR)/gmhbhs.f
	$(FF) $(FFLAGHBHS) -c $< -o $@

$(SRCDIR)/heteroloops.o: $(SRCDIR)/heteroloops.F
	$(FF) $(FFLAG) -c $< -o $@

$(SRCDIR)/heteroloops-dummy.o: $(SRCDIR)/heteroloops-dummy.F
	$(FF) $(FFLAG) -c $< -o $@

clean:
	rm -f *.x
	rm -f $(SRCDIR)/*.o
