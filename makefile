# 
# makefile for GMCALC v1.3.0          Aug 30, 2017
# http://people.physics.carleton.ca/~logan/gmcalc/
#

# LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools-2.14/x86_64-Darwin
LT = /home/ycwu/Workings/SUPPORTS/LoopTools/2.15

SRCDIR := src
FF = gfortran
FFLAG = -I$(LT)/include
FLIBS = -L$(LT)/lib -looptools

SRC = $(wildcard $(SRCDIR)/*.f)
SRCLT = $(SRC) $(SRCDIR/heteroloops.F)
SRCNOLT = $(SRC) $(SRCDIR/heteroloops-dummy.F)
OBJ = $(patsubst $(SRCDIR)/%.f, $(SRCDIR)/%.o, $(SRC))
OBJLT = $(OBJ) $(SRCDIR)/heteroloops.o
OBJNOLT = $(OBJ) $(SRCDIR)/heteroloops-dummy.o

all: gmpoint.x gmscan.x gmmg5.x

all-nolt: gmpoint-nolt.x gmscan-nolt.x gmmg5-nolt.x

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

%.x: %.f $(OBJLT)
	$(FF) $(FFLAG) -o $@ $< $(OBJLT) $(FLIBS)

$(SRCDIR)/%.o: $(SRCDIR)/%.f
	$(FF) $(FFLAG) -c $< -o $@

$(SRCDIR)/heteroloops.o: $(SRCDIR)/heteroloops.F
	$(FF) $(FFLAG) -c $< -o $@

$(SRCDIR)/heteroloops-dummy.o: $(SRCDIR)/heteroloops-dummy.F
	$(FF) $(FFLAG) -c $< -o $@

clean:
	rm -f *.x
	rm -f $(SRCDIR)/*.o
