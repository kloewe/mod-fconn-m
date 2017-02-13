#-----------------------------------------------------------------------------
# File    : makefile
# Contents: build MEX-files
#           NB: tested with Matlab R2011b
#
#           make DEBUG=1
#           make nodedeg
#           make nodedeg_flt
#           make nodedeg_dbl
#
# Author  : Kristian Loewe
#-----------------------------------------------------------------------------
.SUFFIXES:                                   # remove built-in rules
MAKEFLAGS    += -r

SHELL       = /bin/bash
CORRDIR     = ../../corr/src
CPUINFODIR  = ../../cpuinfo/src
FCONNDIR    = ../../fconn/src
STATSDIR    = ../../stats/src

CC          = gcc
CFLAGS      = -std=c99 -march=native \
              -Wall -Wextra -Wno-unused-parameter -Wconversion -Wshadow \
              -pedantic -fPIC
DEFS        = -DPAIRSPLIT -DRECTGRID -DSAFETHREAD

MATLABROOT  = $(dir $(realpath $(shell which matlab)))
MEXCC       = $(realpath $(MATLABROOT))/mex -largeArrayDims
MEXLD       = $(realpath $(MATLABROOT))/mex -largeArrayDims $(MEXOPTSLD)

DEBUG       ?= 0
ifeq ($(DEBUG), 1)
  MEXCC       += -v
  CFLAGS      += -g
  COPTIMFLAGS =  -O2
else
  MEXCC       += -O
  COPTIMFLAGS =  -O3 -funroll-loops
  DEFS        += -DNDEBUG
endif

MEXCC       += $(MEXOPTSCC)
MEXOPTSCC   := CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' $(DEFS)
MEXOPTSLD   = LDOPTIMFLAGS='$(COPTIMFLAGS)'
CFLAGS      += $(COPTIMFLAGS) $(DEFS)

INCS        = -I$(CORRDIR) -I$(CPUINFODIR) -I$(FCONNDIR) -I$(STATSDIR)
SUF         = mexa64
OUTDIR      = ../private

LIBS        = -lm -lpthread

PRGS        = nodedeg edgestats

#-----------------------------------------------------------------------------
# Build Programs
#-----------------------------------------------------------------------------
all: $(PRGS)

$(PRGS): | $(OUTDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# nodedeg
nodedeg_flt = $(OUTDIR)/mxNodeDegFlt.$(SUF)
nodedeg_dbl = $(OUTDIR)/mxNodeDegDbl.$(SUF)
nodedeg: $(nodedeg_flt) $(nodedeg_dbl)
	

$(nodedeg_flt): DREAL = -DREAL=float
$(nodedeg_dbl): DREAL = -DREAL=double

$(nodedeg_flt) $(nodedeg_dbl): \
  CFILES = $(FCONNDIR)/nodedeg.c $(FCONNDIR)/fcmat.c $(CORRDIR)/pcc.c \
    $(CORRDIR)/tetracc.c $(CORRDIR)/binarize.c $(CPUINFODIR)/cpuinfo.c
  PREREQ = $(FCONNDIR)/nodedeg.* $(FCONNDIR)/fcmat.* $(CORRDIR)/pcc.* \
    $(CORRDIR)/tetracc.* $(CORRDIR)/binarize.* $(CPUINFODIR)/cpuinfo.*

$(nodedeg_flt) $(nodedeg_dbl): $(PREREQ)
$(nodedeg_flt) $(nodedeg_dbl): mxNodeDeg.c makefile
	$(MEXCC) $(DREAL) $(INCS) $(LIBS) $< $(CFILES) -outdir $(OUTDIR) -output $@

# mean, std, var, cmp, corr
edgestats_flt = $(OUTDIR)/mxFcmFlt.$(SUF)
edgestats_dbl = $(OUTDIR)/mxFcmDbl.$(SUF)
edgestats: $(edgestats_flt) $(edgestats_dbl)
	

$(edgestats_flt): DREAL = -DREAL=float
$(edgestats_dbl): DREAL = -DREAL=double

$(edgestats_flt) $(edgestats_dbl): \
  CFILES = $(FCONNDIR)/edgestats.c $(FCONNDIR)/fcmat.c $(FCONNDIR)/matrix.c \
    $(CORRDIR)/pcc.c $(CORRDIR)/tetracc.c $(CORRDIR)/binarize.c \
    $(CPUINFODIR)/cpuinfo.c $(STATSDIR)/stats.c
  PREREQ = $(FCONNDIR)/edgestats.* $(FCONNDIR)/fcmat.*  $(FCONNDIR)/matrix.* \
    $(CORRDIR)/pcc.* $(CORRDIR)/tetracc.* $(CORRDIR)/binarize.* \
    $(CPUINFODIR)/cpuinfo.* $(STATSDIR)/stats.*

$(edgestats_flt) $(edgestats_dbl): $(PREREQ)
$(edgestats_flt) $(edgestats_dbl): mxFcm.c makefile
	$(MEXCC) $(DREAL) $(INCS) $(LIBS) $< $(CFILES) -outdir $(OUTDIR) -output $@
