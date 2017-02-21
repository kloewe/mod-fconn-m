#-----------------------------------------------------------------------------
# File    : makefile
# Contents: build MEX-files
# Author  : Kristian Loewe
#-----------------------------------------------------------------------------
.SUFFIXES:                                   # remove built-in rules
MAKEFLAGS   += -r

SHELL       = /bin/bash
DOTDIR      = ../../dot
CORRDIR     = ../../corr
CPUINFODIR  = ../../cpuinfo
FCONNDIR    = ../../fconn
STATSDIR    = ../../stats

CC          = gcc
CFLAGS      = -std=c99 \
              -Wall -Wextra -Wno-unused-parameter -Wconversion -Wshadow \
              -pedantic -fPIC
COPTIMFLAGS = -O2

MATLABROOT  = $(dir $(realpath $(shell which matlab)))
MEXCC       = $(realpath $(MATLABROOT))/mex -largeArrayDims
MEXLD       = $(realpath $(MATLABROOT))/mex -largeArrayDims $(MEXOPTSLD)

SUF         = mexa64
_OBJDIR     = obj/$(shell uname -m)/matlab
OBJDIR      = ../$(_OBJDIR)
_DUMMY      := $(shell mkdir -p $(OBJDIR); mkdir -p $(OBJDIR))
OUTDIR      = ../private
_DUMMY      := $(shell mkdir -p $(OUTDIR); mkdir -p $(OUTDIR))

#-----------------------------------------------------------------------------
# Link MEX-files
#-----------------------------------------------------------------------------
all: nodedeg_flt nodedeg_dbl edgestats_flt edgestats_dbl

nodedeg_flt:                   $(OUTDIR)/mxNodeDegFlt.$(SUF)
$(OUTDIR)/mxNodeDegFlt.$(SUF): $(FCONNDIR)/$(_OBJDIR)/nodedeg_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegFlt.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/nodedeg_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegFlt.o -output $@

edgestats_flt:                 $(OUTDIR)/mxFcmFlt.$(SUF)
$(OUTDIR)/mxFcmFlt.$(SUF):     $(FCONNDIR)/$(_OBJDIR)/edgestats_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_flt.o \
                               $(DOTDIR)/$(_OBJDIR)/dot.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_naive.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_sse2.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avx.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avxfma.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmFlt.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/edgestats_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_flt.o \
                               $(DOTDIR)/$(_OBJDIR)/dot.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_naive.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_sse2.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avx.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avxfma.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmFlt.o -output $@

nodedeg_dbl:                   $(OUTDIR)/mxNodeDegDbl.$(SUF)
$(OUTDIR)/mxNodeDegDbl.$(SUF): $(FCONNDIR)/$(_OBJDIR)/nodedeg_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegDbl.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/nodedeg_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegDbl.o -output $@

edgestats_dbl:                 $(OUTDIR)/mxFcmDbl.$(SUF)
$(OUTDIR)/mxFcmDbl.$(SUF):     $(FCONNDIR)/$(_OBJDIR)/edgestats_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_dbl.o \
                               $(DOTDIR)/$(_OBJDIR)/dot.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_naive.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_sse2.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avx.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avxfma.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmDbl.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/edgestats_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_dbl.o \
                               $(DOTDIR)/$(_OBJDIR)/dot.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_naive.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_sse2.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avx.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_avxfma.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmDbl.o -output $@

#-----------------------------------------------------------------------------
# fconn-m objects
#-----------------------------------------------------------------------------
mxFcmFlt.o:                $(OBJDIR)/mxFcmFlt.o
$(OBJDIR)/mxFcmFlt.o:      $(FCONNDIR)/src/edgestats.h \
                           $(FCONNDIR)/src/fcmat.h $(FCONNDIR)/src/matrix.h \
                           $(STATSDIR)/src/stats.h
$(OBJDIR)/mxFcmFlt.o:      mxFcm.c makefile
	$(realpath $(MATLABROOT))/mex -largeArrayDims \
    CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' \
    -DNDEBUG -DREAL=float \
    -I$(FCONNDIR)/src -I$(STATSDIR)/src \
    -c mxFcm.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxFcm.o $(OBJDIR)/mxFcmFlt.o

mxNodeDegFlt.o:            $(OBJDIR)/mxNodeDegFlt.o
$(OBJDIR)/mxNodeDegFlt.o:  $(FCONNDIR)/src/nodedeg.h $(FCONNDIR)/src/fcmat.h
$(OBJDIR)/mxNodeDegFlt.o:  mxNodeDeg.c makefile
	$(realpath $(MATLABROOT))/mex -largeArrayDims \
    CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' \
    -DNDEBUG -DREAL=float \
    -I$(FCONNDIR)/src -I$(STATSDIR)/src \
    -c mxNodeDeg.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxNodeDeg.o $(OBJDIR)/mxNodeDegFlt.o

mxFcmDbl.o:                $(OBJDIR)/mxFcmDbl.o
$(OBJDIR)/mxFcmDbl.o:      $(FCONNDIR)/src/fcmat.h $(FCONNDIR)/src/matrix.h \
                           $(FCONNDIR)/src/edgestats.h $(STATSDIR)/src/stats.h
$(OBJDIR)/mxFcmDbl.o:      mxFcm.c makefile
	$(realpath $(MATLABROOT))/mex -largeArrayDims \
    CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' \
    -DNDEBUG -DREAL=double \
    -I$(FCONNDIR)/src -I$(STATSDIR)/src \
    -c mxFcm.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxFcm.o $(OBJDIR)/mxFcmDbl.o

mxNodeDegDbl.o:            $(OBJDIR)/mxNodeDegDbl.o
$(OBJDIR)/mxNodeDegDbl.o:  $(FCONNDIR)/src/nodedeg.h $(FCONNDIR)/src/fcmat.h
$(OBJDIR)/mxNodeDegDbl.o:  mxNodeDeg.c makefile
	$(realpath $(MATLABROOT))/mex -largeArrayDims \
    CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' \
    -DNDEBUG -DREAL=double \
    -I$(FCONNDIR)/src -I$(STATSDIR)/src \
    -c mxNodeDeg.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxNodeDeg.o $(OBJDIR)/mxNodeDegDbl.o

#-----------------------------------------------------------------------------
# dot objects
#-----------------------------------------------------------------------------
$(DOTDIR)/$(_OBJDIR)/dot_naive.o:
	cd $(DOTDIR)/src; $(MAKE) -f makefile-mex dot_naive.o

$(DOTDIR)/$(_OBJDIR)/dot_sse2.o:
	cd $(DOTDIR)/src; $(MAKE) -f makefile-mex dot_sse2.o

$(DOTDIR)/$(_OBJDIR)/dot_avx.o:
	cd $(DOTDIR)/src; $(MAKE) -f makefile-mex dot_avx.o

$(DOTDIR)/$(_OBJDIR)/dot_avxfma.o:
	cd $(DOTDIR)/src; $(MAKE) -f makefile-mex dot_avxfma.o

$(DOTDIR)/$(_OBJDIR)/dot.o:
	cd $(DOTDIR)/src; $(MAKE) -f makefile-mex dot.o

#-----------------------------------------------------------------------------
# stats objects
#-----------------------------------------------------------------------------
$(STATSDIR)/$(_OBJDIR)/stats_flt.o:
	cd $(STATSDIR)/src; $(MAKE) -f makefile-mex stats_flt.o

$(STATSDIR)/$(_OBJDIR)/stats_dbl.o:
	cd $(STATSDIR)/src; $(MAKE) -f makefile-mex stats_dbl.o

#-----------------------------------------------------------------------------
# cpuinfo objects
#-----------------------------------------------------------------------------
$(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o:
	cd $(CPUINFODIR)/src; $(MAKE) -f makefile-mex cpuinfo.o

#-----------------------------------------------------------------------------
# corr objects
#-----------------------------------------------------------------------------
$(CORRDIR)/$(_OBJDIR)/pcc_flt.o:
	cd $(CORRDIR)/src; $(MAKE) -f makefile-mex pcc_flt.o

$(CORRDIR)/$(_OBJDIR)/tetracc_flt.o:
	cd $(CORRDIR)/src; $(MAKE) -f makefile-mex tetracc_flt.o

$(CORRDIR)/$(_OBJDIR)/binarize_flt.o:
	cd $(CORRDIR)/src; $(MAKE) -f makefile-mex binarize_flt.o

$(CORRDIR)/$(_OBJDIR)/pcc_dbl.o:
	cd $(CORRDIR)/src; $(MAKE) -f makefile-mex pcc_dbl.o

$(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o:
	cd $(CORRDIR)/src; $(MAKE) -f makefile-mex tetracc_dbl.o

$(CORRDIR)/$(_OBJDIR)/binarize_dbl.o:
	cd $(CORRDIR)/src; $(MAKE) -f makefile-mex binarize_dbl.o

#-----------------------------------------------------------------------------
# fconn objects
#-----------------------------------------------------------------------------
$(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex fcmat_flt.o

$(FCONNDIR)/$(_OBJDIR)/matrix_flt.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex matrix_flt.o

$(FCONNDIR)/$(_OBJDIR)/edgestats_flt.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex edgestats_flt.o

$(FCONNDIR)/$(_OBJDIR)/nodedeg_flt.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex nodedeg_flt.o

$(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex fcmat_dbl.o

$(FCONNDIR)/$(_OBJDIR)/matrix_dbl.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex matrix_dbl.o

$(FCONNDIR)/$(_OBJDIR)/edgestats_dbl.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex edgestats_dbl.o

$(FCONNDIR)/$(_OBJDIR)/nodedeg_dbl.o:
	cd $(FCONNDIR)/src; $(MAKE) -f makefile-mex nodedeg_dbl.o
