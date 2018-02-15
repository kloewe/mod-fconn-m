#-----------------------------------------------------------------------------
# File    : makefile
# Contents: build MEX-files
# Author  : Kristian Loewe
#
# Usage   : make
#           MEX_FLAGS='-v' make
#           make | grep -v 'Warning.*gcc version'
#           make -B | grep -v 'Warning.*gcc version'
#           DEBUG=1 make -B | grep -v 'Warning.*gcc version'
#-----------------------------------------------------------------------------
.SUFFIXES:
MAKEFLAGS   += -r

CC          ?= gcc
CFBASE       = -std=c99 -Wall -Wextra -Wno-unused-parameter -Wconversion \
               -Wshadow -pedantic
DEFS        ?=
MEX_FLAGS   ?=

DEBUG       ?= 0
ifeq ($(DEBUG), 1)
  CFBASE    += -g
  CFOPT     ?= -O0
else
  CFOPT     ?= -O2
  DEFS      += -DNDEBUG
endif
CFLAGS       = $(CFBASE) -fPIC $(DEFS)

MATLABROOT   = $(dir $(realpath $(shell which matlab)))
MEXCC        = $(realpath $(MATLABROOT))/mex -largeArrayDims $(MEX_FLAGS) \
               CFLAGS='$(CFLAGS)'
MEXLD        = $(realpath $(MATLABROOT))/mex -largeArrayDims $(MEX_FLAGS)
SUF          = mexa64

_OBJDIR      = obj/$(shell uname -m)/matlab
OBJDIR       = ../$(_OBJDIR)
OUTDIR       = ../private
_DUMMY      := $(shell mkdir -p $(OBJDIR); mkdir -p $(OUTDIR))

#-----------------------------------------------------------------------------

eq = $(if $(or $(1),$(2)),$(and $(findstring $(1),$(2)),\
                                $(findstring $(2),$(1))),1)

#-----------------------------------------------------------------------------

CORRDIR      = ../../corr
CPUINFODIR   = ../../cpuinfo
DOTDIR       = ../../dot
FCONNDIR     = ../../fconn
STATSDIR     = ../../stats

INCS         = -I$(FCONNDIR)/src -I$(STATSDIR)/src

MEX_FILES    = mxNodeDegFlt mxFcmFlt \
               mxNodeDegDbl mxFcmDbl
#-----------------------------------------------------------------------------
# Link MEX files
#-----------------------------------------------------------------------------
all: $(MEX_FILES)

mxNodeDegFlt:                  $(OUTDIR)/mxNodeDegFlt.$(SUF)
$(OUTDIR)/mxNodeDegFlt.$(SUF): $(FCONNDIR)/$(_OBJDIR)/nodedeg_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegFlt.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/nodedeg_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegFlt.o -output $@

mxFcmFlt:                      $(OUTDIR)/mxFcmFlt.$(SUF)
$(OUTDIR)/mxFcmFlt.$(SUF):     $(FCONNDIR)/$(_OBJDIR)/edgestats_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmFlt.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/edgestats_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmFlt.o -output $@

mxNodeDegDbl:                  $(OUTDIR)/mxNodeDegDbl.$(SUF)
$(OUTDIR)/mxNodeDegDbl.$(SUF): $(FCONNDIR)/$(_OBJDIR)/nodedeg_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegDbl.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/nodedeg_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxNodeDegDbl.o -output $@

mxFcmDbl:                      $(OUTDIR)/mxFcmDbl.$(SUF)
$(OUTDIR)/mxFcmDbl.$(SUF):     $(FCONNDIR)/$(_OBJDIR)/edgestats_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
                               $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                               $(OBJDIR)/mxFcmDbl.o makefile
	$(MEXLD)                     $(FCONNDIR)/$(_OBJDIR)/edgestats_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o \
                               $(FCONNDIR)/$(_OBJDIR)/matrix_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                               $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                               $(STATSDIR)/$(_OBJDIR)/stats_all.o \
                               $(DOTDIR)/$(_OBJDIR)/dot_all.o \
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
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=float $(INCS) -I$(DOTDIR)/src \
    -c mxFcm.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxFcm.o $(OBJDIR)/mxFcmFlt.o

mxNodeDegFlt.o:            $(OBJDIR)/mxNodeDegFlt.o
$(OBJDIR)/mxNodeDegFlt.o:  $(FCONNDIR)/src/nodedeg.h $(FCONNDIR)/src/fcmat.h
$(OBJDIR)/mxNodeDegFlt.o:  mxNodeDeg.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=float $(INCS) \
    -c mxNodeDeg.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxNodeDeg.o $(OBJDIR)/mxNodeDegFlt.o

mxFcmDbl.o:                $(OBJDIR)/mxFcmDbl.o
$(OBJDIR)/mxFcmDbl.o:      $(FCONNDIR)/src/fcmat.h $(FCONNDIR)/src/matrix.h \
                           $(FCONNDIR)/src/edgestats.h $(STATSDIR)/src/stats.h
$(OBJDIR)/mxFcmDbl.o:      mxFcm.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=double $(INCS) -I$(DOTDIR)/src \
    -c mxFcm.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxFcm.o $(OBJDIR)/mxFcmDbl.o

mxNodeDegDbl.o:            $(OBJDIR)/mxNodeDegDbl.o
$(OBJDIR)/mxNodeDegDbl.o:  $(FCONNDIR)/src/nodedeg.h $(FCONNDIR)/src/fcmat.h
$(OBJDIR)/mxNodeDegDbl.o:  mxNodeDeg.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=double $(INCS) \
    -c mxNodeDeg.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxNodeDeg.o $(OBJDIR)/mxNodeDegDbl.o

#-----------------------------------------------------------------------------
# cpuinfo objects
#-----------------------------------------------------------------------------
MKARGS_CPUINFO = -C $(CPUINFODIR)/src -f makefile-mex cpuinfo.o
REMAKE_CPUINFO = $(shell $(MAKE) -q $(MKARGS_CPUINFO) > /dev/null 2>&1; \
                 echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CPUINFO),1),$(findstring B,$(MAKEFLAGS))),)
$(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o: FORCE
	$(MAKE) $(MKARGS_CPUINFO)
FORCE:
endif

#-----------------------------------------------------------------------------
# dot objects
#-----------------------------------------------------------------------------
$(DOTDIR)/$(_OBJDIR)/dot_all.o:
	$(MAKE) -C $(DOTDIR)/src -f makefile-mex dot_all.o

#-----------------------------------------------------------------------------
# stats objects
#-----------------------------------------------------------------------------
$(STATSDIR)/$(_OBJDIR)/stats_all.o:
	$(MAKE) -C $(STATSDIR)/src -f makefile-mex stats_all.o

#-----------------------------------------------------------------------------
# corr objects
#-----------------------------------------------------------------------------
MKARGS_CORR1F = -C $(CORRDIR)/src -f makefile-mex pcc_flt.o
REMAKE_CORR1F = $(shell $(MAKE) -q $(MKARGS_CORR1F) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR1F),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/pcc_flt.o: FORCE
	$(MAKE) $(MKARGS_CORR1F)
FORCE:
endif

MKARGS_CORR2F = -C $(CORRDIR)/src -f makefile-mex tetracc_flt.o
REMAKE_CORR2F = $(shell $(MAKE) -q $(MKARGS_CORR2F) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR2F),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/tetracc_flt.o: FORCE
	$(MAKE) $(MKARGS_CORR2F)
FORCE:
endif

MKARGS_CORR3F = -C $(CORRDIR)/src -f makefile-mex binarize_flt.o
REMAKE_CORR3F = $(shell $(MAKE) -q $(MKARGS_CORR3F) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR3F),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/binarize_flt.o: FORCE
	$(MAKE) $(MKARGS_CORR3F)
FORCE:
endif
#-----------------------------------------------------------------------------
MKARGS_CORR1D = -C $(CORRDIR)/src -f makefile-mex pcc_dbl.o
REMAKE_CORR1D = $(shell $(MAKE) -q $(MKARGS_CORR1D) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR1D),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/pcc_dbl.o: FORCE
	$(MAKE) $(MKARGS_CORR1D)
FORCE:
endif

MKARGS_CORR2D = -C $(CORRDIR)/src -f makefile-mex tetracc_dbl.o
REMAKE_CORR2D = $(shell $(MAKE) -q $(MKARGS_CORR2D) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR2D),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o: FORCE
	$(MAKE) $(MKARGS_CORR2D)
FORCE:
endif

MKARGS_CORR3D = -C $(CORRDIR)/src -f makefile-mex binarize_dbl.o
REMAKE_CORR3D = $(shell $(MAKE) -q $(MKARGS_CORR3D) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR3D),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/binarize_dbl.o: FORCE
	$(MAKE) $(MKARGS_CORR3D)
FORCE:
endif

#-----------------------------------------------------------------------------
# fconn objects
#-----------------------------------------------------------------------------
$(FCONNDIR)/$(_OBJDIR)/fcmat_flt.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex fcmat_flt.o

$(FCONNDIR)/$(_OBJDIR)/matrix_flt.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex matrix_flt.o

$(FCONNDIR)/$(_OBJDIR)/edgestats_flt.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex edgestats_flt.o

$(FCONNDIR)/$(_OBJDIR)/nodedeg_flt.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex nodedeg_flt.o

$(FCONNDIR)/$(_OBJDIR)/fcmat_dbl.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex fcmat_dbl.o

$(FCONNDIR)/$(_OBJDIR)/matrix_dbl.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex matrix_dbl.o

$(FCONNDIR)/$(_OBJDIR)/edgestats_dbl.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex edgestats_dbl.o

$(FCONNDIR)/$(_OBJDIR)/nodedeg_dbl.o:
	$(MAKE) -C $(FCONNDIR)/src -f makefile-mex nodedeg_dbl.o
