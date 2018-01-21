# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linuxegcs

CXX           =
ObjSuf        = o
SrcSuf        = cxx
SrcSufB       = C
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 
DictSuf       = Dict.o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

ifeq ($(ARCH),linux)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxrh42)
# Linux with gcc 2.7.2.x (RedHat 4.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxdeb)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxdeb2)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxsuse6)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxegcs)
# Linux with egcs (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxkcc)
# Linux with the KAI compiler
CXX           = KCC
CXXFLAGS      = -fPIC +K0
LD            = KCC
LDFLAGS       = -O
SOFLAGS       =
endif

ifeq ($(ARCH),linuxppcegcs)
# MkLinux with egcs/glibc
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared -Wl,-soname,
endif

ifeq ($(ARCH),linuxia64gcc)
# Itanium Linux with gcc 2.9x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxia64sgi)
# Itanium Linux with sgiCC
CXX           = sgiCC
CXXFLAGS      = -O -Wall -fPIC
LD            = gsgiCC
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxalphaegcs)
# Alpha Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxarm)
# ARM Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------

TGVKellyO           =TGVKelly.$(ObjSuf) \
                       TGVKellyDict.$(ObjSuf)
TGVKellyS           =TGVKelly.$(SrcSuf) \
                       TGVKellyDict.$(SrcSuf)
TGVKellySO          =libTGVKelly.$(DllSuf)

TGenGeoO           =TGenGeo.$(ObjSuf) \
                       TGenGeoDict.$(ObjSuf)
TGenGeoS           =TGenGeo.$(SrcSuf) \
                       TGenGeoDict.$(SrcSuf)
TGenGeoSO          =libTGenGeo.$(DllSuf)

TGenSBSGeoO           =TGenSBSGeo.$(ObjSuf) \
                       TGenSBSGeoDict.$(ObjSuf)
TGenSBSGeoS           =TGenSBSGeo.$(SrcSuf) \
                       TGenSBSGeoDict.$(SrcSuf)
TGenSBSGeoSO          =libTGenSBSGeo.$(DllSuf)


TGenBaseO           =TGenBase.$(ObjSuf) \
                       TGenBaseDict.$(ObjSuf)
TGenBaseS           =TGenBase.$(SrcSuf) \
                       TGenBaseDict.$(SrcSuf)
TGenBaseSO          =libTGenBase.$(DllSuf)

TGenDVCSO           =TGenDVCS.$(ObjSuf) \
                      TGenDVCSDict.$(ObjSuf)
TGenDVCSS           =TGenDVCS.$(SrcSuf) \
                       TGenDVCSDict.$(SrcSuf)
TGenDVCSSO      =libTGenDVCS.$(DllSuf)

TGenPi0O           =TGenPi0.$(ObjSuf) \
                      TGenPi0Dict.$(ObjSuf)
TGenPi0S           =TGenPi0.$(SrcSuf) \
                       TGenPi0Dict.$(SrcSuf)
TGenPi0SO      =libTGenPi0.$(DllSuf)

MAINO               = SBSDVCS.$(ObjSuf)
MAINS               = SBSDVCS.$(SrcSufB)
MAIN                = SBSDVCS$(ExeSuf)

#MAINO               = SBSPi0.$(ObjSuf)
#MAINS               = SBSPi0.$(SrcSufB)
#MAIN                = SBSPi0$(ExeSuf)

OBJS          = $(TGVKellyO) $(TGenGeoO) $(TGenSBSGeoO) $(TGenBaseO) $(TGenDVCSO) $(TGenPi0O) 

SLIB          = $(TGVKellySO) $(TGenGeoSO) $(TGenSBSGeoSO) $(TGenBaseSO) $(TGenDVCSSO) $(TGenPi0SO)

PROGRAMS      = $(MAIN)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

#all:            $(PROGRAMS) doc
all:            $(PROGRAMS)

lib:            $(ROOTGLIBS) $(OBJS)

#doc:
#	root.exe -b -q make_html.C

%.$(ObjSuf):%.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<

%Dict.$(SrcSuf):%.h
	rootcint -f $@ -c $< 

$(MAIN): $(MAINO) $(OBJS)
	$(LD) $(CXXFLAGS) $(LIBS) $^ -o $@ 

lib%.$(DllSuf): %Dict.$(ObjSuf) %.$(ObjSuf)
	$(LD) $(ROOTGLIBS) $(SOFLAGS) $(CXXFLAGS) $^ -o $@

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(SLIB) $(OBJS) *Dict.* *.def *.o


























