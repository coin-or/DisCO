#==========================================================================#
# This file is part of the BiCePS Linear Integer Solver (BLIS).            #
#                                                                          #
# BLIS is distributed under the Common Public License as part of the       #
# COIN-OR repository (http://www.coin-or.org).                             #
#                                                                          #
# Authors:                                                                 #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#                                                                          #
# Conceptual Design:                                                       #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#          Laszlo Ladanyi, IBM T.J. Watson Research Center                 #
#          Matthew Saltzman, Clemson University                            #
#                                                                          #
#                                                                          #
# Copyright (C) 2001-2009, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13
#          Yan Xu                     Lehigh    2006-08-27

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libBlis                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libBlis.la

# List all source files for this library, including headers
libBlis_la_SOURCES = \
	BlisBranchObjectInt.cpp \
	BlisBranchObjectInt.h \
	BlisBranchObjectBilevel.cpp \
	BlisBranchObjectBilevel.h \
	BlisBranchStrategyMaxInf.cpp \
	BlisBranchStrategyMaxInf.h \
	BlisBranchStrategyPseudo.cpp \
	BlisBranchStrategyPseudo.h \
	BlisBranchStrategyRel.cpp \
	BlisBranchStrategyRel.h \
	BlisBranchStrategyStrong.cpp \
	BlisBranchStrategyStrong.h \
	BlisBranchStrategyBilevel.cpp \
	BlisBranchStrategyBilevel.h \
	BlisConGenerator.cpp \
	BlisConGenerator.h \
	BlisConstraint.cpp \
	BlisConstraint.h \
	Blis.h \
	BlisHelp.cpp \
	BlisHelp.h \
	BlisHeuristic.h \
	BlisHeurRound.cpp \
	BlisHeurRound.h \
	BlisLicense.h \
	BlisMessage.cpp \
	BlisMessage.h \
	BlisModel.cpp \
	BlisModel.h \
	BlisNodeDesc.h \
	BlisObjectInt.cpp \
	BlisObjectInt.h \
	BlisParams.cpp \
	BlisParams.h \
	BlisPresolve.h \
	BlisPresolve.cpp \
	BlisPseudo.cpp \
	BlisPseudo.h \
	BlisSolution.h \
	BlisSubTree.h \
	BlisTreeNode.cpp \
	BlisTreeNode.h \
	BlisVariable.h

# This is for libtool (on Windows)
libBlis_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                        blis binary                                   #
########################################################################

LDADD = $(BCPSOBJDIR)/src/libBcps.la \
	$(ALPSOBJDIR)/src/libAlps.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(OSIOBJDIR)/src/libOsi.la \
	$(CLPOBJDIR)/src/libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la 

bin_PROGRAMS = blis
blis_SOURCES = BlisMain.cpp 
blis_LDADD = libBlis.la $(LDADD)

LIBS += $(ADDLIBS) \
	`cat $(BCPSOBJDIR)/bcps_addlibs.txt` \
	`cat $(ALPSOBJDIR)/alps_addlibs.txt` \
	`cat $(CGLOBJDIR)/cgl_addlibs.txt` \
	`cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(CLPOBJDIR)/clp_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`

blis_LDFLAGS = -rpath $(libdir)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = \
        -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
        -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(ALPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(ALPSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BCPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BCPSOBJDIR)/inc`

# Include MPI header directory if MPI is used
if COIN_HAS_MPI
AM_CPPFLAGS +=  -I`$(CYGPATH_W) $(MPIINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../inc/config_blis.h \
	BlisBranchObjectInt.h \
	BlisBranchStrategyMaxInf.h \
	BlisBranchStrategyPseudo.h \
	BlisBranchStrategyRel.h \
	BlisBranchStrategyStrong.h \
	BlisConGenerator.h \
	BlisConstraint.h \
	Blis.h \
	BlisHelp.h \
	BlisHeuristic.h \
	BlisHeurRound.h \
	BlisLicense.h \
	BlisMessage.h \
	BlisModel.h \
	BlisNodeDesc.h \
	BlisObjectInt.h \
	BlisParams.h \
	BlisPresolve.h \
	BlisPseudo.h \
	BlisSolution.h \
	BlisSubTree.h \
	BlisTreeNode.h \
	BlisVariable.h

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = BlisConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
