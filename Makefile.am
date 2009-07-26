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

## $Id: Makefile.am 788 2006-06-01 18:57:08Z andreasw $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = doxydoc/doxygen.conf \
	Blis/src/Blis.h\
	Blis/src/BlisBranchObjectBilevel.cpp\
	Blis/src/BlisBranchObjectBilevel.h \
	Blis/src/BlisBranchObjectInt.cpp \
	Blis/src/BlisBranchObjectInt.h \
	Blis/src/BlisBranchStrategyBilevel.cpp \
	Blis/src/BlisBranchStrategyBilevel.h \
	Blis/src/BlisBranchStrategyMaxInf.cpp \
	Blis/src/BlisBranchStrategyMaxInf.h \
	Blis/src/BlisBranchStrategyPseudo.cpp \
	Blis/src/BlisBranchStrategyPseudo.h \
	Blis/src/BlisBranchStrategyRel.cpp \
	Blis/src/BlisBranchStrategyRel.h \
	Blis/src/BlisBranchStrategyStrong.cpp \
	Blis/src/BlisBranchStrategyStrong.h \
	Blis/src/BlisConfig.h \
	Blis/src/BlisConGenerator.cpp \
	Blis/src/BlisConGenerator.h \
	Blis/src/BlisConstraint.cpp \
	Blis/src/BlisConstraint.h \
	Blis/src/BlisHelp.cpp \
	Blis/src/BlisHelp.h \
	Blis/src/BlisHeuristic.h \
	Blis/src/BlisHeurRound.cpp \
	Blis/src/BlisHeurRound.h \
	Blis/src/BlisLicense.h \
	Blis/src/BlisMain.cpp \
	Blis/src/BlisMessage.cpp \
	Blis/src/BlisMessage.h \
	Blis/src/BlisModel.cpp \
	Blis/src/BlisModel.h \
	Blis/src/BlisNodeDesc.h \
	Blis/src/BlisObjectInt.cpp \
	Blis/src/BlisObjectInt.h \
	Blis/src/BlisParams.cpp \
	Blis/src/BlisParams.h \
	Blis/src/BlisPresolve.cpp \
	Blis/src/BlisPresolve.h \
	Blis/src/BlisPseudo.cpp \
	Blis/src/BlisPseudo.h \
	Blis/src/BlisSolution.h \
	Blis/src/BlisSubTree.h \
	Blis/src/BlisTreeNode.cpp \
	Blis/src/BlisTreeNode.h \
	Blis/src/BlisVariable.h

########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = $(subdirs)

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd Blis; $(MAKE) test

unitTest: test

tests:
	for dir in $(subdirs); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; $(MAKE) test) \
	  fi; \
	done

unitTests: tests

doxydoc:
	cd $(srcdir); doxygen doxydoc/doxygen.conf

.PHONY: test unitTest tests unitTests doxydoc

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
