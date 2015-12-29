#==========================================================================#
# This file is part of the Discrete Conic Optimization (DisCO) solver.     #
#                                                                          #
# DisCO is distributed under the Eclipse Public License.                   #
#                                                                          #
# Authors:                                                                 #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#          Aykut Bulut, Lehigh University                                  #
#                                                                          #
# Conceptual Design:                                                       #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#          Laszlo Ladanyi, IBM T.J. Watson Research Center                 #
#          Matthew Saltzman, Clemson University                            #
#          Aykut Bulut, Lehigh University                                  #
#                                                                          #
#                                                                          #
# Copyright (C) 2001-2015, Lehigh University, Yan Xu, Ted Ralphs and       #
# Aykut Bulut. All Rights Reserved.                                        #
#==========================================================================#

AUTOMAKE_OPTIONS = foreign

########################################################################
#                           Common Flags                               #
########################################################################

LDADD = $(DISCOLIB_LIBS)
DEPENDENCIES = $(DISCOLIB_DEPENDENCIES)
AM_CPPFLAGS = $(DISCOLIB_CFLAGS)


if OA_SOCO_SOLVER
  AM_CPPFLAGS += -D__OA__
endif

if COLA_SOCO_SOLVER
  AM_CPPFLAGS += -D__COLA__
endif

if MOSEK_SOCO_SOLVER
  AM_CPPFLAGS += -D__OSI_MOSEK__
endif

if CPLEX_SOCO_SOLVER
  AM_CPPFLAGS += -D__OSI_CPLEX__
endif

########################################################################
#                               libDisco                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libDisco.la

# List all source files for this library, including headers
libDisco_la_SOURCES = \
	DcoBranchObjectInt.cpp \
	DcoBranchObjectInt.hpp \
	DcoBranchObjectBilevel.cpp \
	DcoBranchObjectBilevel.hpp \
	DcoBranchStrategyMaxInf.cpp \
	DcoBranchStrategyMaxInf.hpp \
	DcoBranchStrategyPseudo.cpp \
	DcoBranchStrategyPseudo.hpp \
	DcoBranchStrategyRel.cpp \
	DcoBranchStrategyRel.hpp \
	DcoBranchStrategyStrong.cpp \
	DcoBranchStrategyStrong.hpp \
	DcoBranchStrategyBilevel.cpp \
	DcoBranchStrategyBilevel.hpp \
	DcoConGeneratorBase.hpp \
	DcoConGeneratorBase.cpp \
	DcoConGenerator.cpp \
	DcoConGenerator.hpp \
	DcoConicConGenerator.cpp \
	DcoConicConGenerator.hpp \
	DcoConstraint.cpp \
	DcoConstraint.hpp \
	DcoConicConstraint.cpp \
	DcoConicConstraint.hpp \
	Dco.hpp \
	DcoHelp.cpp \
	DcoHelp.hpp \
	DcoHeuristic.hpp \
	DcoHeurRound.cpp \
	DcoHeurRound.hpp \
	DcoLicense.hpp \
	DcoMessage.cpp \
	DcoMessage.hpp \
	DcoModel.cpp \
	DcoModel.hpp \
	DcoNodeDesc.hpp \
	DcoObjectInt.cpp \
	DcoObjectInt.hpp \
	DcoParams.cpp \
	DcoParams.hpp \
	DcoPresolve.hpp \
	DcoPresolve.cpp \
	DcoPseudo.cpp \
	DcoPseudo.hpp \
	DcoSolution.hpp \
	DcoSubTree.hpp \
	DcoTreeNode.cpp \
	DcoTreeNode.hpp \
	DcoVariable.hpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libDisco_la_LIBADD = $(DISCOLIB_LIBS)
endif

# This is for libtool (on Windows)
libDisco_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                        disco binary                                   #
########################################################################

bin_PROGRAMS = disco
disco_SOURCES = DcoMain.cpp
disco_LDADD = libDisco.la $(DISCOLIB_LIBS)
disco_DEPENDENCIES = libDisco.la $(DISCOLIB_DEPENDENCIES)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# AM_CPPFLAGS = $(DISCOLIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	DcoBranchObjectInt.hpp \
	DcoBranchStrategyMaxInf.hpp \
	DcoBranchStrategyPseudo.hpp \
	DcoBranchStrategyRel.hpp \
	DcoBranchStrategyStrong.hpp \
	DcoConGeneratorBase.hpp \
	DcoConGenerator.hpp \
	DcoConicConGenerator.hpp \
	DcoConstraint.hpp \
	DcoConicConstraint.hpp \
	Dco.hpp \
	DcoHelp.hpp \
	DcoHeuristic.hpp \
	DcoHeurRound.hpp \
	DcoLicense.hpp \
	DcoMessage.hpp \
	DcoModel.hpp \
	DcoNodeDesc.hpp \
	DcoObjectInt.hpp \
	DcoParams.hpp \
	DcoPresolve.hpp \
	DcoPseudo.hpp \
	DcoSolution.hpp \
	DcoSubTree.hpp \
	DcoTreeNode.hpp \
	DcoVariable.hpp

########################################################################
#                     Installing manifest (MSVC++)                     #
########################################################################

if COIN_CXX_IS_CL
install-exec-hook:
	test -s disco.exe.manifest && \
	mt -manifest disco.exe.manifest -outputresource:disco.exe;\
	cp disco.exe $(bindir)/disco.exe
endif

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

install-exec-local:
	$(install_sh_DATA) config_dco.h $(DESTDIR)$(includecoindir)/DcoConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/DcoConfig.h
