#============================================================================#
# This file is part of the Discrete Conic Optimization (DisCO) Solver.       #
#                                                                            #
# DisCO is distributed under the Eclipse Public License as part of the       #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Authors:                                                                   #
#          Aykut Bulut, Lehigh University                                    #
#          Yan Xu, Lehigh University                                         #
#          Ted Ralphs, Lehigh University                                     #
#                                                                            #
# Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and       #
#                          Ted Ralphs.                                       #
# All Rights Reserved.                                                       #
#============================================================================#


AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	examples/VRP/AUTHORS \
	examples/VRP/INSTALL \
	examples/VRP/LICENSE \
	examples/VRP/Makefile.in \
	examples/VRP/README \
	examples/VRP/VrpCommonTypes.h \
	examples/VRP/VrpConstants.h \
	examples/VRP/VrpCutGenerator.cpp \
	examples/VRP/VrpCutGenerator.h \
	examples/VRP/VrpHeurTSP.cpp \
	examples/VRP/VrpHeurTSP.h \
	examples/VRP/VrpMain.cpp \
	examples/VRP/VrpModel.cpp \
	examples/VRP/VrpModel.h \
	examples/VRP/VrpNetwork.cpp \
	examples/VRP/VrpNetwork.h \
	examples/VRP/VrpParams.cpp \
	examples/VRP/VrpParams.h \
	examples/VRP/VrpSolution.h \
	examples/VRP/VrpVariable.h

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )


clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest doxydoc

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = disco.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/DisCO

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@:$(DESTDIR)$(pkgconfiglibdir) \
	$(PKG_CONFIG) --libs disco > $(addlibsdir)/disco_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libDisco.lib @DISCOLIB_LIBS_INSTALLED@" > $(addlibsdir)/disco_addlibs.txt
else
	echo -L@abs_lib_dir@ -lDisco @DISCOLIB_LIBS_INSTALLED@ > $(addlibsdir)/disco_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/disco_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
