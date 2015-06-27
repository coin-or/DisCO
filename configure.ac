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

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Disco],[trunk],[https://github.com/aykutbulut/DisCO/issues/new])

AC_COPYRIGHT([
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
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/DcoModel.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(Disco)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Disco)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# We filter out the -pedantic-errors flag for C and C++ compiler in
# case the user wants to use MPI, since otherwise the compilation
# fails
if test "${with_mpi_lib+set}" = set; then
  CFLAGS=`echo $CFLAGS | sed -e s/-pedantic-errors//g`
  CXXFLAGS=`echo $CXXFLAGS | sed -e s/-pedantic-errors//g`
fi

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_CHECK_PACKAGE(CoinDepend, [alps bcps cgl cglconic clp coinutils osi-clp osimosek osicplex cola osiconic], [DiscoLib])
if test $coin_has_coindepend != yes ; then
  AC_MSG_ERROR([Required packages not available.])
fi

AC_COIN_CHECK_PACKAGE(Sample, [coindatasample])

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In this case, the source code is taken from the DisCO examples directory

AC_CONFIG_LINKS(test/disco.par:src/DcoMain.cpp)

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

# ToDo: don't know why didn't work.
#AC_COIN_VPATH_LINK(examples/DcoMain.cpp)
#AC_COIN_VPATH_LINK(examples/disco.par)
#AC_COIN_VPATH_LINK(examples/flugpl.mps)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi Clp Cgl Alps Bcps CglConic Cola OsiConic OsiCplex OsiMosek)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
		 examples/VRP/Makefile
                 src/Makefile
		 src/disco.par
                 test/Makefile
		 disco.pc
		 disco-uninstalled.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_dco.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
