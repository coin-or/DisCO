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


#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

m4_define([gitversion], m4_esyscmd([sh -c "git describe --always | cut -d '=' -f 2 | head --bytes -1"]))
AC_INIT([DisCO],[gitversion],[https://github.com/aykutbulut/DisCO/issues/new])

AC_DEFINE([DISCO_GIT_VERSION], [gitversion], [DisCO git version])

AC_COPYRIGHT([
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
])

##AC_REVISION([m4_esyscmd([./configure.commit])])




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

AC_COIN_CHECK_PACKAGE(CoinDepend, [alps bcps cgl cglconic clp coinutils osi-clp osiconic], [DiscoLib])
if test $coin_has_coindepend != yes ; then
  AC_MSG_ERROR([Required packages not available.])
fi
AC_COIN_CHECK_PACKAGE(Cola, [cola],     [DiscoLib])
AC_COIN_CHECK_PACKAGE(Ipopt, [ipopt],     [DiscoLib])
AC_COIN_CHECK_PACKAGE(OsiIpopt, [osiipopt],     [DiscoLib])
AC_COIN_CHECK_PACKAGE(Sample, [coindatasample])

#############################################################################
#                              check MPI libraries                          #
#############################################################################

AC_COIN_CHECK_USER_LIBRARY([Mpi], [MPI], [mpi.h], [MPI_Irecv], [], [AlpsLib])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available
AC_COIN_CHECK_PACKAGE(Cplex,  [osicplex],  [DiscoLib])
AC_COIN_CHECK_PACKAGE(Mosek,  [osimosek],  [DiscoLib])

#############################################################################
#                      Check for the SOCO solver                            #
#############################################################################
AC_MSG_CHECKING([for DisCO default solver])
AC_ARG_WITH([soco-solver],
[AC_HELP_STRING([--with-soco-solver@<:@=socosolver@:>@],
                [specify the SOCO solver in small
                letters(default socosolver=oa)])],
                [dco_soco_solver=$withval],[dco_soco_solver=oa])

found_package=true
case $dco_soco_solver in
   oa)
    AC_MSG_RESULT(oa)
     if [test $coin_has_coindepend = unavailable || test $coin_has_coindepend = skipping]; \
     then
         found_package=false
     fi
   ;;
   cola)
    AC_MSG_RESULT(Cola)
     if [test $coin_has_cola = unavailable || test $coin_has_cola = skipping]; \
     then
         found_package=false
     fi
   ;;
   ipopt)
    AC_MSG_RESULT(Ipopt)
     if [test $coin_has_osiipopt = unavailable || test $coin_has_osiipopt = skipping]; \
     then
         found_package=false
     fi
   ;;
   mosek)
     AC_MSG_RESULT(Mosek)
     if test $coin_has_mosek = false; \
     then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cplex = false; \
     then
         found_package=false
     fi
   ;;
   *)
   AC_MSG_ERROR(Unknown SOCO solver!)
   ;;
esac

if test $found_package = false; then
  AC_MSG_ERROR(Couldn't find the $dco_soco_solver package!)
fi

AM_CONDITIONAL(OA_SOCO_SOLVER,[test $dco_soco_solver = oa &&
test $coin_has_coindepend = yes])
AM_CONDITIONAL(COLA_SOCO_SOLVER,[test $dco_soco_solver = cola &&
test $coin_has_cola = yes])
AM_CONDITIONAL(IPOPT_SOCO_SOLVER,[test $dco_soco_solver = ipopt &&
test $coin_has_osiipopt = yes])
AM_CONDITIONAL(MOSEK_SOCO_SOLVER,[test $dco_soco_solver = mosek &&
test $coin_has_mosek = yes])
AM_CONDITIONAL(CPLEX_SOCO_SOLVER,[test $dco_soco_solver = cplex &&
test $coin_has_cplex = yes])

AM_CONDITIONAL(MOSEK_EXIST,[test $coin_has_mosek = yes])
AM_CONDITIONAL(CPLEX_EXIST,[test $coin_has_cplex = yes])

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In this case, the source code is taken from the DisCO examples directory

AC_CONFIG_LINKS(test/DcoMain.cpp:src/DcoMain.cpp)

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
                 test/disco.par
                 test/Makefile
                 disco.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_dco.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
