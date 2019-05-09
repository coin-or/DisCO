# DisCO
DisCO (Discrete Conic Optimization) is a solver for Mixed Integer Second Order
Conic Optimization (MISOCO) problems. It is developed on top of COIN-OR
High-Performance Parallel Search ([CHiPPS][8]) framework.

DisCO implements various branch-and-cut algorithms to solve MISOCO
problems. DisCO depends on many other projects. It depends [OsiConic][1] on
communicating with its relaxation solvers. It depends on [CglConic][9] to cut
infeasible solutions.

DisCO shares some commit history with COIN-OR project
[CHiPPS-Blis][13]. However, master branch is re-written from scratch.

DisCO is distributed under the Eclipse Public License and is
freely redistributable. All source code and documentation is Copyright
2001-2018 by Lehigh University, Aykut Bulut and Ted Ralphs. This
README may be distributed freely.

## Master Branch Build Status
[![Build Status](https://travis-ci.org/coin-or/DisCO.svg?branch=master)](https://travis-ci.org/coin-or/DisCO)

## Master Branch Code Quality Report
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a7c5abda08874549bee0324603ecaca1)](https://www.codacy.com/app/aykutbulut/DisCO?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=aykutbulut/DisCO&amp;utm_campaign=Badge_Grade)

## Download version 1.0
[ ![Download](https://api.bintray.com/packages/coin-or/download/DisCO/images/download.svg?version=1.0.0) ](https://bintray.com/coin-or/download/DisCO/1.0.0/link)

## Cite DisCO
[![DOI](https://zenodo.org/badge/35625228.svg)](https://zenodo.org/badge/latestdoi/35625228)

## Basic Installation

### Building on Linux

Most Linux distributions come with all the required tools installed. To obtain
the source code, the first step is to get the installer that will then
fetch the source for ALPS and all its dependencies. *You do not need to
clone this repository first, just do the following!* Open a terminal and execute
```
git clone https://www.github.com/coin-or/coinbrew
```
Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `coinbrew`
subdirectory. To execute the script, you can now do
```
cd coinbrew
./coinbrew
```
for the interactive version. Once you run the script in interactive mode,
you will be prompted to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.
```
./coinbrew fetch build DisCO --no-prompt
```
Options that would have been passed to the `configure` script under the old
build system can simply be added to the command-line. For example, to build
with debugging symbols, do
```
./coinbrew fetch build DisCO --no-prompt --enable-debug
```
To get help with additional options available in running the script, do
```
./coinbrew --help
```
The above procedures will build all required dependencies and DisCO itself.
Afterwards, the binaries will be installed in the directory `build/bin`
and the libraries in the directory `build/lib`. If you wish to
install in a different directory, such as `/usr/local`, then run the commands
```
./coinbrew fetch build DisCO --no-prompt --prefix=/path/to/install/dir
./coinbrew install DisCO --no-prompt
```
After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `LD_LIBRARY_PATH` if you want to link to COIN libraries. 

### Building on Windows (MSys2/CYGWIN and MinGW/MSVC)

By far, the easiest way to build on Windows is with the GNU autotools and the
GCC compilers. The first step is to install either
   * [Msys2](https://msys2.github.io/)
   * [CYGWIN](http://cygwin.org/)
   * [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

If you don't already have CYGWIN installed and don't want to fool around with
WSL (which is a great option if you already know your way around Unix), it is
recommended to use MSys2, since it provides a minimal toolset that is easy to
install. To get MSys2, either download the installer
[here](https://msys2.github.io/) or download and unzip MSys2 base from
[here](http://kent.dl.sourceforge.net/project/msys2/Base/x86_64/msys2-base-x86_64-20190512.tar.xz) 
(this is an out-of-date version, there may be a better place to get an archive
version). 

Following any of the above steps, you should have the `bash` command
(with Msys2, be sure to run `msys2_shell.bat` 
or manually add `msys64\usr\bin`, `msys64\mingw32\bin`, and
`msys64\mingw64\bin` to your Windows path).   

Once you have bash installed and in your `PATH`, open a Windows terminal and
type 
```
bash
pacman -S make wget tar patch dos2unix diffutils git svn
```
Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `coinbrew`
subdirectory. To execute the script, you can now do
```
cd coinbrew
./coinbrew
```
for the interactive version. Once you run the script in interactive mode,
you will be prompted to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.
```
./coinbrew fetch build DisCO --no-prompt
```
Options that would have been passed to the `configure` script under the old
build system can simply be added to the command-line. For example, to build
with debugging symbols, do
```
./coinbrew fetch build DisCO --no-prompt --enable-debug
```
To get help with additional options available in running the script, do
```
./coinbrew --help
```

To use the resulting binaries and/or libraries, you will need to add the
full path of the directory `build\bin` to your Windows executable
search `PATH`, or, alternatively, copy the conents of the build directory to 
`C:\Program Files (x86)\DisCO` and add the directory
`C:\Program Files (x86)\DisCO\bin` 
to your Windows executable search `PATH`. You may also consider adding
`C:\Program Files (x86)\DisCO\lib` to the `LIB` path and 
`C:\Program Files (x86)\DisCO\include` to the `INCLUDE` path. 

It is possible to use almost the exact same commands to build with the Visual
Studio compilers. Before doing any of the above commands in the Windows
terminal, first run the `vcvarsall.bat` script for your version of Visual
Studio. Note that you will also need a compatible Fortran compiler if you want
to build any projects requiring Fortran (`ifort` is recommended, but not
free). Then follow all the steps above, but replace the `build` command
with

```
./coinbrew fetch build DisCO --no-prompt --enable-msvc
```

### Building on OS X

OS X is a Unix-based OS and ships with many of the basic components needed to
build COIN-OR, but it's missing some things. For examples, the latest versions
of OS X come with the `clang` compiler but no Fortran compiler. You may also
be missing the `wget` utility and `subversion` and `git` clients (needed for
obtaining source code). The easiest way to get these missing utilitites is to
install Homebrew (see http://brew.sh). After installation, open a terminal and
do

```
brew install gcc wget svn git
```

To obtain the source code, the first step is to get the installer that will
then fetch the source for DisCO and all its dependencies. *You do not need to
clone DisCO first, just do the following!* Open a terminal and execute
```
git clone https://www.github.com/coin-or/coinbrew
```
Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `coinbrew`
subdirectory. To execute the script in interactive mode, do

```
cd coinbrew
./coinbrew
```
Once you run the script,you will be prompted interactively to select a project to 
fetch and build. the rest should happen automagically. Alternatively, the following 
command-line incantation will execute the procedure non-interactively.
```
./coinbrew fetch build DisCO --no-prompt
```
With this setup, `clang` will be used for compiling C++ by default and
`gfortran` will be used for Fortran. Since `clang` uses the GNU standard
library, `gfortran` is compatible.

If you want to use the `gcc` compiler provided by Homebrew, then replace the
`build` command above with
```
./coinbrew build DisCO --no-prompt CC=gcc-x CXX=g++-x
```
where `x` is the version of `gcc` you have installed.
Options that would have been passed to the `configure` script under the old
build system can simply be added to the command-line. For example, to build
with debugging symbols, do
```
./coinbrew fetch build DisCO --no-prompt --enable-debug
```
To get help with additional options available in running the script, do
```
./coinbrew --help
```
If you wish to install in a different directory, such as `/usr/local`, then run
the commands
```
./coinbrew fetch build DisCO --no-prompt --prefix=/path/to/install/dir
./coinbrew install DisCO --no-prompt
```
After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `DYLD_LIBRARY_PATH` if you want to link to COIN libraries. 

### Specifying Algorithms and Solvers

#### Specifying Algorithms and Solvers for DisCO

DisCO implements an OA algorithm and it is the default behavior you will
get. If you want to use DisCO with a typical branch and bound algorithm (only
integrality constraints are relaxed in nodes and corresponding problems are
solved with a SOCO solver) you need to specify this during configure. DisCO
depends on [OsiConic][1] in communicating with its solver. There are three
solvers available, [Ipopt][2], [Mosek][3] and [Cplex][4]. [OsiIpopt][5],
[OsiMosek][6] and [OsiCplex][7] implement OsiConic interface for the
corresponding solvers.

To compile DisCO with Mosek/Cplex you should first compile [OSI][10] with
Mosek/Cplex. Check [this][11] to compile OSI with these commercial
solvers. Then, you should compile OsiMosek/OsiCplex. You should clone these
from their repo separately and run ```configure``` and ```make```.  Configure
script will find all dependencies if the ```--prefix``` is set right. Once
OsiMosek/OsiCplex is compiled and installed, You can configure DisCO as
follows.

```shell
./coinbrew build DisCO --no-prompt --with-soco-solver=mosek
```

To use Cplex, just replace ```mosek``` with ```cplex``` in the command
above. Similarly for Ipopt just use ```ipopt```.

#### Specifying SOCO Solver for CglConic

DisCO depends on [CglConic][9] to warm start OA algorithm and generate cuts
during branch and bound. CglConic solves the SOCO relaxation in the root node
to warm start OA algorithm. CglConic needs an SOCO solver for this.  You can
specify the solver to be used. To specify this solver, you need to
give ```--with-ipm-solver``` flag to configure script. For example, following
command configures CglConic with Mosek.

```shell
./coinbrew build DisCO --no-prompt --with-ipm-solver=mosek
```

This command assumes you compiled and installed OsiMosek. If not, check
[here](#specifying-algorithms-and-solvers-for-disco) for instructions.
Similarly you can use ```cplex```, ```ipopt``` or ```cola``` instead
of ```mosek```. If no SOCO solver is specified
(no ```--with-ipm-solver=xxxxx``` flag is given) CGL-Conic will use Ipopt.

### Parallel Build Using MPI

To use DisCO in parallel, you need to compile DisCO and [CHiPPS][8] with an MPI
library. The parallel build and MPI library to use are specified during the
configure step. The easiest way to acheive this is using the `coinbrew` script. 

```shell
./coinbrew build DisCO --no-prompt --enable-static \
--disable-shared --with-mpi-incdir=/usr/local/mpich/include \
--with-mpi-lib='-L/usr/local/mpich/lib  -lmpich' MPICC=mpicc MPICXX=mpiCC
```

This command builds DisCO and CHiPPS for parallel runs with the given MPI
options. This command compiles DisCO and its dependents in parallel with 3
processors. DisCO and all its dependents are installed
to ```/install_directory``` and intermediate files generated during build are
stored in ```/build_directory```. The command will pass Mosek options and
flags (```--with-mosek-incdir``` and ```--with-mosek-lib```) to OSI and OSI's
Mosek interface is compiled. Note that, if you want to use Mosek, you need
to clone and compile OsiMosek separately after running this command. See
[here](#specifying-algorithms-and-solvers). You should update the
options in the command above acording to your system and MPI
implementation you intend to use. DisCO is tested with both mpich2 and
openmpi.

### Compiling Individual Projects Independently

All dependent projects and DisCO can be compiled separately if done in the
right order. Moreover, you do not need to compile dependent projects if they
are already installed somewhere in your system.

Make sure all dependencies are accessible through pkg-config. Then DisCO's
configure script will find them using pkg-config. Alternatively DisCO
configure script can locate other projects if ```--prefix``` configure flag is
set right (points to the install directory of dependent projects). Assume
other projects are installed at ```/install_dir```. Then following command will
find them and configure DisCO accordingly.

```shell
./configure --prefix=/install_dir && make install
```

### Build Examples

#### Use Cplex or Mosek in CglConic.

CglConic will use Mosek to warm start OA method.

```shell
./coinbrew build DisCO --with-mosek-incdir=/home/aykut/opt/mosek/8/tools/platform/linux64x86/h --with-mosek-lib="-L/home/aykut/opt/mosek/8/tools/platform/linux64x86/bin -lmosek64"
```

#### Build with OA and MPI

```shell
./coinbrew build DisCO --prefix=/home/aykut/research/conic/software/disco/build_openmpi_opt MPICC=/usr/lib64/openmpi/bin/mpicc MPICXX=/usr/lib64/openmpi/bin/mpic++ --with-mpi-lib=/usr/lib64/openmpi/lib/libmpi.so --with-mpi-incdir=/usr/include/openmpi-x86_64 --with-mosek-incdir=/home/aykut/opt/mosek/8/tools/platform/linux64x86/h --with-mosek-lib="-L/home/aykut/opt/mosek/8/tools/platform/linux64x86/bin -lmosek64"
```

This will compile OSI with Mosek and CHiPPS with OpenMPI. You need to
reconfigure CglConic if you want to use Mosek to solve root node SOCO
relaxation as explained above. By default CglConic will use Ipopt. Once
compiled you can use disco in parallel using mpirun as follows.

```shell
mpirun -n 8 build/bin/disco Alps_instance input.mps
```

## Using DisCO ##

DisCO can read problems in Mosek's extended MPS format (it can handle CSECTION
in mps files, see http://docs.mosek.com/7.1/capi/The_MPS_file_format.html) for
SOCO problems. Once you compiled DisCO you can use is as follows.

```shell
path_to_disco/disco input.mps
```

This uses the default parameters (cut generation/branching/search etc.). You
can modify its default behavior by specifying parameters. See src/disco.par.in
file for available parameters. For example following solves input problem using
strong branching with generating gomory cuts (available in OA algorithm only)
in the root node.

```shell
path_to_disco/disco input.mps Alps_instance input.mps Dco_branchStrategy 3 Dco_cutGomoryStrategy 1
```

## Current Testing Status ##

Following is a list of DisCO options tested. Check [Aykut's dissertation][12]
for more computational results of DisCO with various algorithms and solvers.

  * Operating Systems
    - Linux: Well tested.
    - Mac OS: I did not test DisCO in OS. In theory this should work.
    - Windows: Not tested.
  * Algorithms
    - OA: Well tested and works fine.
    - Ipopt: Works fine. There are missing functions in the [OsiIpopt][5] interface.
      Tested on CBLIB 2014 and random problems. Ipopt fails to converge on some
      instances.  We beleive this is due to nonsmooth formulation of the conic
      constraints in [OsiIpopt][5] interface.
    - Cola: Well tested, works fine.
    - Mosek: Works fine. There are missing functions in the interface. It is
      complete enough to work with DisCO. Well tested on CBLIB 2014 and random
      problems. Mosek might fail on numerically challanging instances.
    - Cplex: Missing functions in the interface. Interface is complete enough
      to work with DisCO. Tested on CBLIB 2014 and random problems. Cplex
      rarely fails on some instances.
  * Branching/Cutting
    - When OA algorithm is used and Ipopt is chosen as an IPM solver in
      CglConic, Ipopt might fail, on some problems, at the root node. You can
      use Mosek or Cplex for this if it is available to you.
  * MPI testing,
    - MPICH2: Tested and works fine.
    - OpenMPI: Tested up to 128 processors and works fine. Great parallelization
      performance when the tree is well balanced.

## Documentation

You can refer to documentations of the dependent projects. DisCO uses doxygen
for documentation purposes. Running ```make doxydoc``` in ```Disco``` directory
inside the directory you specified in option ```--build-dir``` will produce
doxygen documentation of DisCO.

Check [Aykut's dissertation][12] for detailed explanations of the algorithms
implemented.


## SUPPORT

### Authors

Aykut Bulut (aykutblt@gmail.com), Ted Ralphs (ted@lehigh.edu).

### Bug Reports

Bug reports should be reported at the DisCO repository at
https://github.com/coin-or/DisCO/issues/new

[1]: https://github.com/aykutbulut/OSI-CONIC
[2]: https://github.com/coin-or/Ipopt
[3]: https://mosek.com/
[4]: https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
[5]: https://github.com/aykutbulut/OsiIpopt
[6]: https://github.com/aykutbulut/OSI-MOSEK
[7]: https://github.com/aykutbulut/OsiCplex
[8]: https://github.com/coin-or/CHiPPS-ALPS
[9]: https://github.com/aykutbulut/CGL-CONIC
[10]: https://github.com/coin-or/Osi
[11]: https://github.com/coin-or/Osi#dynamically-loading-commercial-solver-libraries
[12]: https://preserve.lehigh.edu/etd/2981/
[13]: https://github.com/coin-or/CHiPPS-BLIS
