## Master Branch Build Status
[![Build Status](https://travis-ci.org/coin-or/DisCO.svg?branch=master)](https://travis-ci.org/coin-or/DisCO)

## Master Branch Code Quality Report
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a7c5abda08874549bee0324603ecaca1)](https://www.codacy.com/app/aykutbulut/DisCO?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=aykutbulut/DisCO&amp;utm_campaign=Badge_Grade)

## Download version 0.95
[ ![Download](https://api.bintray.com/packages/coin-or/download/DisCO/images/download.svg?version=0.95) ](https://bintray.com/coin-or/download/DisCO/0.95/link)

## Cite DisCO
[![DOI](https://zenodo.org/badge/36100320.svg)](https://zenodo.org/badge/latestdoi/36100320)

# DisCO
DisCO (Discrete Conic Optimization) is a solver for Mixed Integer Second Order
Conic Optimization (MISOCO) problems. It is developed on top of COIN-OR
High-Performance Parallel Search ([CHiPPS][8]) framework.

DisCO implements various branch-and-cut algorithms to solve MISOCO
problems. DisCO depends on many other projects. It depends [OsiConic][1] on
communicating with its relaxation solvers. It depends on [CglConic][9] to cut
infeasible solutions.

DisCO shares some commit history with COIN-OR project
[CHiPPS-Blis][8]. However, master branch is re-written from scratch.

DisCO is distributed under the Eclipse Public License and is
freely redistributable. All source code and documentation is Copyright
2001-2018 by Lehigh University, Aykut Bulut and Ted Ralphs. This
README may be distributed freely.

## 1. Installation ##

### 1.1 Basic Installation ###

DisCO depends on many other projects. You should compile the dependent projects
if they are not installed in your system. The easiest way of installing DisCO
is using BuildTools' fetch and build script. For this you can use the following
commands in Linux environment. After cloning DisCO, use

```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch > /dev/null
bash BuildTools/get.dependencies.sh build
```

This compiles DisCO with default options. DisCO use Outer Approximation (OA)
algorithm by default. This algorithm relaxes integrality constraints and conic
constraints. Performs a branch and bound search to find a solution that satisfy
both of these constraints. Check [Aykut's dissertation][12] for a detailed
description of the OA algorithm implemented.

There are other algorithms implemented in DisCO. You can use a typical branch
and bound algorithm where at each node only integrality constraints are
relaxed. See next section for instructions on building DisCO with other
algorithms and solvers.

### 1.2 Specifying Algorithms and Solvers ###

#### 1.2.1 Specifying Algorithms and Solvers for DisCO

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
./configure --with-soco-solver=mosek
```

To use Cplex, just replace ```mosek``` with ```cplex``` in the command
above. Similarly for Ipopt just use ```ipopt```.

#### 1.2.2 Specifying SOCO Solver for CglConic

DisCO depends on [CglConic][9] to warm start OA algorithm and generate cuts
during branch and bound. CglConic solves the SOCO relaxation in the root node
to warm start OA algorithm. CglConic needs an SOCO solver for this.  You can
specify the solver to be used. To specify this solver, you need to
give ```--with-ipm-solver``` flag to configure script. For example, following
command configures CglConic with Mosek.

```shell
./configure --with-ipm-solver=mosek
```

This command assumes you compiled and installed OsiMosek. If not, check
[Section 1.2.1](#specifying-algorithms-and-solvers-for-disco) for instructions.
Similarly you can use ```cplex```, ```ipopt``` or ```cola``` instead
of ```mosek```. If no SOCO solver is specified
(no ```--with-ipm-solver=xxxxx``` flag is given) CGL-Conic will use Ipopt.

### 1.3 Parallel Build Using MPI ###

To use DisCO in parallel, you need to compile DisCO and [CHiPPS][8] with an MPI
library. Parallel build and MPI library to use are specified during the
configure step. Easiest way to acheive this is using BuildTools' fetch and
build script. You need to clone BuildTools first
(see [Section 1.1](#basic-installation)). All options to be used can be
specified as arguments to fetch and build script, they are passed to the
configure scripts of the idividual dependents. An example command is as
follows.

```shell
sh BuildTools/get.dependencies.sh build --prefix=/install_directory --with-mosek-incdir=/opt/mosek/8/tools/platform/linux64x86/h --with-mosek-lib="-L/opt/mosek/8/tools/platform/linux64x86/bin -lmosek64" --parallel-jobs=3 --build-dir=/build_directory
```

This command builds DisCO and CHiPPS for parallel runs with the given MPI
options. This command compiles DisCO and its dependents in parallel with 3
processors. DisCO and all its dependents are installed
to ```/install_directory``` and intermediate files generated during build are
stored in ```/build_directory```. The command will pass Mosek options and
flags (```--with-mosek-incdir``` and ```--with-mosek-lib```) to OSI and OSI's
Mosek interface is compiled. Note that, if you want to use Mosek, you need
to clone and compile OsiMosek separately after running this command. See
[Section 1.2](#specifying-algorithms-and-solvers). You should update the
options in the command above acording to your system and MPI
implementation you intend to use. DisCO is tested with both mpich2 and
openmpi.

### 1.4 Compiling Individual Projects Independently ###

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

### 1.6 Build Examples

#### 1.6.1 Build with Default Options.

```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch > /dev/null
bash BuildTools/get.dependencies.sh build
```

This will use OA algorithm. Ipopt will be used in CglConic to solve the SOCO
relaxation of the root node to warm start OA algorithm.

#### 1.6.2 Use Cplex or Mosek in CglConic.

CglConic will use Mosek to warm start OA method. You need to clone and compile
OsiMosek separately. After installing OsiMosek, re-configure. compile and
install DisCO.

```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch > /dev/null
bash BuildTools/get.dependencies.sh build --with-mosek-incdir=/home/aykut/opt/mosek/8/tools/platform/linux64x86/h --with-mosek-lib="-L/home/aykut/opt/mosek/8/tools/platform/linux64x86/bin -lmosek64"
# clone OsiMosek
# compile and install OsiMosek
# configure, compile and install DisCO again to use Mosek
```

This will build all dependencies and OSI with mosek interface. Next you need to download and compile OsiMosek. Use following for this.

```shell
git clone https://github.com/aykutbulut/OSI-MOSEK.git
cd OSI-MOSEK
./configure --prefix=use_install_dir_from_above
make install
```

After OsiMosek is compiled and installed, you need to re-configure CglConic to
use Mosek. Go to CglConic build directory (typically ./build/CglConic/) and run
following.

```shell
make distclean
../../CglConic/configure --with-ipm-solver=mosek
make install
```

Finally we need to re-link disco to the newly compiled CglConic libraries. You
can remove DisCO binary only for this, and re-compile DisCO. Go to DisCO build
directory (typically ./build/Disco) and run the following.

```shell
rm src/disco
rm ../bin/disco
make install
```

These steps are similar for Cplex.

#### 1.6.3 Build with OA and MPI

```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch > /dev/null
bash BuildTools/get.dependencies.sh build --prefix=/home/aykut/research/conic/software/disco/build_openmpi_opt MPICC=/usr/lib64/openmpi/bin/mpicc MPICXX=/usr/lib64/openmpi/bin/mpic++ --with-mpi-lib=/usr/lib64/openmpi/lib/libmpi.so --with-mpi-incdir=/usr/include/openmpi-x86_64 --with-mosek-incdir=/home/aykut/opt/mosek/8/tools/platform/linux64x86/h --with-mosek-lib="-L/home/aykut/opt/mosek/8/tools/platform/linux64x86/bin -lmosek64"
```

This will compile OSI with Mosek and CHiPPS with OpenMPI. You need to
reconfigure CglConic if you want to use Mosek to solve root node SOCO
relaxation as explained above. By default CglConic will use Ipopt. Once
compiled you can use disco in parallel using mpirun as follows.

```shell
mpirun -n 8 path_to_build_dir/bin/disco Alps_instance input.mps
```

#### 1.6.4 Use Cplex/Mosek at Nodes

For this compile OSI with Cplex/Mosek and compile OsiCplex/OsiMosek as
explained in [Section 1.6.2](#use-cplex-or-mosek-in-cglconic). Once you
compiled OsiCplex/OsiMosek you need to re-configure Disco to use
Cplex/Mosek. Go to DisCO build directory (typically ```./build/Disco```) and
use the following commands.

```shell
make distclean
../../configure --disable-dependency-tracking --prefix=your_build_dir --with-soco-solver=cplex
make install
```

You can compile with cplex/mosek and MPI similarly.

## 2. Using DisCO ##

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

## 3. Current Testing Status ##

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

## 4. Documentation ##

You can refer to documentations of the dependent projects. DisCO uses doxygen
for documentation purposes. Running ```make doxydoc``` in ```Disco``` directory
inside the directory you specified in option ```--build-dir``` will produce
doxygen documentation of DisCO.

Check [Aykut's dissertation][12] for detailed explanations of the algorithms
implemented.


## 5. SUPPORT ##

### 5.1 Authors ###

Aykut Bulut (aykutblt@gmail.com), Ted Ralphs (ted@lehigh.edu).

### 5.2 Bug Reports ###

Bug reports should be reported at the DisCO repository at
https://github.com/aykutbulut/DisCO/issues/new

[1]: https://github.com/aykutbulut/OSI-CONIC
[2]: https://projects.coin-or.org/Ipopt
[3]: https://mosek.com/
[4]: https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
[5]: https://github.com/aykutbulut/OsiIpopt
[6]: https://github.com/aykutbulut/OSI-MOSEK
[7]: https://github.com/aykutbulut/OsiCplex
[8]: https://projects.coin-or.org/CHiPPS
[9]: https://github.com/aykutbulut/CGL-CONIC
[10]: https://projects.coin-or.org/Osi
[11]: https://projects.coin-or.org/Osi#Dynamicallyloadingcommercialsolverlibraries
[12]: https://preserve.lehigh.edu/etd/2981/
