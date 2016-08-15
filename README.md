# DisCO [![Build Status](https://travis-ci.org/aykutbulut/DisCO.svg?branch=master)](https://travis-ci.org/aykutbulut/DisCO)

DisCO (Discrete Conic Optimization) is a solver for Mixed Integer Second Order
Conic Optimization (MISOCO) problems. It is developed on top of COIN-OR
High-Performance Parallel Search ([CHiPPS][8]) framework.

DisCO does branch and bound search to solve MISOCO problems. DisCO depends on
many other projects. It depends [OsiConic][1] on communicating with its
relaxation solvers. It depends on [CglConic][9] to cut infeasible solutions.

DisCO shares some commit history with COIN-OR project [CHiPPS-Blis][8]. DisCO
version 0.1 has this heritage. Current master is re-written from scratch and
does not have Blis commits.

DisCO is distributed under the Eclipse Public License and is
freely redistributable. All source code and documentation is Copyright
2014-2016 by Lehigh University, Aykut Bulut and Ted Ralphs. This
README may be distributed freely.

## 1. Installation ##

### 1.1 Basic Installation ###

DisCO depends on many other projects. You should compile the dependant projects
if they are not installed in your system. The easiest way of installing DisCO
is using BuildTools fetch and build script. For this you can use the following
commands in Linuxenvironment. After cloning DisCO, use

```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch > /dev/null
bash BuildTools/get.dependencies.sh build
```

This compiles DisCO with Outer Approximation (OA) algorithm. This algorithm
relaxes integrality constraints and conic constraints. Performs a branch and
bound search to find a solution that satisfy both of these constraints.

There are other algorithms implemented in DisCO. You can use a typical branch
and bound algorithm where at each node only integrality constraints are
relaxed. For this you need to provide a Second Order Conic Optimization (SOCO)
solver. For now DisCO supports 3 solvers, Ipopt, Mosek and Cplex. To use DisCO
with Ipopt you need to add ```--with-soco-solver=ipopt``` flag to configure
script. This can be acheived with the following command.
```shell
./configure --with-soco-solver=ipopt
```

Afterward you can call ```make install``` to build and install DisCO.

### 1.2 For advanced users ###

Make sure all dependencies are accessible through pkg-config. Then DisCO's
configure script will find them through pkg-config. Alternatively DisCO
configure script can locate other projects if ```--prefix``` configure flag is
set right. Assume other projects are installed at ```install_dir```. Then use

```shell
./configure --prefix=install_dir && make install
```

### 1.3 Specifying an Algorithm/Solver ###

DisCO implements an Outer Approximation algorithm and it is the default
behavior you will get. If you want to use DisCO with a typical branch and bound
algorithm (only integrality constraints are relaxed in nodes and corresponding
problems are solved with a SOCO solver) you need to specify this during
configure. DisCO depends on [OsiConic][1] in communicating with its
solver. There are three solvers available, [Ipopt][2], [Mosek][3] and
[Cplex][4]. [OsiIpopt][5], [OsiMosek][6] and [OsiCplex][7] implement OsiConic
interface for the corresponding solvers.

To compile DisCO with Mosek/Cplex you should first compile OSI with
Mosek/Cplex. Then you should compile OsiMosek/OsiCplex. Please check OSI and
OsiMosek/OsiCplex documentations for details.

Once OsiMosek is compiled and installed, You can configure DisCO as follows.

```shell
./configure --with-soco-solver=mosek
```

For cplex just replace ```mosek``` with ```cplex```. Similarly for Ipopt just
use ```ipopt```.

To specify the solver, you need to give ```--with-ipm-solver``` flag to
configure script. For example, following command configures OSI-Conic with
Mosek.

```shell
./configure --with-ipm-solver=mosek
```

Similarly you can use ```cplex```, ```ipopt``` or ```cola``` instead of
```mosek```. If no IPM solver is specified CGL-Conic will use Ipopt.


### 1.4 Compiling with MPI ###

To compile with MPI you need to give the following options to configure
for DisCO and dependands projects.

```shell
./configure --disable-shared --enable-static --with-mpi-lib=/usr/lib/libmpich.so.3.2 --with-mpi-incdir=/usr/lib/mpich/include MPICC=mpicc.mpich MPICXX=mpic++.mpich"
```

You should update directory locations and executable names acording to your
system and MPI implementation you intend to use. DisCO is tested with and works for
both mpich2 and openmpi.

## 2. Current Testing Status ##

  * Operating Systems
    - Linux: Well tested.
    - OS: I did not test DisCO in OS. In theory this should be possible.
  * Algorithms
    - OA: Well tested and works fine.
    - Ipopt: Works fine. There are missing functions in the interface. Needs
      testing.
    - Cola: Well tested, works fine.
    - Mosek: Works fine. There are missing functions in the interface. It is
      complete enough to work with DisCO. Needs extensive testing.
    - Cplex: Missing functions in the interface. Interface is complete enough
      to work with DisCO. Needs extensive testing.
  * Branching/Cutting
    - When OA algorithm is used and Ipopt is chosen as an IPM solver in
      CglConic, Ipopt might fail on solving at the root node for approximation
      purposes.
  * MPI testing,
    - MPICH2: Hangs (or seems hanging on some instances of CBLIB). Works fine
      in most of the instances. Needs more testing for performance assesment.
    - OpenMPI: Hangs (or seems hanging on some instances of CBLIB). Works fine
      in most of the instances. Needs more testing for performance assesment.


## 3. Using DisCO ##

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

## 4. Documentation ##

You can refer to documentations of the dependant projects. DisCO uses doxygen
for documentation purposes. ```make doxygen``` will produce a documentation
of DisCO.

## 5. SUPPORT ##

### 5.1 Authors ###

Aykut Bulut (aykutblt@gmail.com), Ted Ralphs (tkralphs@lehigh.edu).

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
