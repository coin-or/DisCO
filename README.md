# DisCO [![Build Status](https://travis-ci.org/aykutbulut/DisCO.svg?branch=master)](https://travis-ci.org/aykutbulut/DisCO)

DisCO is distributed under the Eclipse Public License and is
freely redistributable. All source code and documentation is Copyright
2014-2016 by Lehigh University, Aykut Bulut and Ted Ralphs. This
README may be distributed freely.

DisCO shares some commit history with COIN-OR project [CHiPPS-Blis][8]. DisCO
version 0.1 has this heritage. Current master is re-written from scratch and
does not have Blis commits.

## 1. WHAT IS DisCO? ##

DisCO (Discrete Conic Optimization) is a solver for Mixed Integer Second Order
Conic Optimization (MISOCO) problems. It is developed on top of COIN-OR
High-Performance Parallel Search ([CHiPPS][8]) framework.

DisCO does branch and bound search to solve MISOCO problems. DisCO depends on
many other projects. It depends [OsiConic][1] on communicating with its
relaxation solvers. It depends on [CglConic][9] to cut infeasible solutions.

## 2. INSTALLATION ##

### 2.1 Basic Installation ###

s if they are not installed in your system. The easiest way of installing DisCO is using BuildTools fetch and build script. For this you can use the following commands in Linuxenvironment. After cloning DisCO, use
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
script. This can be acheived with the following command.  ```shell ./configure
--with-soco-solver=ipopt ```

Afterward you can call ```make install`` to build and install DisCO.

### 2.2 For advanced users ###

Make sure all dependencies are accessible through pkg-config. Then DisCO's
configure script will find them through pkg-config. Alternatively DisCO
configure script can locate other projects if ```--prefix``` configure flag is
set right. Assume other projects are installed at ```install_dir```. Then use

```shell
./configure --prefix=install_dir && make install
```

### 2.3 Specifying an Algorithm/Solver ###

DisCO implements an Outer Approximation algorithm and it is the default
behavior you will get. If you want to use DisCO with a typical branch and bound
algorithm (only integrality constraints are relaxed in nodes and corresponding
problems are solved with a SOCO solver) you need to specify this during
configure. DisCO depends on [OsiConic][2] in communicating with its
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

## 3. Current Testing Status ##

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

## 4. Documentation ##

You can refer to documentations of the dependant projects. DisCO uses doxygen for documentation purposes. ```make doxygen``` should produce a documentation of DisCO.

## 5. SUPPORT ##

### 5.1 Authors ###

Aykut Bulut (aykutblt@gmail.com)
Ted Ralphs (tkralphs@lehigh.edu)

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
