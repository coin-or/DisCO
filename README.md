DisCO [![Build Status](https://travis-ci.org/aykutbulut/DisCO.svg?branch=master)](https://travis-ci.org/aykutbulut/DisCO)
=========================
DisCO is distributed under the Eclipse Public License and is
freely redistributable. All source code and documentation is Copyright
2001-2013 by Lehigh University, Aykut Bulut, Yan Xu, and Ted Ralphs. This
README may be distributed freely.

## 1. WHAT IS DisCO? ##
DisCO (Discrete Conic Optimization) is a solver for mixed-integer second order
conic optimization problems. It is developed on top of COIN-OR High-Performance
Parallel Search (CHiPPS) framework.


## 2. INSTALLATION ##
### 2.1 Basic Installation ###
DisCO depends on many other projects. You should compile the dependant projects if they are not installed in your system. The easiest way of installing DisCO is using BuildTools fetch and build script. For this you can use the following commands in Linuxenvironment. After cloning DisCO, use
```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch > /dev/null
bash BuildTools/get.dependencies.sh build
```

This compiles DisCO with Outer Approximation (OA) algorithm. This algorithm relaxes integrality constraints and conic constraints. Performs a branch and bound search to find a solution that satisfy both of these constraints.

There are other algorithms implemented in DisCO. You can use a typical branch and bound algorithm where at each node only integrality constraints are relaxed. For this you need to provide a Second Order Conic Optimization (SOCO) solver. For now DisCO supports 3 solvers, Ipopt, Mosek and Cplex. To use DisCO with Ipopt you need to add ```--with-soco-solver=ipopt``` flag to configure script. This can be acheived with the following command.
```shell
./configure --with-soco-solver=ipopt
```

Afterward you can call ```make install`` to build and install DisCO.

### 2.2 For advanced users ###
Make sure all dependencies are accessible through pkg-config. Then DisCO's configure script will find them through pkg-config. Alternatively DisCO configure script can locate other projects if --prefix configure flag is set right. Assume other projects are installed at install_dir. Then use
```shell
./configure --prefix=install_dir && make install
```


## 3. CURRENT TESTING STATUS ##
   - OA: well tested and works fine.
   - Cola: well tested.
   - Cplex: Missing functions in the interface.
   - Ipopt: works fine for branching strategy 0 (Dco_branchStrategy 0). Branching strategies should be updated considering IPM solvers.

## 4. SUPPORT ##
### 4.1 Authors ###
Aykut Bulut (aykut@lehigh.edu)
Ted Ralphs (tkralphs@lehigh.edu), Project Manager

### 4.2 Bug Reports ###
You can create new issues using https://github.com/aykutbulut/DisCO/issues/new.
