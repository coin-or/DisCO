DisCO Version 0.1 README
=========================
DisCO is distributed under the Eclipse Public License and is
freely redistributable. All source code and documentation is Copyright
2001-2013 by Lehigh University, Aykut Bulut, Yan Xu, and Ted Ralphs. This
README may be distributed freely.

WHAT IS DisCO?
=============
DisCO (Discrete Conic Optimization) is a solver for mixed-integer second order
conic optimization problems. It is developed on top of COIN-OR High-Performance
Parallel Search (CHiPPS) framework.


INSTALLATION
============
After cloning DisCO, use
```shell
sh get_dependencies.sh && sh compile.sh
```

For advanced users:
Make sure all dependencies are accessible through pkg-config. Then DisCO's configure script will find them through pkg-config. Alternatively DisCO configure script can locate other projects if --prefix configure flag is set right. Assume other projects are installed at install_dir. Then use
```shell
./configure --prefix=install_dir && make install
```


CURRENT TESTING STATUS
======================
1. Configurations
   - Cola: well tested.
   - Cplex: Missing functions in the interface.
   - OA: Fails on some instances. Check CBLIB problems robust and classical.
   - Dco_branchStrategy reliability is broken

SUPPORT
=======
1. Authors
Aykut Bulut (aykut@lehigh.edu)
Yan Xu (yax2@lehigh.edu)
Ted Ralphs (tkralphs@lehigh.edu), Project Manager
Laci Ladanyi (ladanyi@us.ibm.com)
Matt Saltzman (mjs@clemson.edu)

2. Bug Reports
Bug reports should be reported on the CHiPPS development web site at
https://github.com/aykutbulut/DisCO/issues/new
