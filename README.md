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
Please see the INSTALL file for a guide to installation.


CURRENT TESTING STATUS
======================
1. Configurations
   - Cola: well tested.
   - Mosek: Does not terminate on some instances with some specific branching
strategy.
   - Cplex: Missing functions in the interface.
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
