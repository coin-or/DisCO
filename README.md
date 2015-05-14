# DisCO
Discrete Conic Optimization

Solves discrete conic optimization problems. Depends on other COIN-OR tools. DisCO usea OsiConicSolverInterface to communicate with the continuous solver. We have Mosek and COLA that implements this interface and DisCO can use both to solve relaxation problems.

## Install
To install "../configure && make && make install" should work. Uses pkg-config to find dependencies.

## Use
To solve problems call DisCO with mps file input. "./disco input.mps" should work. For now we support Mosek type mps files only.

