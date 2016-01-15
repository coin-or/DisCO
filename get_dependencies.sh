#!/bin/bash

# get projects from svn repo of coin-or server
#svn co https://projects.coin-or.org/svn/BuildTools/stable/0.8 BuildTools
svn co https://projects.coin-or.org/svn/CoinUtils/stable/2.10/CoinUtils
svn co https://projects.coin-or.org/svn/Osi/trunk/Osi
svn co https://projects.coin-or.org/svn/Clp/stable/1.16/Clp
svn co https://projects.coin-or.org/svn/Cgl/stable/0.59/Cgl
svn co https://projects.coin-or.org/svn/CHiPPS/Alps/trunk/Alps
svn co https://projects.coin-or.org/svn/CHiPPS/Bcps/trunk/Bcps
svn co https://projects.coin-or.org/svn/Ipopt/releases/3.12.4/Ipopt Ipopt
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Blas/stable/1.4 ThirdParty/Blas
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Lapack/stable/1.5 ThirdParty/Lapack
# IPOPT dependencies
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/HSL/stable/1.5    ThirdParty/HSL
##cd ThirdParty/HSL && ./get.HSL && cd ../..
wait
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Metis/stable/1.3  ThirdParty/Metis
cd ThirdParty/Metis && ./get.Metis && cd ../..
wait
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Mumps/stable/1.5  ThirdParty/Mumps
cd ThirdParty/Mumps && ./get.Mumps && cd ../..
wait
# data
svn co https://projects.coin-or.org/svn/Data/Sample/stable/1.2 Data/Sample

# get projects from github
git clone https://github.com/aykutbulut/OSI-CONIC OsiConic
git clone https://github.com/aykutbulut/OSI-MOSEK OsiMosek
git clone https://github.com/aykutbulut/OsiCplex OsiCplex
git clone https://github.com/aykutbulut/OsiIpopt OsiIpopt
git clone https://github.com/aykutbulut/COLA Cola
git clone https://github.com/aykutbulut/CGL-CONIC CglConic
