#!/bin/bash
mkdir -p build
cd build
build_dir=$PWD
inc_dir=${build_dir%%/}/include
lib_dir=${build_dir%%/}/lib
pkg_dir=${lib_dir%%/}/pkgconfig
PKG_CONFIG_PATH=${pkg_dir}:$PKG_CONFIG_PATH
export CXXFLAGS="-std=c++11 -O2"
# configure and install Blas
mkdir -p ThirdParty/Blas
cd ThirdParty/Blas
../../../ThirdParty/Blas/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Lapack
mkdir -p Lapack
cd Lapack
../../../ThirdParty/Lapack/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install ASL
mkdir -p ASL
cd ASL
../../../ThirdParty/ASL/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install HSL
mkdir -p HSL
cd HSL
../../../ThirdParty/HSL/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Metis
mkdir -p Metis
cd Metis
../../../ThirdParty/Metis/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Mumps
mkdir -p Mumps
cd Mumps
../../../ThirdParty/Mumps/configure --prefix=$build_dir
make -j 10 install
cd ../..
# configure and install CoinUtils
mkdir -p CoinUtils
cd CoinUtils
../../CoinUtils/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi
mkdir -p Osi
cd Osi
../../Osi/configure --prefix=$build_dir --with-mosek-incdir=/usr/local/mosek/7.1/tools/platform/linux64x86/h --with-mosek-lib="-L/usr/local/mosek/7.1/tools/platform/linux64x86/bin -lmosek64 -lmosekxx7_1 -lmosekjava7_1 -lmosekscopt7_1 -liomp5" --with-cplex-incdir=/usr/local/cplex/include/ilcplex --with-cplex-lib="-L/usr/local/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread"
make -j 10 install
cd ..
# configure and install Clp
mkdir -p Clp
cd Clp
../../Clp/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Cgl
mkdir -p Cgl
cd Cgl
../../Cgl/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Alps
mkdir -p Alps
cd Alps
../../Alps/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Bcps
mkdir -p Bcps
cd Bcps
../../Bcps/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi-Conic
mkdir -p OsiConic
cd OsiConic
../../OsiConic/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi-Mosek
mkdir -p OsiMosek
cd OsiMosek
../../OsiMosek/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi-Cplex
mkdir -p OsiCplex
cd OsiCplex
../../OsiCplex/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install ipopt
mkdir -p Ipopt
cd Ipopt
../../Ipopt/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install OsiIpopt
mkdir -p OsiIpopt
cd OsiIpopt
../../OsiIpopt/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Cola
mkdir -p Cola
cd Cola
../../Cola/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install CglConic
mkdir -p CglConic
cd CglConic
../../CglConic/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and compile disco
mkdir -p DisCO
cd DisCO
../../configure --prefix=$build_dir --with-soco-solver=oa
make -j 10 install
