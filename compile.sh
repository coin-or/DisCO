#!/bin/bash
mkdir build
cd build
build_dir=$PWD
inc_dir=${build_dir%%/}/include
lib_dir=${build_dir%%/}/lib
pkg_dir=${lib_dir%%/}/pkgconfig
PKG_CONFIG_PATH=${pkg_dir}:$PKG_CONFIG_PATH
export CXXFLAGS="-std=c++11 -g"
# configure and install Blas
mkdir -p ThirdParty/Blas
cd ThirdParty/Blas
../../../ThirdParty/Blas/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Lapack
mkdir Lapack
cd Lapack
../../../ThirdParty/Lapack/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install ASL
mkdir ASL
cd ASL
../../../ThirdParty/ASL/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install HSL
mkdir HSL
cd HSL
../../../ThirdParty/HSL/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Metis
mkdir Metis
cd Metis
../../../ThirdParty/Metis/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Mumps
mkdir Mumps
cd Mumps
../../../ThirdParty/Mumps/configure --prefix=$build_dir
make -j 10 install
cd ../..
# configure and install CoinUtils
mkdir CoinUtils
cd CoinUtils
../../CoinUtils/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi
mkdir Osi
cd Osi
../../Osi/configure --prefix=$build_dir --with-mosek-incdir=/usr/local/mosek/7.1/tools/platform/linux64x86/h --with-mosek-lib="-L/usr/local/mosek/7.1/tools/platform/linux64x86/bin -lmosek64 -lmosekxx7_1 -lmosekjava7_1 -lmosekscopt7_1 -liomp5" --with-cplex-incdir=/usr/local/cplex/include/ilcplex --with-cplex-lib="-L/usr/local/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread"
make -j 10 install
cd ..
# configure and install Clp
mkdir Clp
cd Clp
../../Clp/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Cgl
mkdir Cgl
cd Cgl
../../Cgl/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Alps
mkdir Alps
cd Alps
../../Alps/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Bcps
mkdir Bcps
cd Bcps
../../Bcps/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi-Conic
mkdir OsiConic
cd OsiConic
../../OsiConic/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi-Mosek
mkdir OsiMosek
cd OsiMosek
../../OsiMosek/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi-Cplex
mkdir OsiCplex
cd OsiCplex
../../OsiCplex/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install ipopt
mkdir Ipopt
cd Ipopt
../../Ipopt/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install OsiIpopt
mkdir OsiIpopt
cd OsiIpopt
../../OsiIpopt/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Cola
mkdir Cola
cd Cola
../../Cola/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install CglConic
mkdir CglConic
cd CglConic
../../CglConic/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and compile disco
mkdir DisCO
cd DisCO
../../configure --prefix=$build_dir
make -j 10 install
