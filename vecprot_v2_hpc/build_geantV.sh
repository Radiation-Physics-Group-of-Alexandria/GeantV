#!/bin/bash
INSTALL_PREFIX=`pwd`/install
GEANT_DIR=geant
export CC=gcc-4.9
export CXX=g++-4.9
export FC=gfortran-4.9

cd $GEANT_DIR
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
-DUSE_VECGEOM_NAVIGATOR=ON \
-DUSE_ROOT=ON \
-DUSE_HPC=ON \
-DBUILD_ZMQ=ON \
-DVecGeom_DIR=$INSTALL_PREFIX/lib/CMake/VecGeom \
-DUSE_NUMA=OFF \
&& make CMSApp -j4 
