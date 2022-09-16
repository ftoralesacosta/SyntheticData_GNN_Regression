#!/bin/bash
gzip -cd hdf5-1.13.2.tar.gz | tar xvf -
cd hdf5-1.13.2

./configure --enable-cxx
ldconfig
#./config -flags for options. e.g. --enable-parallel
make
make install

export HDF5_DIR=${PWD}/hdf5
export PATH=$PATH:${HDF5_DIR}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib

echo $HDF5_DIR
echo "HDF5 Install Done"
