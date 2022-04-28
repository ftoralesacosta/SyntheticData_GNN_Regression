#!/bin/bash
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.1/src/hdf5-1.13.1.tar.gz
gzip -cd hdf5-1.13.1.tar.gz | tar xvf -
cd hdf5-1.13.1
#./config -flags for options. e.g. --enable-parallel
make
make install
export HDF5=${PWD}/hdf5
export PATH=$PATH:${HDF5}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5}/lib

cd ../
echo "HDF5 install done"
