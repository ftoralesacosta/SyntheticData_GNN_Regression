#!/bin/bash
# wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.1/src/hdf5-1.13.1.tar.gz
# gzip -cd hdf5-1.13.1.tar.gz | tar xvf -
# cd hdf5-1.13.1
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.22/src/hdf5-1.8.22.tar.gz
gzip -cd hdf5-1.8.22.tar.gz | tar xvf -
cd hdf5-1.8.22
./configure --enable-cxx
ldconfig
#./config -flags for options. e.g. --enable-parallel
make
make install

export HDF5_DIR=${PWD}/hdf5
export PATH=$PATH:${HDF5_DIR}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib

cd ../

echo $HDF5_DIR
echo "HDF5 Install Done"
