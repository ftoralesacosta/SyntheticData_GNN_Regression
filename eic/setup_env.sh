#!/bin/bash
source /opt/detector/setup.sh
export LD_LIBRARY_PATH=$PWD/development/lib:$LD_LIBRARY_PATH
export PATH=$PWD/development/bin:$PATH
export DETECTOR_PATH=$PWD/athena
export JUGGLER_DETECTOR=athena
export JUGGLER_INSTALL_PREFIX=$PWD/development

#HDF5
# export HDF5_DIR=${PWD}/athena/hdf5-1.8.22/hdf5
# export PATH=$PATH:${HDF5_DIR}/bin
# echo $HDF5
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib

#HDF5                                                                                                                                                                    
# export HDF5_DIR=${PWD}/../to_hdf5/hdf5-1.8.22/hdf5                                                                                                                           
export HDF5_DIR=/clusterfs/ml4hep_nvme2/ftoralesacosta/generate_data/to_hdf5/hdf5-1.8.22/hdf5
#FIXME: relative path not working, temp hard coded path
#Probably need to move directories around for ease
export PATH=$PATH:${HDF5_DIR}/bin                                                                                                                                        
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib
echo $HDF5_DIR                                                                                                                                                               
