#!/bin/bash

EIC_DIR=/p/lustre2/ftorales/generate_data/eic
# EIC_DIR=$PWD

source /opt/detector/setup.sh                                                                          

PARTICLE=$1
export LD_LIBRARY_PATH=$EIC_DIR/development/lib:$LD_LIBRARY_PATH                                           
export PATH=$EIC_DIR/development/bin:$PATH                                                                 
export DETECTOR_PATH=$EIC_DIR/athena                                                                       
export JUGGLER_DETECTOR=athena                                                                         
export JUGGLER_INSTALL_PREFIX=$EIC_DIR/development                                                         
                                                                                                       
#HDF5                                                                                                  
export HDF5_DIR=${EIC_DIR}/../to_hdf5/hdf5-1.8.22/hdf5                                                     
export PATH=$PATH:${HDF5_DIR}/bin                                                                      
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib                                                

cd $EIC_DIR/reconstruction_benchmarks/
pwd
bash benchmarks/clustering/full_cal_clusters.sh -p pion+ -n 1 --pmin 0.0 --pmax 100. -t piplus_out


