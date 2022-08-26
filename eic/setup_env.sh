#!/bin/bash
source /opt/detector/setup.sh

EIC_DIR = $PWD
#EIC_DIR=/p/lustre2/ftorales/generate_data/eic
#EIC_DIR=/p/lustre2/marratia/generate_data/eic
#EIC_DIR=/p/lustre1/dongwi1/analysis/hip/generate_data/eic
#printf "EIC_DIR = ${EIC_DIR} \n"

export LD_LIBRARY_PATH=${EIC_DIR}/development/lib:$LD_LIBRARY_PATH
export PATH=${EIC_DIR}/development/bin:$PATH
export DETECTOR_PATH=${EIC_DIR}/athena
export JUGGLER_DETECTOR=athena
export JUGGLER_INSTALL_PREFIX=${EIC_DIR}/development

#HDF5                                                                                                                                                                    
export HDF5_DIR=${EIC_DIR}/../to_hdf5/hdf5-1.8.22/hdf5                                                                        
export PATH=$PATH:${HDF5_DIR}/bin                                                                                             
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib

#TODO: Reinstate ability to run without specifying $EIC_DIR
