#!/bin/bash
source /opt/detector/setup.sh

EIC_DIR=${EIC_DIR:-${PWD}} #If not defined, set it to pwd
printf "EIC_DIR from setup_env.sh = ${EIC_DIR} \n"

export LD_LIBRARY_PATH=${EIC_DIR}/development/lib:$LD_LIBRARY_PATH
export PATH=${EIC_DIR}/development/bin:$PATH
export DETECTOR_PATH=${EIC_DIR}/athena
export JUGGLER_DETECTOR=athena
export JUGGLER_INSTALL_PREFIX=${EIC_DIR}/development

#HDF5                                                                                                                                                                    
export HDF5_DIR=${EIC_DIR}/../to_hdf5/hdf5-1.8.22/hdf5                                                                        
export PATH=$PATH:${HDF5_DIR}/bin                                                                                             
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib
