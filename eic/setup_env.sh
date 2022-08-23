#!/bin/bash
source /opt/detector/setup.sh

printf "\n CURRENT WORKING DIRECTORY = {$PWD}"
EIC_DIR=/p/lustre2/ftorales/generate_data/eic
printf "\n EIC_DIR set to ${EIC_DIR}" 

export LD_LIBRARY_PATH=${EIC_DIR}/development/lib:$LD_LIBRARY_PATH
export PATH=${EIC_DIR}/development/bin:$PATH
export DETECTOR_PATH=${EIC_DIR}/athena
export JUGGLER_DETECTOR=athena
export JUGGLER_INSTALL_PREFIX=${EIC_DIR}/development

#HDF5                                                                                                                                                                    
export HDF5_DIR=${EIC_DIR}/../to_hdf5/hdf5-1.8.22/hdf5                                                                        
export PATH=$PATH:${HDF5_DIR}/bin                                                                                             
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib

#export HDF5_DIR=/clusterfs/ml4hep_nvme2/ftoralesacosta/generate_data/to_hdf5/hdf5-1.8.22/hdf5 
#For Reference
