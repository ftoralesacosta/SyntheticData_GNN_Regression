#!/bin/bash

export DETECTOR=hadron_endcap
export DETECTOR_PATH=${EIC_SHELL_PREFIX}/share/${DETECTOR}

export LD_LIBRARY_PATH=${EIC_SHELL_PREFIX}/lib:$LD_LIBRARY_PATH
export PATH=${EIC_SHELL_PREFIX}/bin:$PATH

#HDF5                                                                                                                                                                    
export HDF5_DIR=${EIC_SHELL_PREFIX}/../../to_hdf5/hdf5-1.13.2/hdf5                                                                        
export PATH=$PATH:${HDF5_DIR}/bin                                                                                             
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib