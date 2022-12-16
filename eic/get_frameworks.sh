#!/bin/bash

#EIC container paths
export LD_LIBRARY_PATH=$EIC_SHELL_PREFIX/lib:$LD_LIBRARY_PATH
export PATH=$EIC_SHELL_PREFIX/bin:$PATH 

#Beamline
git clone https://github.com/eic/ip6.git
cd ip6
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EIC_SHELL_PREFIX
make install -j8
cd ../..

# Hadron endcap geometry
git clone https://github.com/eiccodesign/eic_geometry.git hadron_endcap
cd hadron_endcap
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EIC_SHELL_PREFIX
make install -j8
cd ../..

source setup_env.sh