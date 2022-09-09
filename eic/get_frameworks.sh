#!/bin/bash

#EIC container paths
export LD_LIBRARY_PATH=$PWD/development/lib:$LD_LIBRARY_PATH
export PATH=$PWD/development/bin:$PATH 

# Hadron endcap geometry
git clone https://github.com/eiccodesign/eic_geometry.git hadron_endcap
cd hadron_endcap
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../development
make install
cd ../..

#Beamline
git clone https://eicweb.phy.anl.gov/EIC/detectors/ip6.git
cd ip6
git reset --hard 45ed8e8a77862dd727aabe5aad943e4adccec020
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../development
make install
cd ../..

git clone https://eicweb.phy.anl.gov/EIC/juggler.git
cd juggler/
git reset --hard fa539623403c86259f9538bfce6a0defef517d95
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../development -DCMAKE_CXX_STANDARD=20
make install
cd ../..
pwd

source setup_env.sh
