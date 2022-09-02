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
cp -r ip6/ip6 athena/

git clone https://eicweb.phy.anl.gov/EIC/juggler.git
cd juggler/
git reset --hard eda9e3ed1265f1a3da2f8a197eb7bd26cb7b9f77
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../development -DCMAKE_CXX_STANDARD=20
make install
cd ../..
pwd

source setup_env.sh

git clone https://eicweb.phy.anl.gov/EIC/benchmarks/reconstruction_benchmarks.git
cd reconstruction_benchmarks/
git checkout ai_codesign #includes changes from Miguel and Fernando
git branch # just to check this is on the ai_codesign branch
./bin/get_calibrations
