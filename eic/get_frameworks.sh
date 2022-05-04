#!/bin/bash

#EIC container paths
export LD_LIBRARY_PATH=$PWD/development/lib:$LD_LIBRARY_PATH
export PATH=$PWD/development/bin:$PATH 

#ATHENA Framework
git clone https://eicweb.phy.anl.gov/EIC/detectors/athena.git
cd athena
git reset --hard f31bac7b8f8ae0f2088321191a1e116640c48296
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../development
make install
cd ../..

#Beamline
git clone https://eicweb.phy.anl.gov/EIC/detectors/ip6.git
git reset --hard f31bac7b8f8ae0f2088321191a1e116640c48296
cd ip6
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../development
make install
cd ../..
cp -r ip6/ip6 athena/

git clone https://eicweb.phy.anl.gov/EIC/juggler.git
git reset --hard ece2557345c0b600b0c564f61fa5f4ae6420b285
cd juggler/
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
