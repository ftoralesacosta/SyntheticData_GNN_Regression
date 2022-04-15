First. Get the container: https://eicweb.phy.anl.gov/containers/eic_container/
and follow instructions

Then: 

> mkdir development


> export LD_LIBRARY_PATH=$PWD/development/lib:$LD_LIBRARY_PATH


> export PATH=$PWD/development/bin:$PATH

Then get and install the detector (ATHENA):

> git clone https://eicweb.phy.anl.gov/EIC/detectors/athena.git

> cd athena


> mkdir build && cd build


> cmake .. -DCMAKE_INSTALL_PREFIX=../../development


> make install


> cd ../..

Thne install beamline: 

> git clone https://eicweb.phy.anl.gov/EIC/detectors/ip6.git


> cd ip6


> mkdir build && cd build


> cmake .. -DCMAKE_INSTALL_PREFIX=../../development


> make install


> cd ../..

Then copy beamline compact files to detector
> cp -r ip6/ip6 athena/

Then install juggler (deals with reconstruction)

> git clone https://eicweb.phy.anl.gov/EIC/juggler.git


> cd juggler/


> mkdir build && cd build


> cmake .. -DCMAKE_INSTALL_PREFIX=../../development -DCMAKE_CXX_STANDARD=20


> make install


> cd ../..

Set environmental variables needed to generated data (run reconstruction_benchmarks). You might want to add this to a your bashrc
 
#!/bin/bash                                                                     
source /opt/detector/setup.sh

export LD_LIBRARY_PATH=$PWD/development/lib:$LD_LIBRARY_PATH

export PATH=$PWD/development/bin:$PATH

export DETECTOR_PATH=$PWD/athena

export JUGGLER_DETECTOR=athena

export JUGGLER_INSTALL_PREFIX=$PWD/development



#get reconstruction code and run (for example the full calo cluster benchmark)

git clone https://eicweb.phy.anl.gov/EIC/benchmarks/reconstruction_benchmarks/-/tree/ai_codesign

Note that the one above has a few changes we made that are not in the main branch of reconstruction_bechmarks. 


> cd reconstruction_benchmarks/


> ./bin/get_calibrations                                      #You only need this once


> bash benchmarks/clustering/full_cal_clusters.sh -p "electron" -n 100 --pmin 0.99 --pmax 1.01 -t test1

The last command will use single-particle gun to generate 100 electrons with a momentum in range 0.99-1.01 GeV. 
Then, digitization and reconstruction will be run. 
Output will contain sim_TAGNAME.root and rec_TAGNAME.root files, which correspond to the G4-level and reconstructed level respectively (digitization, clustering). 
The files will contain a ROOT TTree with the generated particles, hits at the G4 level, hits after digitization, and clusters. 
We will mainly be working with hits either before or after digitization.
