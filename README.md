** This Branch aims to automate the process of setting up the environment for ATHENA simulation and reconstruction. There are scripts to download and mount the EIC singularity container, as well as installing __specific__ commits of `Athena`, `reconstruction_benchmarks`, `ip6`, and `juggler` repositories. Lastly, it introduces HDF5 for the final data format.
*** Prerequisites: linux (requires singularity v3+) and MacOS (requires docker) 

1. Download the EIC container
> ./get_eic-container.sh

This will download the container, enter the container, and create a few required directories for the next step.

2. Grab specific commits of the repositories mentioned above
> source get_frameworks.sh

This downloads the specific commits used to generate data. It builds them and and then sources the setup_env.sh script to set the proper paths
Make sure you're still in the container when running this.

3. Try the simulation
> bash benchmarks/clustering/full_cal_clusters.sh -p "pion+" -n 100 --pmin 19.99 --pmax 20.01 -t pionplus_20Gev_test

Also make sure to still be in the container when running this.
This last command will use single-particle gun to generate 100 electrons with a momentum in range 0.99-1.01 GeV. 
Then, digitization and reconstruction will be run. 
Output will contain sim_TAGNAME.root and rec_TAGNAME.root files, which correspond to the G4-level and reconstructed level respectively (digitization, clustering). 
The files will contain a ROOT TTree with the generated particles, hits at the G4 level, hits after digitization, and clusters. 
We will mainly be working with hits either before or after digitization.

List of particles for the particle gun are found in /eic/athena/reconstruction_benchmarks/benchmarks/imaging_ecal/scripts/gen_particles.py

    "pion0": (111, 0.1349766),       # pi0                                                                  
    "pion+": (211, 0.13957018),      # pi+                                                                  
    "pion-": (-211, 0.13957018),     # pi-                                                                  
    "kaon0": (311, 0.497648),        # K0                                                                   
    "kaon+": (321, 0.493677),        # K+                                                                   
    "kaon-": (-321, 0.493677),       # K-                                                                   
    "proton": (2212, 0.938272),      # proton                                                               
    "neutron": (2112, 0.939565),     # neutron                                                              
    "electron": (11, 0.51099895e-3), # electron                                                             
    "positron": (-11, 0.51099895e-3),# positron                                                             
    "photon": (22, 0),               # photon                                                               
    "muon-": (13, 0.1056583),         # muon-                                                               
    "muon+": (-13, 0.1056583),         # muon+  

4. Download and install HDF5
While still inside the container, navigate to `to_hdf5`, and `source` the grab hdf5 script.
> make root_to_hdf5

> ./root_to_hdf5 [input root file] [new_hdf5_file.hdf5]
