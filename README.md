### This Repo aims to automate the process of setting up the environment for ATHENA simulation and reconstruction. There are scripts to download and mount the EIC Singularity container, as well as installing [our hadron endcap geometry](https://github.com/eiccodesign/eic_geometry) and the [IP6 (beamline) repository](https://github.com/eic/ip6). Lastly, it introduces HDF5 for the final data format.
#### Prerequisites: linux (requires singularity v3+) and MacOS (requires docker) 

1. Setup the environment for loading the container
> source get_eic-container.sh

This will download the nightly container and creates a few required directories for the next step. It also sets various environment variables for using the container 

2. Enter the container downloaded in Step 1
> ./eic-shell

3. Run get_frameworks.sh. This grabs the geometry repositories mentioned above, builds them, and sets a handful of important environment variables
> source get_frameworks.sh

This downloads the specific commits used to generate data. It builds them and and then sources the setup_env.sh script to set a handful of important environment variables inside the container.
Make sure you're still in the container when running this.
### Note: You need to use `source setup_env.sh` everytime you want to generate data. It is done automatically the first time in `get_frameworks.sh`.
4. Try the simulation
> $DETECTOR_PATH/scripts/run_sim_hepmc.sh -part "pi+" -n 10 -p 20

Also make sure to still be in the container when running this.
This last command will use a HepMC generator (`$DETECTOR_PATH/hepmc_generation/gen_particles.cxx`) that fires a single particle gun along the proton beam axis to generate 10 pions with a momentum 20 GeV.
Then, digitization and reconstruction will be run with `$DETECTOR_PATH/scripts/hadron_endcap_reco.py`
Output will contain sim_${info_string}.edm4hep.root and rec_${info_string}.edm4hep.root files, which correspond to the G4-level and reconstructed level respectively (digitized only at this point). ${info_string} can be set within run_sim_hepmc.sh and is set by default to include the particle name, energy, and theta values by default. 
The files will contain a ROOT TTree with the generated particles, hits at the G4 level, and hits after digitization and reconstruction.
We will mainly be working with hits after reconstruction.

Examples of particles for the particle gun are:

    "pi0": (111, 0.1349766),       # pi0                                                                  
    "pi+": (211, 0.13957018),      # pi+                                                                  
    "pi-": (-211, 0.13957018),     # pi-                                                                  
    "K0": (311, 0.497648),         # K0                                                                   
    "K+": (321, 0.493677),         # K+                                                                   
    "K-": (-321, 0.493677),        # K-                                                                   
    "proton": (2212, 0.938272),    # proton                                                               
    "neutron": (2112, 0.939565),   # neutron                                                              
    "e-": (11, 0.51099895e-3),     # electron                                                             
    "e+": (-11, 0.51099895e-3),    # positron                                                             
    "gamma": (22, 0),              # photon                                                               
    "mu-": (13, 0.1056583),        # muon-                                                               
    "mu+": (-13, 0.1056583),       # muon+  

The HepMC generator contains parameters for minimum and maximum values of phi and theta and generates values uniformly between these. The momentum distribution can be changed with the `dist` parameter. By default, momentum is fixed (`dist = 0`), but it can also be distributed uniformly with +- 50% variation (`dist = 1`), as a Gaussian with sigma = 0.1*mean (`dist = 2`), and as log uniform, where the energy is chosen randomly in the set {2, 4, 8, 16, 32, 64} GeV.

### Note on geometry: Please read the README in the [eic_geometry](https://github.com/eiccodesign/eic_geometry) repo if you have questions about the simulation geometry. 

5. Download and install HDF5
While still inside the container, navigate to `to_hdf5`, and `source` the grab hdf5 script. Then run the next two commands to compile the code.
> make root_to_hdf5

> ./root_to_hdf5 [input root file] [new_hdf5_file.hdf5]

### Once setup
To re-enter the container, one just need to run two commands:
> ./eic-shell

> source setup_env.sh

### Note on Running Batch Jobs:
Please read the README in the `llnl_batch` directory. Some minor edits simply specifying the output directory and $EIC_DIR environment variable are needed. After that, one should be able to push job arrays to SLURM.

### To Do:
1. Look into https://cloud.sylabs.io/library for hosting the image. Ideally, something that supports the `Singularity pull` command would be best for easy scripting.
2. Script the generation together with the HDF5 conversion for batch jobs
