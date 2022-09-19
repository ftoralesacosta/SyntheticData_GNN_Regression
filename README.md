### This Repo aims to automate the process of setting up the environment for ATHENA simulation and reconstruction. There are scripts to download and mount the EIC Singularity container, as well as installing [our hadron endcap geometry](https://github.com/eiccodesign/eic_geometry) and __specific__ commits of `ip6` and `juggler` repositories. Lastly, it introduces HDF5 for the final data format.
#### Prerequisites: linux (requires singularity v3+) and MacOS (requires docker) 


Prerequisites: linux (requires singularity v3+) and MacOS (requires docker) 

Note: Due to a breaking update and a lapse in versioning of the EIC image, we must use a backed up Singularity image hosted ourselves. Google Drive is the current solution.

1. Setup the environment for loading the container
> ./get_eic-container.sh

This will download the nightly container (which we no longer use) and creates a few required directories for the next step. It also sets various environment variables for using the container 

2. Download the EIC singularity image from this [GDrive Link](https://drive.google.com/file/d/10WuqchbaVqLZthWtGjth2QMlSfEthw_t/view?usp=sharing)

**Make sure the image is in the `eic` directorty.**

If running on a headless machine, one can try using `gdown` package available through `pip`. For example, try  ```gdown 10WuqchbaVqLZthWtGjth2QMlSfEthw_t```


3. Enter the container downloaded in Step 2  
> ./enter_container.sh

The enter_container script importantly specifies the container to use before running the `./eic-shell` command.
One can set this environment variable themselves with `export SIF=$PWD/working_image.sif`.

4. Run get_frameworks.sh. This grabs specific commits of the repositories mentioned above, builds them, and sets a handful of important environment variables
> source get_frameworks.sh

This downloads the specific commits used to generate data. It builds them and and then sources the setup_env.sh script to set a handful of important environment variables inside the container.
Make sure you're still in the container when running this.

5. IMPORTANT (Temporary as of 8/29/2022): Remove references to `ecal` in `/generate_data/eic/athena/athena.xml`, Line 123. It should look like this:
```
120   <documentation level="10">                                               
121   ## Central calorimetry                                                   
122   </documentation>                                                         
123   <!-- <include ref="compact/ecal.xml"/> -->                               
124   <include ref="compact/hcal.xml"/>  
```
6. Try the simulation
> bash benchmarks/clustering/full_cal_clusters.sh -p "pion+" -n 100 --pmin 19.99 --pmax 20.01 -t pionplus_20Gev_test

Also make sure to still be in the container when running this.
This last command will use a HepMC generator (`$DETECTOR_PATH/hepmc_generation/gen_particles.cxx`) that fires a single particle gun along the proton beam axis to generate 10 pions with a momentum 20 GeV.
Then, digitization and reconstruction will be run with `$DETECTOR_PATH/scripts/hadron_endcap_reco.py`
Output will contain sim_${info_string}.edm4hep.root and rec_${info_string}.edm4hep.root files, which correspond to the G4-level and reconstructed level respectively (digitization, clustering). ${info_string} can be set within run_sim_hepmc.sh and is set by default to include the particle name, energy, and theta values by default. 
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

The HepMC generator contains parameters for minimum and maximum values of phi and theta and generates values uniformly between these. The momentum distribution can be changed with the `dist` parameter. By default, momentum is fixed (`dist = 0`), but it can also be distributed uniformly with +- 50% variation (`dist = 1`) and as a Gaussian with sigma = 0.1*mean (`dist = 2`). 

6. Download and install HDF5
While still inside the container, navigate to `to_hdf5`, and `source` the grab hdf5 script. Then run the next two commands to compile the code.
> make root_to_hdf5

> ./root_to_hdf5 [input root file] [new_hdf5_file.hdf5]

### Once setup
To re-enter the container, one just need to run two commands:
> ./enter_container.sh

> source setup_env.sh

### Note on Running Batch Jobs:
Please read the README in the `llnl_batch` directory. Some minor edits simply specifying the output directory and $EIC_DIR environment variable are needed. After that, one should be able to push job arrays to SLURM.

### To Do:
1. Look into https://cloud.sylabs.io/library for hosting the image. Ideally, something that supports the `Singularity pull` command would be best for easy scripting.
2. Script the generation together with the HDF5 conversion for batch jobs
