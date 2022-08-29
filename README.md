### This Repo aims to automate the process of setting up the environment for ATHENA simulation and reconstruction. There are scripts to download and mount the EIC Singularity container, as well as installing __specific__ commits of `Athena`, `reconstruction_benchmarks`, `ip6`, and `juggler` repositories. Lastly, it introduces HDF5 for the final data format.
#### Prerequisites: linux (requires singularity v3+) and MacOS (requires docker) 

###Note: Due to a breaking update and a lapse in versioning of the EIC image, we must use a backed up Singularity image hosted ourselves. Google Drive is the current solution.

> ./get_eic-container.sh
1. Setup the environment for loading the container

This will download the nightly container (which we no longer use) and creates a few required directories for the next step. It also sets various environment variables for using the container 

2. Download the EIC singularity image from this [GDrive Link](https://drive.google.com/file/d/10WuqchbaVqLZthWtGjth2QMlSfEthw_t/view?usp=sharing)

**Make sure the image is in the `eic` directorty.**
If running on a headless machine, one can try using `gdown` package available through `pip`, e.g.  ```gdown 10WuqchbaVqLZthWtGjth2QMlSfEthw_t```

3. IMPORTANT (Temporary as of 8/29/2022): Remove references to `ecal` in `/generate_data/eic/athena/athena.xml`, Line 123. It should look like this:
```
120   <documentation level="10">                                               
121   ## Central calorimetry                                                   
122   </documentation>                                                         
123   <!-- <include ref="compact/ecal.xml"/> -->                               
124   <include ref="compact/hcal.xml"/>  
```

4. Enter the container downloaded in Step 2  
> ./enter_container.sh

The enter_container script importantly specifies the container to use before running the `./eic-shell` command.
One can set this environment variable themselves with `export SIF=$PWD/working_image.sif`.

5. Grab specific commits of the repositories mentioned above, builds them, and sets a handful of important environment variables
> source get_frameworks.sh

This downloads the specific commits used to generate data. It builds them and and then sources the setup_env.sh script to set a handful of important environment variables inside the container.
Make sure you're still in the container when running this.

__note:__ The reconstruction benchomarks is a separate branch, rather than a specific commit like the other frameworks. It is the [ai_codesign](https://eicweb.phy.anl.gov/EIC/benchmarks/reconstruction_benchmarks/-/tree/ai_codesign) branch. The script automatically checks this branch out.

6. Try the simulation
> bash benchmarks/clustering/full_cal_clusters.sh -p "pion+" -n 100 --pmin 19.99 --pmax 20.01 -t pionplus_20Gev_test

Also make sure to still be in the container when running this.
This last command will use single-particle gun to generate 100 pions with a momentum in range 19.99-20.01 GeV. 
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
While still inside the container, navigate to `to_hdf5`, and `source` the grab hdf5 script. Then run the next two commands to compile the code.
> make root_to_hdf5

> ./root_to_hdf5 [input root file] [new_hdf5_file.hdf5]

### Once setup
To re-enter the container, one just need to run two commands:
> ./enter_container.sh

> source setup_env.sh

### Note on Running Batch Jobs:
Two scripts exist in the eic directory which allow you to push jobs to Livermore Clusters. Please take the time to edit both of these in order to set the appropriate paths for storing the log files from slurm (push2clust.sh: L12). Secondly, you will need to edit the contJob.cmd, line 30, to your preferred variable values. These scripts should even if the contJob.cmd is untouched (use default values).
To push jobs you the cluster you simply execute the push2clust.sh as follows:
> ./push2clust pion+

You can of course change the particle name, but it must match the aforementioned particle name syntax exactly.

TODO: Make this more user-friendly, input is welcome.

### Note and Editing
Most of the frameworks should not be edited, with the exception of the reconstruction `reconstruction_benchmarks` repository. We are using a branch specifically for the ai_codesign project, and changes to the particle gun as well as changes to the structure of the output root file can be done by editing the python files in that repository.

### To Do:
1. Look into https://cloud.sylabs.io/library for hosting the image. Ideally, something that supports the `Singularity pull` command would be best for easy scripting.
2. Script the generation together with the HDF5 conversion for batch jobs
