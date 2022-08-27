# LLNL Slurm Batch
This directory contains files based on Aaron's example code on running parrellel SLURM batch jobs. 
It combines progress made by Bishoy with Aaron's example to submit a job array to the SLURM scheduler.
with a nice argument parser written in python.

There are some caveats: Because the jobs must be run in a singularity container, 
it is difficult to pass arguments to the script that actually runs the simulation. 
For this reason, as many of the generatation parameters as possible are encapsulated in the eic/run_gun.sh script.
The run_gun.sh script is called from `slurm-batch.py`. This is a python script that interfaces with the SLURM scheduler to submit an array of jobs. The default is set to 4 hours. Change [this line](https://github.com/eiccodesign/generate_data/blob/c421b91b562e1b69825443c665e9031aa9243f52/llnl_batch/slurm-batch.py#L60) to increase or reduce the time per subjob.


## Steps before Submitting Jobs:

1. The `EIC_DIR` variable should should be the full path including the __eic__ directory, and edited according to the user in each of these files:
- generate_data/eic/`run_gun.sh`
- generate_data/eic/reconstruction_benchmarks/benchmarks/clustering/`full_cal_clusters.sh`

An example for user __ftorales__ on LLNL login is: 

> export EIC_DIR="/p/lustre2/ftorales/generate_data/eic"

This edit sets the full path to the Singularity image and parent directory to the rest of the frameworks, such that the batch scripts can be run in pretty much any directory, run well, and save the output to the correct directory.

2. Specify the `submit_dir` when running jobs by doing __one__ of the following:
- always specify the directory to submit and store data using the `-d` flag when calling `slurm-batch.py`
- change the default directory on [this line](https://github.com/eiccodesign/generate_data/blob/c421b91b562e1b69825443c665e9031aa9243f52/llnl_batch/slurm-batch.py#L24) of the `slurm-batch.py` script
The submission directory will coy the appropriate scripts, and more importantly will also become the location for the output root files. Make sure to see the section on storage available, and to not have you home directory under `/g/` set as the submission directory.

## To Submit the jobs
This command **s**ubmits **n**=10 subjobs to the batch system, each generating 1000 pion+ events.
The default time is 4 hours, so this command should finish comfortably within this timeframe. Finally, it specifies the submission/output directory with the `-d` flag.

> ./slurm-batch.py -j [job name] -s -n 10 -ne 1000 -p pion+ -d /p/lscratchh/ftorales/AI-codesign/

## Testing
This command will run fewer 10 pion+ events, using one process, and will run on the login node by omitting the `-s` flag.
./slurm-batch.py -j [job name] -n 1 -ne 10 -p pion+ 

### Storage
`/p/lustre1` and `/p/lustre2` have 20TB of storage each. `lscratchh` is another option, but is more for temporary storage. Check with `quota -v`. Can check your current jobs with `squeue -u [user]`, and cancel them with `scancel [ID]`. 
###TO DO:
Script the Merging of the files using a separate merge.sh and merge job.
