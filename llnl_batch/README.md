# LLNL Slurm Batch
This directory contains files based on Aarons example code on running parrellel SLURM batch jobs. 
It combines progress made by Bishoy with Aarons example to submit a job array of batch jobs 
with a nice argument parser written in python.

There are some caveats: Because the jobs must be run in a singularity container, 
it is difficult to pass arguments to the script that actually runs the simulation. 
For this reason, as much of the EIC stuff as possible is encapsulated in the eic/run_gun.sh script, but please note that the __particle species__ and __number of events__ are hardcoded in this script.

## Steps before Submitting Jobs:

1. The `EIC_DIR` variable should should be the full path including the `eic` directory.
The variable should be edited in these files:
- generate_data/eic/__eic-shell__
- generate_data/llnl_batch/__run-subjob.sh__
- generate_data/eic/__run_gun.sh__
- generate_data/eic/reconstruction_benchmarks/benchmarks/clustering/__full_cal_clusters.sh__

An example for user ftorales on LLNL login is: 

/p/lustre2/ftorales/generate_data/eic

2. Edit `eic-shell` by adding this line right under the $EIC_DIR declaration:
- export SIF=$EIC_DIR/working_image.sif

This edit sets the full path to the image, such that the batch scripts can be run in pretty much any directory, run well, and save the output to the correct directory.

## To Submit the jobs
This command __s__ubmits __n__=10 subjobs to the batch system launches dependant 
__m__erge job to add all the root files together:
> ./slurm-example.py -j [job name] -s -m -n 10 

### TODO:
- improve the move command in run-subjob.sh to be more robust, and include the gen and hepmc files
- find a way to pass arguments to the generation command (full_cal_clusters.sh) from outside the conainer.
`singularity exec [container] [command] [arguments?]` might be the solution.
