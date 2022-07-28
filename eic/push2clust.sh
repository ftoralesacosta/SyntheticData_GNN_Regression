#!/bin/bash

# Make sure that particle to be simulated has been entered
if [ $# -ne 1 ]; then
    echo "Usage: push2clust.sh <PARTICLE_NAME> with exact particle name"
    exit 33
fi

particle=$1
# Create logfile for book-keeping
# Need to change this for preferred logfile dir
export WORKDIR=/p/lustre1/dongwi1/analysis/hip/generate_data/eic/reconstruction_benchmarks
LOGDIR=${WORKDIR}/log
if [ ! -d "$LOGDIR" ]; then
    mkdir -p $LOGDIR
fi

sbatch -o ${LOGDIR}/slurm-${particle}.out -e ${LOGDIR}/slurm-${particle}.err contJob.cmd ${particle}
