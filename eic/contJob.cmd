#!/bin/bash


#SBATCH --job-name=${particle}g4sim
#SBATCH --nodes=1
#SBATCH --time=180:00:00 
#SBATCH -n 36
#SBATCH -p pbatch
#SBATCH -A mlodd

particle=$1

echo "Checking particle name here ... ${particle}"
export SIF=$PWD/working_image.sif

#================================================
# Create file to execute upon entering eic-shell
#================================================
RUNFILE=run-${particle}Job.sh
#Make sure that similar file does not exist to avoid complications
if [ -f "${RUNFILE}" ]; then
    echo "${RUNFILE} exists and will be removed..."
    rm "${RUNFILE}"
fi

echo "#!/bin/bash" > ${RUNFILE}
echo "source setup_env.sh"  >> ${RUNFILE}
echo "cd reconstruction_benchmarks"  >> ${RUNFILE}
echo "The current directory is as follows: $(pwd)"  >> ${RUNFILE}
echo "bash benchmarks/clustering/full_cal_clusters.sh -p \"${particle}\" -n 100000 --pmin 0.0 --pmax 100. -t ${particle}_1e5Evts_0to100GeV"  >> ${RUNFILE}
chmod 700 ${RUNFILE}
bash eic-shell -- ./${RUNFILE}
