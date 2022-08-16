#! /bin/bash

date

#cd into the submission directory
SUB_DIR=$1
cd $SUB_DIR

particle=$2
N_Events=$3


#python one-liner to format the task ID into something friendlier
FORMATTED_TASK_ID=`python -c 'import os; print("{:03}".format(int(os.getenv("SLURM_ARRAY_TASK_ID"))))'`

TASK_DIRECTORY=job.$FORMATTED_TASK_ID
if [  -d $TASK_DIRECTORY ]
then
    rm -fr $TASK_DIRECTORY
fi

mkdir -p $TASK_DIRECTORY
chmod g+rwsx $TASK_DIRECTORY
cd $TASK_DIRECTORY

#EIC SSSHHHTUFF
EIC_DIR=/p/lustre2/ftorales/generate_data/eic
RECO_DIR=$EIC_DIR/reconstruction_benchmarks
OUTPUT_ROOT_FILE=rec_piplus_out.root
export SIF="${EIC_DIR}/working_image.sif"

#payload
echo 'Hello World!" subjob: ' $FORMATTED_TASK_ID 
echo $TASK_DIRECTORY
pwd

#eic-shell loads the container. run_gun sets env, and runs the gen command
$EIC_DIR/eic-shell -- $EIC_DIR/run_gun.sh 

# mv $RECO_DIR/$OUTPUT_ROOT_FILE ${SUB_DIR}/${TASK_DIRECTORY}/../output/${FORMATTED_TASK_ID}_${particle}_${N_Events}.root
mv $RECO_DIR/$OUTPUT_ROOT_FILE ${1}/${particle}_${N_Events}.root
#mv here is a race condition between competing subjobs.
#can't think of a good solution: full_cal_cluster.sh has
#the files saved its locar dir, and it's difficult to pass a directory
#path as an agrument while outside the container. 

# mv $EIC_DIR/reconstruction_benchmarks/rec*.root ../output/${FORMATTED_TASK_ID}_${particle}_${N_Events}.root

# singularity exec -B cvmfs:/cvmfs -B /p/lustre2/ftorales:/hip_sphenix cvmfs/sphenix.sdcc.bnl.gov/singularity/rhic_sl7_ext.simg sh /hip_sphenix/Singularity/e_Jet_sPHENIX/slurm_batch/my_jetpyth_exe.sh $1
# singularity exec $SIF run_gun.sh ${particle} ${N_Events}
#singularity exec <path-to-container> <command-to-run>
# mv out.root ../output/$FORMATTED_TASK_ID.root"

date
