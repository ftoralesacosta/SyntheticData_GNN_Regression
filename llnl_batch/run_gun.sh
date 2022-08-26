#!/bin/bash

#===========================
# User Specified Directories
#===========================
# -- Fernando -- 
export EIC_DIR=/p/lustre2/ftorales/generate_data/eic
export OUTPUT_DIR=/p/lscratchh/ftorales/AI-codesign
# -- Bishoy -- 
# export EIC_DIR=/p/lustre1/dongwi1/analysis/hip/generate_data/eic
# export OUTPUT_DIR=/p/lscratchh/dongwi1 # OR /p/lustre1/dongwi1/analysis/hip 
# -- Miguel --
# export EIC_DIR=/p/lustre2/marratia/generate_data/eic
# export OUTPUT_DIR=/p/lscratchh/marratia/AI-codesign


#==================================
# Parse Arguments and Set Variables
#==================================
#DEFAULTS before argument Parsing:
export SIF=${EIC_DIR}/working_image.sif
export NEVENTS="10"
export PMIN="0.0"
export PMAX="100.0"
export PARTICLE="pion+"

JOB_ARRAY_ID=`python -c 'import os; print("{:03}".format(int(os.getenv("SLURM_ARRAY_TASK_ID"))))'`
SLURM_JOB_ID=`python -c 'import os; print("{:03}".format(int(os.getenv("SLURM_JOB_ID"))))'` 

function print_the_help {
  echo "USAGE: ${0} -n <nevents> -d <output_dir> -p <particle> "
  echo "  OPTIONS: "
  echo "    -n,--nevents     Number of events to generate"
  echo "    -d,--directory   Directory for storing output root files, as well as temporary job scripts"
  echo "    -j,--jobname     Name of Job, used for labelling large batches of submissions"
  echo "    --pmin           Minimum particle momentum (GeV)"
  echo "    --pmax           Maximum particle momentum (GeV)"
  echo "    -p,--particle    Particle Species"
  echo "                     Allowed Particles: pion0, pion+, pion-, kaon0, kaon+, kaon-, proton, neutron, electron, positron, photon"
  exit
}

POSITIONAL=()
while [[ $# -gt 0 ]] 
do
  key="$1"
  case ${key} in
    -h|--help)
      shift # past argument
      print_the_help
      ;;
    -j|--jobname)
      export JOB_NAME="$2"
      shift #past argument
      shift #past value
      ;;
    -d|--directory)
      export OUTPUT_DIR="$2"
      shift #past argument
      shift #past value
      ;;
    -p|--particle)
      export PARTICLE="$2"
      shift #past argument
      shift #past value
      ;;
    -n|--nevents)
      export NEVENTS="$2"
      shift #past argument
      shift #past value
      ;;
    --pmin)
      export PMIN="$2"
      shift #past argument
      shift #past value
      ;;
    --pmax)
      export PMAX="$2"
      shift #past argument
      shift #past value
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $2"
      print_the_help
      shift # past argument
      break
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#Check the Directories
printf "\n EIC_DIR set to ${EIC_DIR} \n" 
printf "\n Output Files will be moved to ${OUTPUT_DIR} \n" 

#=====================================================
# Crete Directories for storing Root Files and Scripts
#=====================================================
#Directory for run scripts
export TEMPSCRIPT_DIR=${OUTPUT_DIR}/tempscripts
if [ ! -d "${TEMPSCRIPT_DIR}" ]; then
    mkdir -p ${TEMPSCRIPT_DIR}
fi

#Reco Root Directory 
export RECO_DIR=${OUTPUT_DIR}/recosim
if [ ! -d "${RECO_DIR}" ]; then
    mkdir -p ${RECO_DIR}
fi

#Gen Root Directory
export GEN_DIR=${OUTPUT_DIR}/gensim
if [ ! -d "${GEN_DIR}" ]; then
    mkdir -p ${GEN_DIR}
fi

#Gen Hepmc Directory
export HEPMC_DIR=${GEN_DIR}/hepmc
if [ ! -d "${HEPMC_DIR}" ]; then
    mkdir -p ${HEPMC_DIR}
fi

#Make Unique Base Name
NAME_TAG="${PARTICLE}_${SLURM_JOB_ID}_${JOB_ARRAY_ID}"


#================================================
# Create file to execute upon entering eic-shell
#================================================
#Make the script that runs inside the container
export GENSCRIPTNAME="${NAME_TAG}.sh"
if [ -f "${GENSCRIPTNAME}" ]; then
    echo "${GENSCRIPTNAME} exists and will be removed..."
    rm "${GENSCRIPTNAME}"
fi

#Write to Script
echo "#!/bin/bash" > ${GENSCRIPTNAME}
echo -en "\n" >> ${GENSCRIPTNAME}
echo "source ${EIC_DIR}/setup_env.sh"  >> ${GENSCRIPTNAME}
echo "cd ${EIC_DIR}/reconstruction_benchmarks"  >> ${GENSCRIPTNAME}
echo -en "\n" >> ${GENSCRIPTNAME}
echo "bash benchmarks/clustering/full_cal_clusters.sh -p \"${PARTICLE}\" -n ${NEVENTS} --pmin ${PMIN} --pmax ${PMAX} -t ${NAME_TAG}"  >> ${GENSCRIPTNAME}
chmod 700 ${GENSCRIPTNAME}
bash ${EIC_DIR}/eic-shell -- ./${GENSCRIPTNAME}


#=========================================
# Move Output Files to Correct Directories
#=========================================
#Make sure that similar file does not exist to avoid complications
printf "\n Moving files to ${OUTPUT_DIR} \n"
mv ${GENSCRIPTNAME} ${TEMPSCRIPT_DIR}
mv ${EIC_DIR}/reconstruction_benchmarks/rec_${NAME_TAG}*.root ${RECO_DIR}
mv ${EIC_DIR}/reconstruction_benchmarks/sim_${NAME_TAG}*.root ${GEN_DIR}
mv ${EIC_DIR}/reconstruction_benchmarks/gen_${NAME_TAG}*.hepmc ${HEPMC_DIR}

#Make the stored files read/writeable
chmod -R 777 ${RECO_DIR}
chmod -R 777 ${GEN_DIR}
chmod -R 777 ${HEPMC_DIR}
