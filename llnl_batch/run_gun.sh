#!/bin/bash

#===========================
# User Specified Directories
#===========================
export EIC_DIR=/p/lustre1/milton3/eiccodesign_10_23/generate_data/eic
export OUTPUT_DIR=/p/lustre1/milton3/eiccodesign_10_23/generate_data/eic/pi+_data
export DETECTOR_PATH=${EIC_DIR}/local/share/hadron_endcap
export num_events="50"
export theta_min=17.0
export theta_max=17.0
export PARTICLE="pi0"
export ENERGY_MIN="1"
export ENERGY_MAX="150"
export DISTRIBUTION="log10continuous" # Energy distribution options: fixed, uniform, Gaussian, discrete, log10continuous

JOB_ARRAY_ID=${SLURM_ARRAY_TASK_ID}                           
SLURM_JOB_ID=${SLURM_JOB_ID}

function print_the_help {
  echo "USAGE: ${0} -n <nevents> -part <"\"particle\""> -p <momentum> -thmin <theta_min> -thmax <theta_max> -dist <"\"dist\""> "
  echo "  OPTIONS: "
  echo "    -n, --nevents           Number of events to generate"
  echo "    -d, --directory         Directory for storing output root files, as well as temporary job scripts"
  echo "    -j, --jobname           Name of Job used for labelling large batches of submissions"
  echo "    -pmin, --energy_min     Minimum particle momentum (GeV)"
  echo "    -pmax, --energy_max     Maximum particle momentum (GeV)"
  echo "    -thmin, --theta_min     Theta Minimum (deg)"
  echo "    -thmax, --theta_max     Theta Maximum (deg)"
  echo "    -part,--particle        Particle Species"
  echo "                            Allowed Particles: pi0, pi+, pi-, ka0, ka+, ka-, proton, neutron, e-, e+,photon"
  echo "    -dist,--distribution    Energy distribution"
  echo "                            Allowed options: fixed, uniform, gaussian, log10continuous, log10discrete"
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
     -d|--directory)
      export OUTPUT_DIR="$2"
      shift #past argument
      shift #past value
      ;;
    -j|--jobname)
      export JOB_NAME="$2"
      shift #past argument
      shift #past value
      ;;
    -jid)
      export JOB_ID="$2"
      shift #past argument
      shift #past value                                                                                                                         
      ;;
    -part|--particle)
      export PARTICLE="$2"
      shift #past argument
      shift #past value
      ;;
    -n|--nevents)
      export NUM_EVENTS="$2"
      shift #past argument
      shift #past value
      ;;
    -thmin)
      export THETA_MIN="$2"
      shift #past argument
      shift #past value
      ;;
    -thmax)
      export THETA_MAX="$2"
      shift #past argument
      shift #past value
      ;;
    -pmin)
      export ENERGY_MIN="$2"
      shift #past argument
      shift #past value
      ;;
    -pmax)
      export ENERGY_MAX="$2"
      shift #past argument
      shift #past value
      ;;
    -dist)
      export ENERGY_DISTRIBUTION="$2"
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

NAME_TAG="${SLURM_JOB_ID}_${JOB_ARRAY_ID}"
echo "________ ${NAME_TAG}"
#================================================
# Create file to execute upon entering eic-shell
#================================================
#Make the script that runs inside the container
export GENSCRIPTNAME="${NAME_TAG}.sh"
if [ -f "${GENSCRIPTNAME}" ]; then
    echo "${GENSCRIPTNAME} exists and will be removed..."
    rm "${GENSCRIPTNAME}"
fi
cd ${EIC_DIR}
#Write to Script
echo "#!/bin/bash" > ${GENSCRIPTNAME}
echo -en "\n" >> ${GENSCRIPTNAME}
echo "source ${EIC_DIR}/setup_env.sh" >> ${GENSCRIPTNAME}
echo "bash  ${DETECTOR_PATH}/scripts/run_sim_hepmc.sh -part \"${PARTICLE}\" -n ${num_events} -thmin ${theta_min} -thmax ${theta_max} -pmin ${ENERGY_MIN} -pmax ${ENERGY_MAX} -dist \"${DISTRIBUTION}\" -jid ${NAME_TAG}"  >> ${GENSCRIPTNAME}
chmod 700 ${GENSCRIPTNAME}
bash ${EIC_DIR}/eic-shell -- ./${GENSCRIPTNAME}
echo "bash  ${DETECTOR_PATH}/scripts/run_sim_hepmc.sh -part \"${PARTICLE}\" -n ${num_events} -thmin ${theta_min} -thmax ${theta_max} -pmin ${ENERGY_MIN} -pmax ${ENERGY_MAX} -dist \"${DISTRIBUTION}\" -jid ${NAME_TAG}"
 
#=========================================
# Move Output Files to Correct Directories
#=========================================
#Make sure that similar file does not exist to avoid complications
printf "\n Moving files to ${OUTPUT_DIR} \n"
info_string="${PARTICLE}_${DISTRIBUTION}_${ENERGY_MIN}GeV-${ENERGY_MAX}GeV_theta_${theta_min}deg-${theta_max}deg_${NAME_TAG}"

mv ${GENSCRIPTNAME} ${TEMPSCRIPT_DIR}
mv ${EIC_DIR}/reco_${info_string}.edm4hep.root ${RECO_DIR}
mv ${EIC_DIR}/sim_${info_string}*.root ${GEN_DIR}
mv ${EIC_DIR}/gen_${info_string}*.hepmc ${HEPMC_DIR}

#Make the stored files read/writeable
chmod -R 777 ${RECO_DIR}
chmod -R 777 ${GEN_DIR}
chmod -R 777 ${HEPMC_DIR}
