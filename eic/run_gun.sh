#!/bin/bash

export SCRATCHDIR=/p/lscratchh/dongwi1
export EIC_DIR=/p/lustre1/dongwi1/analysis/hip/generate_data/eic
export SIF=${EIC_DIR}/working_image.sif

#Making temp storage for container run scripts
export TEMPFILEDIR=${SCRATCHDIR}/aicodedesign/tempscripts
if [ ! -d "${TEMPFILEDIR}" ]; then
    mkdir -p ${TEMPFILEDIR}
fi

# Storage of reconstruction and generated ROOT files
export SIMWORKDIR=/p/lustre1/dongwi1/analysis/hip
export RECODIR=${SIMWORKDIR}/recosim
if [ ! -d "${RECODIR}" ]; then
    mkdir -p ${RECODIR}
fi
export GENROOTDIR=${SIMWORKDIR}/gensim
if [ ! -d "${GENROOTDIR}" ]; then
    mkdir -p ${GENROOTDIR}
fi


FORMATTED_TASK_ID=`python -c 'import os; print("{:03}".format(int(os.getenv("SLURM_ARRAY_TASK_ID"))))'`

function print_the_help {
  echo "USAGE: ${0} -n <nevents> -t <nametag> -p <particle> "
  echo "  OPTIONS: "
  echo "    -n,--nevents     Number of events"
  echo "    -t,--nametag     name tag"
  echo "    -p,--particle    particle type"
  echo "    --pmin           minimum particle momentum (GeV)"
  echo "    --pmax           maximum particle momentum (GeV)"
  echo "                     allowed types: pion0, pion+, pion-, kaon0, kaon+, kaon-, proton, neutron, electron, positron, photon"
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
    -j|--makejob)
      export FILENAME="$2"
      shift #past argument
      shift #past value
      ;;
    -t|--nametage)
      export G4FILENAME="$2"
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

#================================================
# Create file to execute upon entering eic-shell
#================================================
#Make sure that similar file does not exist to avoid complications
export TEMPFILENAME=${TEMPFILEDIR}/${FILENAME}
if [ -f "${TEMPFILENAME}" ]; then
    echo "${TEMPFILENAME} exists and will be removed..."
    rm "${TEMPFILENAME}"
fi

echo "#!/bin/bash" > ${FILENAME}
echo -en "\n" >> ${FILENAME}
echo "source ${EIC_DIR}/setup_env.sh"  >> ${FILENAME}
echo "cd ${EIC_DIR}/reconstruction_benchmarks"  >> ${FILENAME}
echo -en "\n" >> ${FILENAME}
echo "bash benchmarks/clustering/full_cal_clusters.sh -p \"${PARTICLE}\" -n ${NEVENTS} --pmin ${PMIN} --pmax ${PMAX} -t ${FORMATTED_TASK_ID}${G4FILENAME}"  >> ${FILENAME}
chmod 700 ${FILENAME}
bash eic-shell -- ./${FILENAME}

mv ${FILENAME} ${TEMPFILEDIR}
mv ${EIC_DIR}/reconstruction_benchmarks/rec_*.root ${RECODIR}
mv ${EIC_DIR}/reconstruction_benchmarks/sim_*.root ${GENROOTDIR}
