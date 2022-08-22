#!/bin/bash

export EIC_DIR=/p/lustre2/ftorales/generate_data/eic
# export EIC_DIR=/p/lustre1/dongwi1/analysis/hip/generate_data/eic

export SIF=${EIC_DIR}/working_image.sif


FORMATTED_TASK_ID=`python -c 'import os; print("{:03}".format(int(os.getenv("SLURM_ARRAY_TASK_ID"))))'`

function print_the_help {
  echo "USAGE: ${0} -n <nevents> -t <nametag> -p <particle> "
  echo "  OPTIONS: "
  echo "    -n,--nevents     Number of events"
  echo "    -t,--nametag     name tag" #can use this for unique output dir
  echo "    -j,--filename    name of output file "
  echo "    -p,--particle    particle type"
  echo "    --pmin           minimum particle momentum (GeV)"
  echo "    --pmax           maximum particle momentum (GeV)"
  # echo "    -o,--output      output dir for reco,gen, and hepmc files"
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
if [ -f "${FILENAME}" ]; then
    echo "${FILENAME} exists and will be removed..."
    rm "${FILENAME}"
fi

echo "#!/bin/bash" > ${FILENAME}
echo -en "\n" >> ${FILENAME}
echo "source ${EIC_DIR}/setup_env.sh"  >> ${FILENAME}
echo "cd ${EIC_DIR}/reconstruction_benchmarks"  >> ${FILENAME}
echo -en "\n" >> ${FILENAME}
echo "The current directory is as follows: $(pwd)"  >> ${FILENAME}
echo -en "\n" >> ${FILENAME}
echo "bash benchmarks/clustering/full_cal_clusters.sh -p \"${PARTICLE}\" -n ${NEVENTS} --pmin ${PMIN} --pmax ${PMAX} -t ${G4FILENAME} "  >> ${FILENAME}
chmod 700 ${FILENAME}
bash eic-shell -- ./${FILENAME}

EICDIR
