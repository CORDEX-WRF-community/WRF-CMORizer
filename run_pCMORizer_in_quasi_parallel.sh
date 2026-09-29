#!/bin/bash
#SBATCH --job-name=cmor
#SBATCH --output=cmor%j.out
#SBATCH --error=cmor%j.error
#SBATCH --ntasks=1
#SBATCH --qos=meteo_high
##SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=720:00:00
##SBATCH --exclusive
#SBATCH --mem-per-cpu=8G
##SBATCH --mem=60G
#SBATCH --hint=nomultithread
##SBATCH --nodelist=wncompute051
#SBATCH --mail-user=milovacj@unican.es
#SBATCH --partition=wncompute_meteo
##SBATCH --exclude=wncompute051


# Set the enviroment
source loadenv.ini

# Adjust to your situations
export YEAR=$1 	
export YEAR_next=$(date -u --date="$YEAR-01-01 1 year" '+%Y')
export DOM=$2 
export VARNAME=$3
export PROJECT=$4
export FREQ=$5

if [[ ${PROJECT} == "project_name" ]]; then # Enter the project name used as the main namelist extension.
  export dir_data_in="set-full-path-to-raw-wrfoutput-files"
else
  echo "======================== WARNING  ================================="
  echo "No input data directory specified for the given $PROJECT projec.   "
  echo "Please provide the full path to the raw WRF output files by        "
  echo "definding 'dir_data_in' in the 'run_pCMORizer_in_quasi_parallel.sh'"
  echo "==================================================================="
  exit 1
fi

nvar=1 

# Set cases according to the data in the pCMORizer code
case $FREQ in
    1hr)
        freq_id=1 ;;
    3hr)
        freq_id=2 ;;
    6hr)
        freq_id=3 ;;
    day)
        freq_id=4 ;;
    fx)
        freq_id=7 ;;
    *)
        echo "Unknown frequency." ;;
esac

echo "Changing in pCMORizer.f90 freq_id to $freq_id for the frequency $FREQ"

# Create working and running directories
dir_home=$(pwd)
dir_work=${dir_home}/${PROJECT}/${DOM}/${YEAR}; mkdir -p ${dir_work}; 
dir_run=${dir_work}/${FREQ}/${VARNAME}; mkdir -p ${dir_run};


files=$(find "${dir_data_in}" -maxdepth 1 -type f -name "wrf*_${DOM}_${YEAR}*" -print)
files_next=$(find "${dir_data_in}" -maxdepth 1 -type f -name "wrf*_${DOM}_${YEAR}*" -print)

if [[ -n "${files}" ]]; then
  ln -sf ${dir_data_in}/wrf*_${DOM}_${YEAR}* "${dir_work}/"
  if [[ -n "${files_next}" ]]; then
    ln -sf ${dir_data_in}/wrf*_${DOM}_${YEAR_next}-01-01* "${dir_work}/"
  fi
else
  echo "ERROR: No WRF output files found for ${DOM} and ${YEAR} in:"
  echo "       ${dir_data_in}"
  exit 1
fi

# Adapt general namelist for the seleted variabels
template_nml=${dir_home}/runctrl.current.nml_template_${DOM}_${PROJECT}
if [ ! -f "$template_nml" ]; then
    echo "Error: Template namelist file $template_nml not found."
    exit 1
fi
cd ${dir_run}
cp -f ${template_nml} ${dir_run}/runctrl.current.nml_${DOM}
sed -i "s/__YYYY__/$Y/g" runctrl.current.nml_${DOM}
sed -i "s/__nvar__/$nvar/g" runctrl.current.nml_${DOM}
ln -sf runctrl.current.nml_${DOM} runctrl.current.nml

# Generate namelist for the selected variables
cp -f ${dir_home}/generate_vars_namelist.py ${dir_run}/
cp -f ${dir_home}/CORDEX_CMIP6_variables.csv ${dir_run}/
python generate_vars_namelist.py ${VARNAME}
mv runctrl.vars.${VARNAME}.nml runctrl.vars.nml

# Copy file with the pCMORizer version
cp ${dir_home}/VERSION.txt ${dir_run}/

# Compile the code for the specific varlist
cp -f ${dir_home}/Makefile ${dir_run}/
cp -f ${dir_home}/pCMORizer.f90 ${dir_run}/pCMORizer.f90
sed -i "s/DO ifrq = 1, 1, 1/DO ifrq = $freq_id, $freq_id, 1/g" ${dir_run}/pCMORizer.f90
make veryclean
make
sleep 2
mv pCMORizer pCMORizer.exe
sleep 2

#./pCMORizer.exe                      # run directly in the interface
srun --cpu-bind=cores ./pCMORizer.exe # when sending job and running on nodes
exit 0
