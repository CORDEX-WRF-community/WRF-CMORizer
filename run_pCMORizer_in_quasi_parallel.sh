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
source loadenv.UCAN-IFCA_intel.ini

# Adjust to your situations
export YEAR=$1 	
export YEAR_next=$(date -u --date="$YEAR-01-01 1 year" '+%Y')
export DOM=$2 
export VARNAME=$3
export PROJECT=$4
export FREQ=$5

if [[ ${PROJECT} == "I4C" ]]; then
  export dir_data_in="set-full-path-to-raw-wrfoutput-files"
elif [[ ${PROJECT} == "EUROCORDEX" ]]; then
  export dir_data_in="set-full-path-to-raw-wrfoutput-files"
else
  echo "Provide full path to you raw wrfout files."
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

# Create working directory
dir_home=$(pwd)
dir_work=${dir_home}/${PROJECT}/${DOM}/${FREQ}/${YEAR}
mkdir -p ${dir_work}; cd ${dir_work}
ln -sf ${dir_data_in}/wrf*_${DOM}_${YEAR}* ${dir_work}/
ln -sf ${dir_data_in}/wrf*_${DOM}_${YEAR_next}-01-01* ${dir_work}/

# Create working directory per variable
mkdir -p ${dir_work}/${VARNAME}
cd ${dir_work}/${VARNAME}

# Adapt general namelist for the seleted variabels
cp -f ${dir_home}/runctrl.current.nml_template_${DOM}_${PROJECT} ${dir_work}/${VARNAME}/runctrl.current.nml_${DOM}
sed -i "s/__YYYY__/$Y/g" runctrl.current.nml_${DOM}
sed -i "s/__nvar__/$nvar/g" runctrl.current.nml_${DOM}
ln -sf runctrl.current.nml_${DOM} runctrl.current.nml

# Generate namelist for the selected variables
cp -f ${dir_home}/generate_vars_namelist.py ${dir_work}/${VARNAME}/
cp -f ${dir_home}/CORDEX_CMIP6_variables.csv ${dir_work}/${VARNAME}/
python generate_vars_namelist.py ${VARNAME}
mv runctrl.vars.${VARNAME}.nml runctrl.vars.nml

# Copy file with the pCMORizer version
cp ${dir_home}/VERSION.txt ${dir_work}/${VARNAME}/

# Compile the code for the specific varlist
cp -f ${dir_home}/Makefile ${dir_work}/${VARNAME}/
cp -f ${dir_home}/pCMORizer.f90 ${dir_work}/${VARNAME}/pCMORizer.f90
sed -i "s/DO ifrq = 1, 1, 1/DO ifrq = $freq_id, $freq_id, 1/g" ${dir_work}/${VARNAME}/pCMORizer.f90
make veryclean
make
sleep 2
mv pCMORizer pCMORizer.exe
sleep 2

# running the code
cd ${dir_work}/${VARNAME}/
#./pCMORizer.exe                      # run directly in the interface
srun --cpu-bind=cores ./pCMORizer.exe # when sending job and running on nodes
exit 0
