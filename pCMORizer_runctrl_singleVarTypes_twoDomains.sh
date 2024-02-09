#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --threads-per-core=1
#SBATCH --time=00:15:00
#SBATCH --partition=dc-cpu-devel
#SBATCH --mail-type=all
#SBATCH --mail-user=k.goergen@fz-juelich.de
#SBATCH --account=jjsc39
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# USAGE="export Y=2013 && sbatch --export=ALL,Y=$Y --job-name=pCMORizer$Y pCMORizer_runctrl_singleVarTypes_twoDomains.sh"
# USAGE="export Y=2004 && sbatch --export=ALL,Y=$Y --job-name=pCMORizer$Y pCMORizer_runctrl_singleVarTypes_twoDomains.sh"

source loadenv.JURECA-DC_2023_Intel-PSMPI.ini

let Y=$Y
let nvar=16

dir_work=$(pwd)
mkdir -p ${dir_work}/${Y}
cd "${dir_work}/${Y}"

cp -f ../runctrl.current.nml_template_d01_DA runctrl.current.nml_d01
sed -i "s/__YYYY__/$Y/g" runctrl.current.nml_d01
sed -i "s/__nvar__/$nvar/g" runctrl.current.nml_d01
ln -sf runctrl.current.nml_d01 runctrl.current.nml

cp -f ../runctrl.vars.std_sfc.nml .

cp -f ../pCMORizer . && chmod a+x pCMORizer

#ln -sf "/p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/simres/d01/*/wrfout_d0?_${Y}*" .
#ln -sf /p/largedata/jjsc39/jjsc3900/sim/CORDEX-FPSCONV_EUR-15-ALP-3_SMHI-EC-EARTH_historical_r12_FZJ-IDL-WRF381DA_v00aJurecaDcCpuProdPrjTt19952005/simres/d01/*/wrfout_d01_${Y}* .
ln -sf /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/tmp/test_input_UNICAN/wrfout_d01_${Y}* .
#ln -sf /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/tmp/test_input_FZJ/wrfout_d01_${Y}* .

#srun --exact --cpu-bind=threads --distribution=block:cyclic:fcyclic --ntasks=$nvar ./pCMORizer > ./log.txt 2>&1
srun --ntasks=$nvar ./pCMORizer > ./log.txt 2>&1

exit 0
