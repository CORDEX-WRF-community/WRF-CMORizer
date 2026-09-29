#!/bin/bash

#SBATCH --job-name="pCMORizer_checksums_distributed"
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00
#SBATCH --disable-turbomode
#SBATCH --output=tmp.%x.out.%j
#SBATCH --error=tmp.%x.err.%j
#SBATCH --partition=dc-cpu
#SBATCH --account=jjsc39

# sbatch ./pCMORizer_checksums_distributed.sh /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized/CORDEX-FPSCONV

# preparation: chmod ugo-w $(find CORDEX-FPSCONV/ -type f)

inFiles=$@
echo "${inFiles[@]}"
declare -i MAX_PARALLEL=${SLURM_NTASKS} #128
echo "MAX_PARALLEL: $MAX_PARALLEL"
declare -i tmp_parallel_counter=0
declare -i tmp_file_counter=0

function calcChecksum(){
  cd $1
  pwd
  #ls -l *.nc
  sha256sum *.nc > checksum.sha256
  if [[ $? -ne 0 ]] ; then
    echo "WARNING: checksumming failed: $1"
    return 1
  fi
  chmod ugo-w checksum.sha256
}

for inFile in $(find ${inFiles[@]} -name "*.nc" -exec dirname {} \; | uniq -d | sort)
#for inFile in $(find ${inFiles[@]} -name "*.nc" -exec dirname {} \; | uniq -d | sort | grep "psl")
do
  echo "DEBUG: inFile ${inFile}"
  calcChecksum ${inFile} &
  sleep 1
  (( tmp_parallel_counter++ ))
  (( tmp_file_counter++ ))
  if [ $tmp_parallel_counter -ge $MAX_PARALLEL ]; then
    wait
    tmp_parallel_counter=0
    date
  fi
done

wait

exit 0
