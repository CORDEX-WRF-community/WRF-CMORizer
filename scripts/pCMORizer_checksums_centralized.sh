#!/bin/bash

#SBATCH --job-name="pCMORizer_checksums_centralized"
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

# sbatch ./$0 /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized  # 1 level above CORDEX-FPSCONV

declare -i MAX_PARALLEL=${SLURM_NTASKS}  #128
declare -i tmp_parallel_counter=0
declare -i tmp_file_counter=0

function calcChecksum(){
  sha256sum $1 >> checksum_ALP-3.sha256
  if [[ $? -ne 0 ]] ; then
    echo "WARNING: checksumming failed: $1"
    return 1
  fi
}

cd $1 && pwd
for inFile in $(find . -type f -wholename "*.nc" | sort | grep "EUR-15") ; do
  echo "DEBUG: inFile ${inFile}"
  calcChecksum ${inFile} &
  (( tmp_parallel_counter++ ))
  (( tmp_file_counter++ ))
  if [ $tmp_parallel_counter -ge $MAX_PARALLEL ]; then
    wait
    tmp_parallel_counter=0
  fi
done

wait

exit 0
