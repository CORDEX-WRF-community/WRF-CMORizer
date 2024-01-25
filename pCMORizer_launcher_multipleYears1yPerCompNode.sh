#!/bin/bash

# start each year on a seperate compute node
#for i in {2006,2011,2012,2013,2014} ; do
for i in {2000..2014} ; do

  echo $i

  # launch d01 and d02 at the same time on same compute node, but only one set of vars -- previous default
  #export Y=$i && sbatch --export=ALL,Y=$Y --job-name=pCMORizer$Y pCMORizer_runctrl_singleVarTypes_twoDomains.sh

  # launch d01 and d02 at the same time on different compute nodes and all vars for one year at the same time on one compute node -- 2023 default
  export Y=$i && export DOM=d01 && sbatch --export=ALL,Y=$Y,DOM=$DOM --job-name=pCMORizer$Y$DOM pCMORizer_runctrl_allVarTypes_singleDomain.sh
  #export Y=$i && export DOM=d02 && sbatch --export=ALL,Y=$Y,DOM=$DOM --job-name=pCMORizer$Y$DOM pCMORizer_runctrl_allVarTypes_singleDomain.sh

done
