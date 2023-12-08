#!/bin/bash

# selectively copy set of files of variables from CMORized DRS directory tree
# needed becasue the official VLs require temporal aggregations of specific vars
# as we CMORize variables in a first step at the highest temporal aggregation 
# level, too many variables would be transferred, as a violation of the VL and 
# archival protocol
# run on jsc-cordex, pull data in, rather than push as JSC HPC systems cannot push

# not yeat ready to be used, variable list needs to be expanded based on FPSCONV protocol

# transfer complete experiemnt, one domain at a time
dir_src="/p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_BB/postpro/CMORized/CORDEX-FPSCONV/output/ALP-3/FZJ/ECMWF-ERAINT/evaluation/r1i1p1/FZJ-WRF381BB/fpsconv-x1n2-v1"
dir_target=""

# this is the variable list
declare -a list_aggr_var=( "1hr/tas" "1hr/pr" \
                           "3hr/*" )

# pull into this directory
cd $dir_target && pwd

for ifile in ${!list_aggr_var[@]}
do
  echo "index = ${ifile}, value = ${list_aggr_var[$ifile]}"
  mkdir -p ${list_aggr_var[$ifile]}
  cd ${dir_target}/${list_aggr_var[$ifile]}
  scp -i ~/.ssh/sshkey_rsa_KGo_FPS-CPCM_sftp_JSC kgoergen@jsc-cordex.fz-juelich.de:${dir_src}/${list_aggr_var[$ifile]}/*.nc .
  if [[ $? -ne 0 ]] ; then
    echo "WARNING: scp failed"
  fi
done

#higher level dir
#pull checksum file
#ssh sha256sum

exit 0
