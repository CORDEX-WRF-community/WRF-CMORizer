#!/bin/bash

# nohup ./$0 >> nohup.out &

# selectively copy a set of CMORized files of variables from a DRS directory tree
# from a source to a target destination
# needed because the official VLs require temporal aggregations of specific vars
# but as we CMORize variables in a first step at the highest temporal aggregation 
# level, too many variables would be transferred, as a violation of the VL and 
# archival protocol
# run on jsc-cordex, pull data in, rather than push as JSC HPC systems cannot push
# alternatively: run on judac, push data from judac to jsc-codex

date
echo $USER
echo $0
hostname
uname -a

# transfer complete experiemnt, one domain at a time
dir_src="/p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_BB/postpro/CMORized/CORDEX-FPSCONV/output/ALP-3/FZJ/ECMWF-ERAINT/evaluation/r1i1p1/FZJ-WRF381BB/fpsconv-x1n2-v1"
dir_target="/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/ALP-3/FZJ/ECMWF-ERAINT/evaluation/r1i1p1/FZJ-WRF381BB/fpsconv-x1n2-v1"
dir_fn_src_checksum="/p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_BB/postpro/CMORized/checksum_ALP-3.sha256"
dir_target_checksum="/mnt/CORDEX_FPS_CPCM"

# CORDEX-FPSCONV VL v2017-11-30 complete
# many vars are processed at shorter output interval, others are not processed available at all
# pull/push what is needed or less, not more
# check whather a certain variable exists or not beforehand
declare -a list_aggr_var=( 1hr/tas 1hr/pr )
#declare -a list_aggr_var=( 1hr/tas 1hr/ts 1hr/pr 1hr/wsgsmax 1hr/uas 1hr/ua100m 1hr/vas 1hr/va100m 1hr/huss 1hr/rsds 1hr/rlds 1hr/hfls 1hr/hfss 1hr/rsus 1hr/rlus 1hr/evspsbl 1hr/evspsblpot 1hr/clt 1hr/mrros 1hr/mrro 1hr/prw 1hr/clwvi 1hr/clivi 1hr/zmla 1hr/clgvi 1hr/clhvi \
#3hr/ua1000 3hr/ua925 3hr/ua850 3hr/ua700 3hr/ua500 3hr/ua200 3hr/va1000 3hr/va925 3hr/va850 3hr/va700 3hr/va500 3hr/va200 3hr/wa1000 3hr/wa925 3hr/wa850 3hr/wa700 3hr/wa500 3hr/wa200 3hr/hus1000 3hr/hus925 3hr/hus850 3hr/hus700 3hr/hus500 3hr/hus200 3hr/ta1000 3hr/ta925 3hr/ta850 3hr/ta700 3hr/ta500 3hr/ta200 \
#6hr/zg1000 6hr/zg925 6hr/zg850 6hr/zg700 6hr/zg500 6hr/zg200 6hr/ps 6hr/psl \
#day/tasmax day/tasmin day/wsgsmax day/sfcWindmax day/snc day/snd day/mrsol day/mrso day/wsgmax100m day/ic_lightning day/cg_lighning day/total_lightning day/cape day/cin \
#fx/areacella fx/orog fx/sftlf fx/sftgif fx/mrsofc fx/rootd )

# get the checksum
cd $HOME
scp -p -i ~/.ssh/id_ed25519_JSC_HPCsys_goergen1_transfer goergen1@jureca.fz-juelich.de:${dir_fn_src_checksum} .

# find out how much there is to transfer
ssh -i ~/.ssh/id_ed25519_JSC_HPCsys_goergen1_transfer goergen1@jureca.fz-juelich.de find $dir_src -type f -wholename "*.nc" | sort | wc -l

# create top-level top-level directory
mkdir -p $dir_target
cd $dir_target && pwd

# go through the variables
for ifile in ${!list_aggr_var[@]}
do
  echo "index = ${ifile}, value = ${list_aggr_var[$ifile]}"
  # check whether the source dir exists at all
  # if it exists, then there should be netCDF data contained to be uploaded
  if [[ `ssh -i ~/.ssh/id_ed25519_JSC_HPCsys_goergen1_transfer goergen1@jureca.fz-juelich.de test -d ${dir_src}/${list_aggr_var[$ifile]} && echo exists` ]] ; then
    # create new dir, cd into, scp files
    mkdir -p ${dir_target}/${list_aggr_var[$ifile]}
    cd ${dir_target}/${list_aggr_var[$ifile]}
    #scp -p -i ~/.ssh/id_ed25519_JSC_HPCsys_goergen1_transfer goergen1@jureca.fz-juelich.de:${dir_src}/${list_aggr_var[$ifile]}/*.nc .
    rsync -hanvP -e "ssh -i $HOME/.ssh/id_ed25519_JSC_HPCsys_goergen1_transfer" goergen1@jureca.fz-juelich.de:${dir_src}/${list_aggr_var[$ifile]}/ ./ # further options: c takes longer, z (compressed netCDF)
    if [[ $? -ne 0 ]] ; then
      echo "WARNING: scp/rsync failed"
    fi
    # just get the relevent portion of checksum file we are working on
    # issues warning message anyhow
    grep "${list_aggr_var[$ifile]}/" $HOME/$(basename $dir_fn_src_checksum) > ${HOME}/to_check.txt
    cd $dir_target_checksum
    sha256sum --check ${HOME}/to_check.txt
    if [[ $? -ne 0 ]] ; then
      echo "WARNING: sha256sum failed"
    fi
  fi
done

# check number of transferred files
cd $dir_target
echo "number of files in new target dir ="
find . -type f -wholename "*.nc" | sort | wc -l

exit 0
