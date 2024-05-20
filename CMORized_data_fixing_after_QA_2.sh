#!/bin/bash

# nohup time ./CMORized_data_fixing_after_QA_2.sh >> CMORized_data_fixing_after_QA_2.log &
# this is a quick fix for the fx files which are missing a variable "char rotated_pole"
# the python script from heimo would need to be adjusted

date
echo $0
echo $USER
uname -a

declare -a exp_list=( 
#"/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/ALP-3/FZJ/ECMWF-ERAINT/evaluation/r1i1p1/FZJ-WRF381BB/fpsconv-x1n2-v1"
#"/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/EUR-15/FZJ/ECMWF-ERAINT/evaluation/r1i1p1/FZJ-WRF381BB/fpsconv-x0n1-v1"
#"/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/ALP-3/FZJ-IDL/SMHI-EC-EARTH/historical/r12i1p1/FZJ-IDL-WRF381CA/fpsconv-x1n2-v1"
#"/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/EUR-15/FZJ-IDL/SMHI-EC-EARTH/historical/r12i1p1/FZJ-IDL-WRF381CA/fpsconv-x0n1-v1"
#"/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/ALP-3/FZJ-IDL/SMHI-EC-EARTH/historical/r12i1p1/FZJ-IDL-WRF381DA/fpsconv-x1n2-v1"
#"/mnt/CORDEX_FPS_CPCM/CORDEX-FPSCONV/output/EUR-15/FZJ-IDL/SMHI-EC-EARTH/historical/r12i1p1/FZJ-IDL-WRF381DA/fpsconv-x0n1-v1"
)

for i_pn in "${exp_list[@]}" ; do
  echo $i_pn
  pn_array=(`find ${i_pn} -name "*.nc" -exec dirname {} \; | sort | uniq | grep fx`)
  for i_exp_var_pn in "${pn_array[@]}" ; do
    chmod u+w $i_exp_var_pn
    cd $i_exp_var_pn && pwd
    for i_fn in *.nc ; do
      echo $i_fn
      # add a missing var to fx files
      ncap2 -h -s 'rotated_pole=0' $i_fn ~/tmp1.nc
      ncap2 -h -s 'rotated_pole=rotated_pole.convert(NC_CHAR)' ~/tmp1.nc ~/tmp2.nc
      ncatted -a long_name,rotated_pole,a,c,"coordinates of the rotated North Pole" -h ~/tmp2.nc ~/tmp3.nc
      ncatted -a grid_mapping_name,rotated_pole,a,c,"rotated_latitude_longitude" -h ~/tmp3.nc ~/tmp4.nc
      ncatted -a grid_north_pole_latitude,rotated_pole,a,f,39.25 -h ~/tmp4.nc ~/tmp5.nc
      ncatted -a grid_north_pole_longitude,rotated_pole,a,f,-162. -h ~/tmp5.nc ~/tmp6.nc
      chmod u+w $i_fn
      rm -v $i_fn
      mv -v ~/tmp6.nc $i_fn
      rm -v ~/tmp?.nc
    done
    chmod -R a-w $i_exp_var_pn
  done
done

exit 0

