#!/bin/bash

#SBATCH --job-name="pCMORizer_check_integrity"
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

#-------------------------------------------------------------------------------
# 
# CURRENT CODE OWNER(S): 
#   Klaus GOERGEN (KGo), FZJ/IBG-3, k.goergen@fz-juelich.de
#
# VERSION (version control via git): 
#   2023-11-11 -- revision details are in the git log
#
# STATUS:
#   Release
#
# PURPOSE: 
#   - After CMORization: checks integrity of highest temporal resolution (e.g., 
#     1hr) output to clarify whether (i) the CMORizer ran OK, (ii) there are 
#     no I/O issues, (iii) the filesystem has not caused any issues with the 
#     data.
#   - This shall help to avoid errors downstream during aggregation, QA 
#     checking, and usage of the data.
#
# USAGE: 
#   - sbatch ./$0 <root_dir(s)_of_CMORized_output>
#   - Adjust sbatch settings to your local HPC environment.
#
# EXAMPLES: 
#   Typical usage is on HPC compute node.
#   - sbatch ./pCMORizer_check_integrity.sh /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized/CORDEX-FPSCONV
#   - sbatch ./pCMORizer_check_integrity.sh /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized/CORDEX-FPSCONV /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized_2_was_thought_to_be_OK/CORDEX-FPSCONV
#
# EXAMPLE FOR TESTING: 
#   Command line test also possible.
#   - ./$0 <root_dir(s)_of_CMORized_output>
#   - nohup ./$0 <root_dir(s)_of_CMORized_output> > output.txt &
#   - See manual adjustments in the code to run on front node for single 
#     variable (and year) tests, mainly adjust the number of processes and
#     most likely restruct the search pattern of find to work on less input
#     files.
#
# PROCEDURE: 
#   - Designed to run on HPC compute node to allow for fast processing. Uses
#     ntasks in parallel; loads all threads ("&") and then waits ("wait) until 
#     current batch of processes (e.g., 128 ) has finished. 
#   - Total runtime is below 30min to check 10 years of 1hr CMORized 
#     CORDEX-FPSCONV output with 128 cores on JSC JURECA-DC HPC.
#   - By reading ifiles and calculating statistics (cdo info) the tool detects
#     (i) unreasonable data), (ii) fields with NaN values (missing values in the
#     data due to pressure level intersections with topography are not detected
#     as this is a feature, (iii) issues with the internal file structure (if a 
#     HDF5 based netCDF4-classic file has is internally flawed for whatever 
#     reason the tool issues a warning in the log. 
#   - Adjust the RegEx in the for loop according to your needs. Otherwise all 
#     netCDF files are checked.
#   - Manually check the standard output on the command line (test example) or
#     check the file with the scheduler output for results. Optionally check 
#     statistics for each timestep for each input file analysed.
#
# INPUTS / ADJUSTMENTS:
#   - Parm1: CMORized data path, top level directory. Multiple paths may be 
#            specified after each other. They must not lead to CMORized data
#            of the same experiment. Otherwise the same statistics files of
#            the two experiments overwrite each other in the dirLog/ output.
#            directory.
#   - sbatch settings: Adjust to your needs. The job-name is used to to name
#            the output log file and the directory with the individual 
#            statistics files.
#   - Adjust the compute environment settings to your needs.
#
# OUTPUTS:
#   - Output1: When started with sbatch (normal operation) main output is in
#              scheduler log file, naming is according to script name and slurm
#              job ID ID. Output contains of (i) all files worked on, including 
#              error messages if a bad file has been encountered, (ii) a listing 
#              which file has how many NaN fields, (iii) a detailed listing per 
#              file where NaNs occur (which timestep).
#   - Output2: If there are file strzctural errors (due to filesystem?) then 
#              also the *err slurm output contains messages.
#   - Output3: A newly created directory in tune run-directory with a naming 
#              according to script name and slurm job-ID (or date alternatively)
#              which contains per checked input file a ASCII statistics file
#              of the cdo info output.
#
# REQUIREMENTS: 
#   - bash shell
#   - CDO -- tested with CDO v2.1.1
# 
# RESTRICTION:
#   - Adjust the NaN search pattern to your needs manually in the code. If the 
#     size of the model domain changes, then the pattern needs adjustment.
#   - Assumes slurm is used as scheduler if used on a HPC compute node.
#   - If the CMORizer fails due to corrupted netCDF/HDF5 input file, this tool
#     only helps to find the first occurence of the error. I.e. if there is a 
#     real issue with the inputs in the same year, then this will need more 
#     intensive checking.
# 
#-------------------------------------------------------------------------------

# load environment; needs adjustment for different HPC environments
source loadenv.JURECA-DC_2023_Intel-PSMPI.ini

# multiple input paths are possible, if this is however the same experiment
# (like different CMORization attempts) then there might be naming conflicts
# with the stats files; but different domains of the same experiment might be 
# dealt with at the same time
inFiles=$@
echo "${inFiles[@]}"
declare -i tmp_file_counter=0

# set MAX_PARALLEL to the number of parallel threads, 1-n
# if there are less files to check than MAX_PARALLEL, the number of files i
# determines the number of parallel threads; OS usually takes care of the load
# distribution on several cores
#declare -i MAX_PARALLEL=128  # for front node testing
declare -i MAX_PARALLEL=${SLURM_NTASKS}
echo "MAX_PARALLEL: $MAX_PARALLEL"
declare -i tmp_parallel_counter=0

# a tmp direcory contains statistics outpouts from cdo info for each input file
#dirLog=$(pwd)/tmp.pCMORizer_check_integrity.logs.$(date +%Y%m%d%H%M%S)  # for front node testing
dirLog=$(pwd)/tmp.${SLURM_JOB_NAME}.logs.${SLURM_JOB_ID}
mkdir -p $dirLog

# cdo info is the main functionality
# if it fails, the script continues, a warning is issued withan indication of
# the file that causes an issue
function calcStats(){
  cdo info $1 > "${dirLog}/$(basename $1).cdo_info_ouput.txt"
  if [[ $? -ne 0 ]] ; then
    echo "WARNING: cdo failed (perhaps corrupted file): $1"
    return 1
  fi
}

# loop over all files and initiate the checking
# start as many cdo info jobs per single CPU core until all cores are filed
# or the file list reaches its end, whatever occurs first
# wait until all checks have finishedm, then continue to the next set of files
# to check
date
for inFile in $(find ${inFiles[@]} -type f -wholename "*.nc" | sort)
# for testing only with a subset of the data
#for inFile in $(find ${inFiles[@]} -type f -wholename "*1hr/*/*.nc" | sort)
#for inFile in $(find ${inFiles[@]} -type f -wholename "*1hr/*/*.nc" | sort | grep -E "rsds_")
#for inFile in $(find ${inFiles[@]} -type f -wholename "*1hr/*/*.nc" | sort | grep -E "2002" | grep -E "rsds_")
do
  echo "DEBUG: inFile ${inFile}"
  calcStats ${inFile} &
  (( tmp_parallel_counter++ ))
  (( tmp_file_counter++ ))
  if [ $tmp_parallel_counter -ge $MAX_PARALLEL ]; then
    wait
    tmp_parallel_counter=0
    date
  fi
done

wait

# after all statistics have been calulated, grep finds the missing values from
# the statistics files and writes to standard out, which is the slurm job log
# files
# grep does not pick up missing values in fields where pressure level is below 
# the land surface
echo "========================================"
echo "total number of nc files checked = ${tmp_file_counter}"
echo "========================================"
echo "files which have NaN and for how many fields (timesteps):"
for inFile in $dirLog/* ; do 
  n_nan=$(grep -i -c "197400  197400 :                     nan" $inFile)  # ALP-3
  #n_nan=$(grep -i -c "111870  111870 :                     nan" $inFile)  # EUR-15
  if [[ $n_nan -ge 1 ]] ; then
    echo "${inFile}: ${n_nan}"
  fi
done
echo "========================================"
echo "per file indicate where (which timstep) exactly the NaNs occur"
for inFile in $dirLog/* ; do 
  echo "--------------------"
  echo $inFile
  grep -i "197400  197400 :                     nan" $inFile  # ALP-3
  #grep -i "111870  111870 :                     nan" $inFile  # EUR-15
done

exit 0
