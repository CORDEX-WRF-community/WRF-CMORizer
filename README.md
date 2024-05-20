# pCMORizer

2024-05-20

The pCMORizer is a free open-source software tool to transfer and postprocess raw (regional climate) model ouputs into netCDF files, which comply to the [CMOR](https://cmor.llnl.gov) data standard as used by the [WCRP CORDEX](https://cordex.org) project. The goal is to make the model data standard compliant ([Standards section](#ref_standards)) to be efficiently (i) shared via data servers, such as the [ESGF](https://esgf-data.dkrz.de/projects/esgf-dkrz/) data nodes, and (ii) used in analyses. The tool was developed for the WRF RCM but can also be used with other model outputs after a few adjustments.

The main features of the pCMORizer are (see also [Additional information section](#ref_addInfo)):
- Tool is run *parallel* (this is what the "p" stands for) (Step 1 and Step 2) (multi-year, multi-domain experiment can be fully processed in less than 8h wall clock time) (preferrably on compute node of HPC cluster: e.g., all required variables per simulated year and model domain are run on a single node, MPI parallel).
- Few I/O operations (data is read, processed, and written to CMOR compliant files).
- Low memory footprint (data is read timespep by timestep and sorted into CMOR- / CORDEX data standard compliant, compressed netCDF files).
- Lightweight, highly versatile in terms of operating modes (flexible timespans and variable lists) and workflow embedding, either controlled via Fortran code or more via the run-control scripts, easily adjustable and expandable.

<!-- 
Fortran-based software tool with some ancilliary bash scripts that also use [cdo](https://code.mpimet.mpg.de/projects/cdo/). pCMORizer.f90 and companion tools are used
according to the [CORDEX-CMIP3 archive specification](http://is-enes-data.github.io/cordex_archive_specifications.pdf) as part of the [CORDEX WCRP project](https://cordex.org/experiment-guidelines/how-to-submit-data-rcms/). 
-->

Attention: there is a quick fix for the fx files which need to be applied: `CMORized_data_fixing_after_QA_2.sh`

## Usage<a name="ref_usage"></a>

The pCMORizer can be used in many different "operating modes", i.e., serially or parallel, controlled via the namelist files and the Fortran code, or in conjunction with the auxilliary runcontrol scripts, which also do some file handling.  The most efficient variant dealt with here makes use of 4 different exacutables of `pCMORizer.f90`, each tailored to a specific set of variables to process and uses MPI parallelism to substantially speedup processing.

All aspects of using and adjusting the pCMORizer to our needs are explained below.

### Step 0: Setup

#### Obtain the code

Get the current pCMORizer version (will create the `pCMORizer/` directory):

```
git clone https://icg4geo.icg.kfa-juelich.de/ExternalRepos/pCMORizer.git
```

For bash shell, set the pCMORizer root directory of the tool, this variable is just used in this manual to make things clear:

```
export pCMORizer_DIR=$(pwd)/pCMORizer
```

#### Set the computational environment

Load the compile-time and run-time environment, also loaded by the runcontrol scripts during runtime; the default environment used here is Intel based on the [JURECA-DC HPC](https://www.fz-juelich.de/en/ias/jsc/systems/supercomputers/jureca): 

```
cd $pCMORizer_DIR
source loadenv.JURECA-DC_2023_Intel-PSMPI.ini
```

No further environmental variables needed; all envrionment settings needed during runtime are set within the respective scripts.

In a nutshell JURECA-DC has compute nodes with x2 CPUs with 64 physical CPU cores each (without hyperthreading, which would double the amount of processes to 256, this is not used). 

If you are not on JURECA-DC, make sure the following tools and libraries are available; any fairly recent version as of 2023 will do:

General Linux operating environment, needed with all steps:
- bash (including: uuidgen, sha256sum)
- sbatch/srun HPC scheduling environment

Needed for Step 1:
- Intel Fortran compiler (a GCC built is in principle also possible)
- MPI implementation, e.g., OpenMPI 
- netCDF library
- make

Needed for Step 2:
- cdo
- nco

Needed for Step 3:
- Python3 (including: numpy, xarray, f90nml)

Useful at any time:
- ncdump
- ncview

#### Manual code adjustments and compilation

The `pCMORizer.f90` is the main code of the tool. Some things need to be manually adjusted in `pCMORizer.f90`.

To be adjusted once, depending on the (i) model config during simulation, (ii) filename patterns, (iii) storage pattern:
- T00 (290K or 300K, check your `Registry.COMMON`), line 1855
- Depths of the land surface model, lines 309 and 1774f -- determined by the Noah LSM
- Type of pressure calc, line 366 -- default is as agreed upon by CORDEX-FPSCONV
- Vertical interpolation, line 369 -- dfeault is as agreed upon by CORDEX-FPSCONV
- `inputtimesteptruncate` (line 361) handling of the time info in the files, 00-23 (.F.), or 00-00 (.T.), i.e., does the model output contain data until 23UTC or is the 00UTC of the subsequent day also contained
- Pattern matching of the filenames if they are not WRF default (`wrfout`, `wrfxtrm`, `wrfpress`), filename pattern with find, lines 623ff

The following should be done more elegantly. As the code offers many operating options as a development legacy, the following needs to be adjusted for each type of variable set (see below), then the code needs to be compiled with the respective settings; i.e. one ends up with 4 executables, one for each variable type and a specific name of the executable. Per type of variable set / for each executable change this:
- The maximum number of variables to work on, line 42
- Hourly vs daily (for tasmin/tasmax) data, lines 582f
- Variable namelist to be used, lines 691-694

```
cd $pCMORizer_DIR
vim pCMORizer.f90
```

#### Compile pCMORizer.f90

Adjust the makefile (`FCFLAGS`, `LDFLAGS`, lines 4 and 5) (and compile flags if you are not using Intel, lines 2 and 3), and compile. 

```
cd $pCMORizer_DIR
vim Makefile
```

Compile the code 4 times, see necessary edits above:

```
cd $pCMORizer_DIR
make veryclean
vim pCMORizer.f90
make
rm *.mod
mv pCMORizer <exe_name>
```

Currently 4 executable are predefined: 
```
pCMORizer_std_sfc_39vars.exe
pCMORizer_special_2vars.exe
pCMORizer_std_presslev_36vars.exe
pCMORizer_std_minmax_2vars.exe
```

The code only needs to be adjusted for the different types of variable sets to be processed, and for very basic model configurations. It does not have to be adjusted for years, or experiments (e.g., evaluation vs historical vs rcp).

### Step 1: Run the pCMORizer, Fortran tool

**In Step 1, raw model output is transferred to CMORized, fully standard compliant files, with correct meta-data, DRS, filename conventions. After this step, files do not need to be touched anymore. With the different CORDEX archive protocols, data is needed at specific temporal intervals with a aspecific aggregation, (cell_method time:point or time:mean). In Step 1, all variables are processed at the highest available output frequency, whether this is required by the protocol or not. The reasoning is that after this step, data might be archived or erased and only data of Step 1 are kept for further analysis and derivation of data at lower frequency.**

#### Configuration, adjust to simulation experiment

##### Namelist 1: Meta-data, timespan, processing settings, directories

The main `runctrl.current.nml` namelist configures the global attributes and defines all directories and the year range. The runctrl files are documented and determine all major metadata and domain settings and also point to static files which need to be provided for the processing (`geo_em` files).

The current configuration uses '0000' as start and end year, this means the file name pattern used by the `find` command in `pCMORizer.f90` does not apply any years as a pattern for the input files.

Per experiment, two `runctrl.current.nml` files are provided one for each model domain (EUR-15 and ALP-3 for CORDEX-FPSCONV); see `d01` and `d02`.

In addition, for each experiment, B, C and D runs with the WRF subgroup in CORDEX-FPSCONV, a set of runcontrol file sexists, in line with the archive protocol and ESGF requirements.

If the experiment does not change, no changes are needed here. If a new experiment is to be processd, adjust these files accordingly. *These files are the only location that needs an adjustment for a new experiment.*

The following configurations and setups are currently provided for usage and as examples:
```
runctrl.current.nml_template_d01_BB
runctrl.current.nml_template_d01_CA
runctrl.current.nml_template_d01_DA
runctrl.current.nml_template_d02_BB
runctrl.current.nml_template_d02_CA
runctrl.current.nml_template_d02_DA
```

If necessary edit them:
```
cd $pCMORizer_DIR
vim -R -o runctrl.current.nml_template_d01_DA runctrl.current.nml_template_d02_DA
```

##### Namelist 2: Variable lists

The variable namelists do not need to be changed usually. There are 4 variable lists, one per "variable group" or "type" of variable. If a variable is not needed, this variable has to be erased form the list; the functionality which could be used before to switch a variable ON and OFF is essentially deprecated with the MPI version of the tool (although the switch as such still works).

If the protocol does not change, if no new variables are required, or some variables are not required anymore, no changes are needed here.

```
vim runctrl.vars.std_sfc.nml # 39 vars
vim runctrl.vars.std_minmax.nml # 2 vars
vim runctrl.vars.std_presslev.nml # 36 vars
vim runctrl.vars.special.nml # 2 vars
```

The content of the variable lists depends on the simulation experiment configuration, the model's capabilities and the archive protocol, or the variable list (VL) as requested by the CORDEX experiment.

#### Running the pCMORizer

**The currently supported operating mode, which allows for a very fast operation of the pCMORizer, assumes that the tool is run on a compute node of Linux cluster that has access to a single or some years of raw model output that shall be CMORized after the model run has well finished. This is the CORDEX-FPSCONV usage case.** The tool can also be run for any timespan other than years, and any temporal coverage of input data, i.e., also on data just stored by the model run or already archived raw data; see Additional Information below. In summary, there are many ways of running the pCMORizer tool, depending on processing setup, urgency, amount of data, available resources. 

The run control scripts (i) handle the data and (ii) start the CMORizer, making sure via srun options, that there is no overlap as `srun` submitted two times on the same compute node. Here multiple manual adjustments are needed. See comments in the code (sbatch options, wallclock times, see performance section; linking the correct variable and runctrl namelists).

The most common modes of operation with the MPI default variant of the pCMORizer are decribed below.

##### Operating mode 1: All types of variable set, single domain, one year per node (default)

Within the pCMORizer working directory, adjust the run-control script (see comments in the script):

```
cd $pCMORizer_DIR
vim pCMORizer_runctrl_allVarTypes_singleDomain.sh
```

Launch for single year, here an example for 1998 and domain `d01`; this will process all 4 types of the variables sets in parallel, each with its own executable, one variable per CPU core, using a shared raw model data input file and writing to multiple CMORized output files, one per CPU core:

```
cd`$pCMORizer_DIR
export Y=1998 && export DOM=d01 && sbatch --export=ALL,Y=$Y,DOM=$DOM --job-name=pCMORizer$Y$DOM pCMORizer_runctrl_allvars_single_domain.sh
```

The run-control script generates a domain directory, e.g. `d01/`, inside the `$pCMORizer_DIR` and within one subdirectory per year. All relevant namelist files, executables and input files are linked into to these directories. The output (fully standard compliant DRS-based directory tree) is written as specified in the `runctrl.current.nml` namelist. The script does not clean up after itself, hence the working directory needs to be erased manually afterwards.

The runcontrol script can recursively launch itself once finished, hence working on year by year.

Alternatively, if all data is available to the compute nodes, launch the processing of several years and both domain for all variables. For a 10 year time slice CORDEx-FPSCONV simulation this results in 10 years x 2 domains nodes; with a single pass of the tool, all CMORized variables are generated:

```
cd $pCMORizer_DIR
vim pCMORizer_launcher_multipleYears1yPerCompNode.sh
./pCMORizer_launcher_multipleYears1yPerCompNode.sh
```

> Performance:<br>
> - Wall clock runtime on JURECA-DC with 128 cores per compute node, 79 variables (most of CORDEX-FPSCONV variable list), single year, hourly data processing: ALP-3 domain: **about 5-6.5h** with about **42GB** on average RAM usage; EUR-15 domain: **about 3.5h**.<br>
> - Because the different types of variables have different procssing times (CAPE and CIN take the longest, pressure level the 2nd longest, then surface variables, and min/max variables, where no processing is to be done and which are on a daily basis only), the CPU load throughout the runtime varies and becomes somewhat inefficient towards the end, yet still CPU load is about 88% on average.<br>
> - To shorten the processing time, CAPE and CIN are calculated on a monthly basis and then concatenated to yearly files.
> - The runtime varies depending on from where data are read on JURECA-DC, scratch-based data I/O is about 1h faster than largedata I/O.

In general the pCMORizer can work on any data granularity (yearly, monthly, sub-monthly) and input data time span. Data are always sorted into the output file to where they belong. If a shorter timespan is desired then the input data, then the input data outside this timespan is ignored.

*Operating mode 1 is considered the most efficient for v1.0.0.*

##### Operating mode 2: Single type of variable set, two domains, one year per node 

Within the pCMORizer working directory, adjust the run-control script (see comments in the script):

```
cd $pCMORizer_DIR
vim pCMORizer_runctrl_singleVarTypes_twoDomains.sh
```

Launch for single year, here an example for 1998 and domain `d01` and `d02`; this will process a single type of variable set in parallel, based on a single executable, one variable per CPU core, using a shared raw model data input file and writing to multiple CMORized output files, one per CPU core:

```
cd $pCMORizer_DIR
export Y=1998 && sbatch --export=ALL,Y=$Y --job-name=pCMORizer$Y pCMORizer_runctrl_singleVarTypes_twoDomains.sh
```

As above, the compute environment is loaded automatically. Within the current working directory a `<year>/` directory is created; all scripts and namelists will be contained in the year-directory. The input data is symbolically linked into this directory; after succesful processing, this directory can be erased, manually. The postprocessed / CMORized data are stored as specified (in their own separate directory). (The input directory could also be anywhere.) All needed input files (runcontrol namelist, variable namlist, executable, input files are in the run-dir, i.e., the year directory.)

Again, the launcher can launch multiple years on seperate nodes to each other at the same time.

```
cd $pCMORizer_DIR
vim pCMORizer_launcher_multipleYears1yPerCompNode.sh
./pCMORizer_launcher_multipleYears1yPerCompNode.sh
```

> Performance:<br>
> On [JURECA-DC compute nodes](https://apps.fz-juelich.de/jsc/hps/jureca/index.html) in combination with [JUST](https://www.fz-juelich.de/en/ias/jsc/systems/storage-systems/just) [filesystems](https://apps.fz-juelich.de/jsc/hps/jureca/filesystems.html) the following performance can be reached.<br>
> Use the information below to set your sbatch wallclock times.<br>
Setup: Input data is WRF raw model outputs (`wrfout`), compressed; one netCDF file per day; one compute node processes one year at a time; two domains are processed concurrently per node (EUR-15, ALP-3); each CPU-core processes one variable; different variable types are processed in seperate, manually started and adjusted processing stepsi; ten compute nodes are run via 10 sbatch commands concurretlyi; Intel compile, with `-O3`; logging is OFF. See [Usage section](#ref_usage) above for details; given are the approximate average runtimes accross the 10 compute jobs; all model outputs are provided on spinning disk at the same time.<br>
> - 39 near-surface variables: about 5h25min (single month test: 30min suffices)
> - 36 pressure level variables: about 6h30min
> - 2 min/max variables: about 5min
> - 2 variables with more extensive diagnostics calculations:
>
> I.e., using 10 compute nodes with exclusive allocation, a decade long km-scale experiment can be CMORized in less than 24h wall clock time.

### Step 2: Temporal aggregation, using shell scripting with `cdo` and `nco`

Based on data from Step 1 (e.g., at 1hr intervals), the temporal aggregation as specified in the CORDEX-FPSCONV variable list is done using cdo and nco in a bash script, the script starts as many single threaded processing streams until the `nproc_max` is reached; then it waits until processing has finished to start the next batch; the outputs are standard complant (metadata, time-information, DRS, filenames, etc.) and stored at their respective location in the CMOR directory tree. Only the input directory, number of concurrent processes (can start with single process) and compute environment need to be specified. If an error occurs script exits:

```
cd $pCMORizer_DIR
vim pCMORizer_aggregate_in_time.sh
sbatch ./pCMORizer_aggregate_in_time.sh
```

> Performance:<br>
> JURECA-DC compute node, 10 years, 43 vars, single node takes about 20min. 

If the variable list changes, the script needs adjustment; primarily variable lists inside the script need to be expanded.

### Step 3: Create fx files

In order to create the fx files, a Python3 script is run. Due to the availability of the Python modules it might be used on any Linux machine. The input is the `runctrl.current.nml` namelist files, this controls the operation. Additionally the static input data are needed, like `geo_em` files for the WRF RCM.

```
cd $pCMORizer_DIR
vim pCMORizer_create_fixed_fields.py
python3 pCMORizer_create_fixed_fields.py
```

The output is the fx-files with their correct DRS pathname and filename and meta-data. The root-dir is the `$pCMORizer_DIR`. The files needs to be manually integrated into the main CMORized data directory tree.

### Step 4: Run auxilliary tasks: integrity check and checksums

#### Integrity check

Run this after CMORization in Step 1: checks integrity of highest temporal resolution (e.g., 1hr) output to clarify whether (i) the CMORizer ran OK, (ii) there are no I/O issues, (iii) the filesystem has not caused any issues with the data. See the extensive preamble to the code with further explanations. The log files of the tool easily indicate where potential issues might occur. The tool runs in parallel, on as many resources as available on Linux cluster compute node.

```
cd $pCMORizer_DIR
vim pCMORizer_check_integrity.sh
sbatch ./pCMORizer_check_integrity.sh /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized/CORDEX-FPSCONV`
```

The DRS-based directory tree of the CMORized data remains untouched. 

> Performance:<br>
> About 20min wall clock time with 128 cores on a single compute node for all 1hr FPSCONV variables for 10 years.

Background: It happened that filesystem issues either during writing of original raw model outputs or during long-term storage on spinnning disk lead to corrupted files netCDF, i.e. underlying HDF5 files. This is a major problem and may lead inevitably to processing errors. As this is usually not to be assumed, this check is done after the CMORization; if there is an issue with input files, the processing might have to be repeated, e.g., with restaged input data.

#### Checksums

If data is transferred by whatever means, it might be helpful to have checksums available to ensure data integrity; one directory level above the root dir of the DRS diretory tree a file is created with checksums for each individual netCDF file in the DRS tree. After file transfers this checksum inventory might be used to ensure the exact replication of directory tree and data files. The tool operates like the integrity checker in parallel on a Linux cluster compute node:

```
cd $pCMORizer_DIR
vim pCMORizer_checksums_centralized.sh
sbatch ./pCMORizer_checksums_centralized.sh /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized  # 1 level above CORDEX-FPSCONV/
```

> Performance:<br>
> About 10min wall clock time for 1240 files with 128 core on a single compute node.

An alternative tool exists which creates a checksum file per data directory, if this is preferred. Because it clutters the directory tree it should not be the preferred solution:

```
cd $pCMORizer_DIR
vim pCMORizer_checksums_distributed.sh
sbatch ./pCMORizer_checksums_distributed.sh /p/scratch/cjjsc39/goergen1/sim/tmp_FPSCONV/tmp_DA/postpro/CMORized  # 1 level above CORDEX-FPSCONV/
```

### Step 5: Run the CMORized data through the QA Checker

See [Andreas Doblers repo](https://github.com/doblerone/QA-DKRZ_FPSCONV), with DKRZ QA Checker, adjusted to CORDEx-FPSCONV.

Data processed with pCMORizer v1.0.0 have passed the QA Checker.

**More information to come.**

## Compatibility and standards compliance<a name="ref_standards"></a>

pCMORizer was originally developed to cmorize raw [WRF RCM](https://www2.mmm.ucar.edu/wrf/users/) outputs. It seems generic enough in terms of structure and configuration to also be used for other RCMs, LSMs, HMs. 

At this point, the CORDEX-CMIP5 and the CORDEX-FPSCONV experiments archive specifications and variable lists (VLs), in their (most current) versions as of 2022-09-25, are supported, see the [references](#ref_refs) and linked documents therein.

The tool also covers CORDEX-CMIP6 variable standards.

Some (for the majority of the users not so relevant) variables cannot be derived (yet), see [ToDos](#ref_todos).

The [fps_convection_variables_vfinal_2021update_with_addons.xlsx](fps_convection_variables_vfinal_2021update_with_addons.xlsx) table is the final VL for FPSCONV with addons from the CORDEX-CMIP5 and CORDEX-CMIP6 VLs and gives and overview which variables are covered by the pCMORizer and other tools. If you add new variables, please to not only edit the variable namelists but also update this table.

Adjustments have been made to adhere to QA Checker standards and requirements.

## ToDos<a name="ref_todos"></a>

The current version is in the VERSION.txt file. pCMORizer is a release version which allows to process data for later-on ESGF staging. Not all VL variables need to be available to qualify for ESGF staging with RCM ensemble.

Ongoing modfications:

- Adjust for additional RCMs, LSMs, HMs based on the exsting structure and configuration namelists

Missing variables in v1.0.0 according to CORDEX-FPSCONV VL, additional processing procedures would need to be implemented in `pCMORizer.f90` (mainly for ua100m, va100m); other variables need different model configuration (evspsblpot, lightning, cl):

Fluxes:

- `evspsblpot`

Wind speed:

- `ua100m`
- `va100m`
- `wsgsmax100m` 

fx-fields:

- `mrsofc` 
- `rootd`

Lightning parametrisation

- `ic_lightning` 
- `cg_lighning`
- `total_lightning` 

Microphys:

- `clgvi`
- `clhvi`

Additional CMIP6 VL, excerpt:

- `clh`
- `clm`
- `cll`

Nice to have:

- More efficient calculation of more complex diagnostics, such as `cape` and `cin`
- Compilation using GCC and OpenMPI, fix of some compile time warnings, which lead to runtime abort with GCC
- Improved chunking, although a series of tests has shown no improvement when netCDF file chunking is altered to what the netCDF library does as a default
- All temporal aggregations and regriddings as required by the complete CORDEX-CMIP3 standard

Future Todos:

- Extension to be used with the CORDEX-CMIP6 experiment

## Documentation

This README.md file and comments in the files are the only documentation.

## Release

The current release version as of 2023-11-22 is v1.0.0 (see git tags).

Written and tested by in the CONTRIBUTORS.txt file.

pCMORizer is released under the MIT License.

For details and restrictions, please read the LICENSE.txt file.

## References<a name="ref_refs"></a>

- [CF-conventions](https://cfconventions.org/)
- [CORDEX-CMIP3 downscaling experiment guideleines documents](https://cordex.org/experiment-guidelines/how-to-submit-data-rcms/)
- [Specific reuqirements of the CORDEX FPSCONV initiative](https://www.hymex.org/cordexfps-convection/wiki/doku.php?id=protocol)
- [DKRZ QA checker](https://github.com/IS-ENES-Data/QA-DKRZ)
- [ESGF infrastructure data node for public dissemination](https://esgf-data.dkrz.de/projects/esgf-dkrz/)

## Links<a name="ref_links"></a>

Similar tools and resources:

- [https://github.com/C2SM-RCM/CCLM2CMOR](https://github.com/C2SM-RCM/CCLM2CMOR)

## Additional information<a name="ref_addInfo"></a>

The information below might be a bit outdated but is provided for completeness nevertheless. In order to use the pCMORizer, the above information is sufficient.

### Motivation, features, design principles

At an early stage of the CORDEX-CMIP3 experiment in 2013, there were very few tools available, which would ingest raw RCM output and generate CMORized simulation results for dissemination and/or efficient in-house usage or archival, meeting FAIR principles. Existing solutions at the time often consisted of complex, difficult-to-maintain script-based processing chains, that resulted in lots of time-consuming I/O operations, issues when it came to the generation of more complex variables and diagnostics (as script-based tools could not be used for this purpose any more), or large memory footprints. 

In this context the pCMORizer.f90 was orginally developed as a purely serial tool for WRF RCM outputs. The orignal version was expanded by the community now and then, partly adopted to specific workflows, and further customized. Full operational capability, was only reached step by step, triggered by archival, dissemination, and publication needs. Today, a range of CMORizer tools exist for various RCMs, see the [Links section](#ref_links). The pCMORizer.f90 has undergone three major development stages (or variants), only the most current is maintained and documented here. 

### Maintained MPI-parallel variant

- Two-staged processing: Fortran tool does initial CMORization step, bash+cdo does temporal aggregation
- MPI-parallel (thanks to Heimo TRUHETZ): each MPI task processes one variable, independently
- Model input files are handled via an sbatch run control shell script
- Processing is done on HPC compute node
- Via slurm srun functionality two model domains are processed at the same time on the same node, no overlap between processes
- One processing stream currently uses one compute node: e.g., 2x64core CPUs per node -> 2model domains, 39variables to process each -> processing of two domains, and overall 78 variables at a time
- Individual years (granularity of CORDEX protocol) can be processed one after the other; if your storage permits: submit several processing jobs concurrently, i.e., one year per node, multiple years at the same time
- Apart from some hardcoding issues, see below, there is no need to modify the CMORizer, it is only controlled via a runctrl script and the CMORizer config namelist file
- To balance out wall clock times, variables are split into 4 groups, depending on the time to process them, see below; so to process a complete variable list, the tool has to be run 4 times; if one wants to accespt some inefficiency or if some cores are idle anyway, it would also be possible to process all variables in one go (but they all need to have the same temporal resolution, e.g., hourly sums or instantaneous variables cannot be mixed with daily min max variables)

### Features of the current variant

- Some manual code adjustments are still needed, see below
- Easy to built, use, and adjust, few (standard) software dependencies, with as little components and files to maintain as necessary
- The overall principle is to have a single tool with a small memory footprint, which can take raw model outputs and out them fully standards compliant with as little as possible I/O operations, it should not be necessary to fix any aspect of the resulting netCDF afterwards (compression, attributes, name, etc.)
- Use as part of modelling workflow immediately after a model run, or asynchronous in parallel, or totally independent to process existing model data repositiories
- Any change in the standards can be implmented by code changes, e.g., netCDF variable attribute naming conventions, or by an adjustment or extension of the varibale namelists
- Processing requirements differ very much: Some variables are just passed through (e.g., tas), others need minor calculations (e.g., pr), others require the input of multiple variables with some processing (e.g., zg500), some diagnostics require multiple variableds and extensive processing (e.g., cape)
- Each variable is processed independently, hence the tool can be run for a single variable or a complete set at a time
- A a reference time vector is used; data are read form netCDFF files and sorted in according to time information in the input netCDF file
- Different outut files are possible: yearly, monthly, any timespan, useful for testing and shorter experiments
- An output file can be filled successively and also for different timeslices, there is no check however, whether data for a specific output timestep have already be written and are overwritten, given the input data remain the same this is of no concern.
- The tool is configured through 2 types of namelists: variable namelists, specifying the lists of variables to work on, and experiment and runtime specific settings namelists, which determine, e.g., the root paths for inputs and ouputs or the global attributes of the netCDF file; as the standard variable requirements do not change, these variables barely need any adjustment; for each experiment, e.g., CORDEX-CMIP3 or CORDEX-FPSCONV, one may have a separate runctrl namelist.
- The tool uses relative paths; per tool root directory one can run multiple years concurrently in their respecive year dirctories; one level up, one may have, e.g., 4 different tool diretories to process 4 different variables types concurrently
- In order to account for highly differing wall clock processing times between different types of variables, the CORDEX variables list variables are split into 4 groups; standard near surface variables, sub-daily, daily min/max, 3D pressure level variables, variables/diagnostics which require a lot of computations; with a different model output, this can be adjusted, i.e., all variables could be in a single namelist and the tool could be run over multiple nodes, albeit with a lot of imbalance in between different variable' processing
- Whether input netCDF files contain data of a single day or month or year or range over multiple years does not affect the tool: data are read one timestep at a time, output files are generated as needed; if output is only wanted for a specific timespan (this is an option), then the rest of the input data is ignored, i.e., the time information as required and provided is always matched
- Aside from job submission (with the MPI variant) the tool can be run completely without any run-control shell scripts, solely controlled via its runcontrol namelist and the variable namelist; however there is an interplay between the runcontrol namelist settings, the filename pattern matching in the f90 code, and the data handling by the controlling shell script; e.g., (i) the shell script may stage all input data for a single year, the filename pattern matching, as done in the F90 code using `find` via a system call to generate the filelist to work on, and controlled through the runctrl namelist may be such that it matches all input files with `wrfout` in their filename, then the pCMORizer.f90 would process all input timesteps, that are store din the respective yearly output file; (ii) the controlling shell script may not handle any input raw model output files at all, in this case the path as specified in the runctrl namelist for the input files in combinaiton with the filename pattern matching would genarte the input data filelist
- The flexible handling of input data etc. makes the tool very versatile: not only dioes it feature a small memory footprint, the provisioning/staging/retrieval of inout data can be adjusted to the filesystem environment and especially the storage pace available
- The processing essentially consists of two steps: (i) the f90-tool is used to generate CMOrized variables at the highest temporal resolution for all variables, even if lower temporal resolutions of the tiers in the standard are required; this way all variables from the variable lists are CMORized and raw model outouts might be archived or erased; right now the tool is usually used with 1hr outputs; (ii) the second step consists then of a temporal aggregation of the, e.g., 1hr CMORized data using bash and cdo; originally this aggregation step was also to be included in the f90 tool, but with the MPI-parallelisation  this did not seem feasible and too complicated. If data are not needed at the highest resolution they can be erased; the first CMORization step usually already brings many benefits: lower data volume, standard compliancy, data provenance, more easy data handling.
- The fx variables are also processed seperately using ncgen, ncdump, cdo and bash.
- Right now multiple (2) model domains from a double nest model run are run concurrently on a single compute node; this can be changed to a single domain or instead of domains, output from different component models can be processed; (different filetypes or component models can also be specified in the variable namelists in conjunction with filename pattern matching in the f90 code via the find command system call).
- There is no dependency between individually processed years 

### Other (older) variants

- Older variants are not maintained any more
- Tool was originally designed to be run serially, one variable after another, single process; or concurrently, e.g., by replicating the CMORizer, this also yielded good performance, albeit with far more scripting effort
- The looping constructs and the namelist format in the current tool still resemble the mode of operation ot looping through variable lists, and through temporal aggregations
- Could be run on any machine, local computer as well as HPC compute node, no need for MPI environment
- Processing functionality and basic design has remained unchanged

## Collection of addons

1.

To be integrated into the texts above. 

In case the integrity check shows there is a problem and a single year needs to be repeated, the following preparations are needed before launching the run-ctrl once more:

```
rm d02/
rm $(find ./ -type f -wholename "*_2011*.nc")
```

See the pre-amble of the integrity-checker for details, but if there is a problem with a corrupted (netCDF / HDF5) file, then there is an entry in the err-log file; also, in the out-log file, the `cdo info` command indicates `nan`. Why can this happen? File system issues during raw data writing or later on storage might corrupt a raw model output data file. 

2.

The integrity checker needs to be manually adjusted for the different model domains with respect to the search pattern.
