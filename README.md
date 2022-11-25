# pCMORizer.f90

2022-11-25

The pCMORizer.f90 is an open-source Fortran-based software tool with some ancilliary bash scripts that also use [cdo](https://code.mpimet.mpg.de/projects/cdo/). pCMORizer.f90 and companion tools are used transfer or postprocess raw (regional climate) model ouputs into netCDF files, which comply to the [CMOR](https://cmor.llnl.gov) standard according to the [CORDEX-CMIP3 archive specification](http://is-enes-data.github.io/cordex_archive_specifications.pdf) as part of the [CORDEX WCRP project](https://cordex.org/experiment-guidelines/how-to-submit-data-rcms/). The goal is to make the model data compliant to be disseminated through [ESGF](https://esgf-data.dkrz.de/projects/esgf-dkrz/), and/or generate a dataset which is more efficient to be used in analysis, and local sharing, etc.

## Motivation, features, design principles

At an early stage of the CORDEX-CMIP3 experiment in 2013, there were very few tools available, which would ingest raw RCM output and generate CMORized simulation results for dissemination and/or efficient in-house usage or archival, meeting FAIR principles. Existing solutions at the time often consisted of complex, difficult-to-maintain script-based processing chains, that resulted in lots of time-consuming I/O operations, issues when it came to the generation of more complex variables and diagnostics (as script-based tools could not be used for this purpose any more), or large memory footprints. 

In this context the pCMORizer.f90 was orginally developed as a purely serial tool for WRF RCM outputs. The orignal version was expanded by the community now and then, partly adopted to specific workflows, and further customized. Full operational capability, was only reached step by step, triggered by archival, dissemination, and publication needs. Today, a range of CMORizer tools exist for various RCMs, see the [Links section](#ref_links). The pCMORizer.f90 has undergone three major development stages (or variants), only the most current is maintained and documented here. 

## Maintained MPI-parallel variant

- Two-staged processing: Fortran tool does initial CMORization step, bash+cdo does temporal aggregation
- MPI-parallel (thanks to Heimo TRUHETZ): each MPI task processes one variable, independently
- Model input files are handled via an sbatch run control shell script
- Processing is done on HPC compute node
- Via slurm srun functionality two model domains are processed at the same time on the same node, no overlap between processes
- One processing stream currently uses one compute node: e.g., 2x64core CPUs per node -> 2model domains, 39variables to process each -> processing of two domains, and overall 78 variables at a time
- Individual years (granularity of CORDEX protocol) can be processed one after the other; if your storage permits: submit several processing jobs concurrently, i.e., one year per node, multiple years at the same time
- Apart from some hardcoding issues, see below, there is no need to modify the CMORizer, it is only controlled via a runctrl script and the CMORizer config namelist file
- To balance out wall clock times, variables are split into 4 groups, depending on the time to process them, see below; so to process a complete variable list, the tool has to be run 4 times; if one wants to accespt some inefficiency or if some cores are idle anyway, it would also be possible to process all variables in one go (but they all need to have the same temporal resolution, e.g., hourly sums or instantaneous variables cannot be mixed with daily min max variables)

## Features of the current variant

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

## Compatibility and standards compliance

pCMORizer.f90 was originally developed to cmorize raw [WRF RCM](https://www2.mmm.ucar.edu/wrf/users/) outputs.

It seems generic enough in terms of structure and configuration to also be used for other RCMs, LSMs, HMs. 

At this point, the CORDEX-CMIP3 and the FPSCONV experiments archive specifications and variable lists (VLs), in their (most current) versions as of 2022-09-25, are supported, see the [references](#ref_refs) and linked documents therein.

Some (for the majority of the users not so relevant) variables cannot be derived (yet), see [ToDos](#ref_todos).

## Usage<a name="ref_usage"></a>

### Computational environment / requirements

Needed for the currently maintained variant of the tool and the scripts provided (GCC compile is possible, but not provided at this point).

- Intel Fortran compiler
- netCDF library
- ncdump
- ncgen
- bash
- make
- find
- uuidgen
- sbatch/srun HPC schedulaing environment
- cdo

### Configuration, compute environment

For bash shell, set the pCMORizer root directory of the tool, this variable is just used in this manual to make things clear:

```
export pCMORizer_DIR=~/pCMORizer
```

Load the compile and runtime environment, loaded by the runcontrol script during runtime:

```
cd $pCMORizer_DIR
source loadenv.JURECA-DC_2020_Intel-PSMPI.ini
```

### Manual code adjustments and compilation

More information / fixes to come. These things have to be manually adjusted in `pCMORizer.f90`.

- T00 (290K or 300K, check your Registry.COMMON)
- depths of LSM
- max number of vars
- which type of pressure calc
- handling of the time info in the files, 00-23, or 00-00
- pattern matching of the filenames if they are not WRF default, filename pattern with find
- hourly vs daily (for tasmin/tasmax)

If the code has been adjusted to (i) your filestructure, (ii) model config during simulation, (iii) filename patterns, (iv) storage pattern, the only thing that needs adjustment is the 1hr vs daily setup.

```
vim pCMORizer.f90
make veryclean
make
```

### Configuration, CMORizer operation

The main `runctrl.current.nml` namelist configures the global attributes and defines all directories and the year range. The current configuration uses '0000' as start and end year, this means the file name pattern used used by the `find` command in `pCMORizer.f90` does not apply any years as a pattern for the input files.

Here two `runctrl.current.nml` files are provided one for each model domain; the runctrl files are documented and determine all major metadata and domain settings and also point to static files which need to be provided for the processing (`geo_em` files).

If the experiment does not change, no changes are needed here. If a new experiment is to be processd, adjust these files accordingly.

```
vim -R -o runctrl.current.nml_template_d01_DA runctrl.current.nml_template_d02_DA
```

The variable namelists do not need to be changed usually. There are 4 variable lists, one per "variable group". If a variable is not needed, this variable has to be erased form the list; the functionality which could be used before to switch a variable ON and OFF is essentially deprecated with the MPI version of the tool (although the switch as such still works).

If the protocol does not change, if no new variables are required, or some variables are not required anymore, no changes are needed here.

```
vim runctrl.vars.std_sfc.nml
vim runctrl.vars.std_minmax.nml
vim runctrl.vars.std_presslev.nml
vim runctrl.vars.special.nml
```

The run control script (i) handles the data and (ii) starts the CMORizer, making sure via srun options, that there is no overlap as `srun` submitted two times on the same compute node. Here multiple manual adjustments are needed. See comments in the code (sbatch options, wallclock times, see performance section; linking the correct variable and runctrl namelists).

```
vim pCMORizer_runctrl.sh
```

### Running the CMORizer

There are many ways of running the CMORizer tool, depending on processing setup, urgency, amount of data, available resources. 

The currently maintained variant uses two major steps:
1. CMORize all data at highest temporal resolution, this is what the pCMORizer.f90 does
2. Temporally aggregate data as required by the protocol using cdo and bash (this was planned to also go into f90) but with the MPI restructuring it seemed quicker to implement it like this

**Step 1: CMORization at highest temporal output interval, all variables**

Two options are curretnly supported.

1. Run the tool for a complete variable set at a time, but just run a single year (or timespan), the compute environment is loaded automatically, within the current working directory a `<year>/` directory is created; all scripts and namelists will be contained in the year-directory. The input data is symbolically linked into this directory; after succesful processing, this directory can be erased, manually. The postprocessed / CMORized data are stored as specified (in their own separate directory). (The input directory could also be anywhere.) All needed input files (runcontrol namelist, variable namlist, executable, input files are in the run-dir, i.e., the year directory.

```
cd $pCMORizer_DIR
export Y=1998 && sbatch --export=ALL,Y=$Y --job-name=pCMORizer$Y pCMORizer_runctrl.sh
```

More file-handling could be added to this script. In our case, all data is on spinning disk and readily available to the tool.

After successful postprocessing:

```
cd $pCMORizer_DIR
rm -rf 1998
```

2. Just run the above start sequence x10 to process 10 years on 10 nodes at the same time, for the same variable type. 10 year-directories are created, the tool runs completely independent within these directories.

```
cd $pCMORizer_DIR
./launcher_pCMORizerMultipleYears_sameVars_1yPerCompNode.sh
```

**Step 2: Temporal aggregation -- under development and implementation, 2022-11-16**

...

## Performance<a name="ref_perf"></a>

On [JURECA-DC compute nodes](https://apps.fz-juelich.de/jsc/hps/jureca/index.html) in combination with [JUST](https://www.fz-juelich.de/en/ias/jsc/systems/storage-systems/just) [filesystems](https://apps.fz-juelich.de/jsc/hps/jureca/filesystems.html) the following performance can be reached. 

Use the information below to set your sbatch wallclock times.

**Step 1: CMORization at highest temporal output interval, all variables**

Setup: Input data is WRF raw model outputs (`wrfout`), compressed; one netCDF file per day; one compute node processes one year at a time; two domains are processed concurrently per node (EUR-15, ALP-3); each CPU-core processes one variable; different variable types are processed in seperate, manually started and adjusted processing stepsi; ten compute nodes are run via 10 sbatch commands concurretlyi; Intel compile, with `-O3`; logging is OFF. See [Usage section](#ref_usage) above for details; given are the approximate average runtimes accross the 10 compute jobs; all model outputs are provided on spinning disk at the same time.

- 39 near-surface variables: about 5h25min (single month test: 30min suffices)
- 36 pressure level variables: about 6h30min
- 2 min/max variables: about 5min
- 2 variables with more extensive diagnostics calculations:

I.e., using 10 compute nodes with exclusive allocation, a decade long km-scale experiment can be CMORized in less than 24h wall clock time.

## ToDos<a name="ref_todos"></a>

The current version is in the VERSION.txt file. Though the tool is operational some ToDos are pending towards v1.0.

Most urgent for 2022 ESGF upload:

- bash and cdo based temporal aggregation
- bash and cdo based fx variables processing
- Adjustments in conjunction with QA check

Further modfications:

- Adjust for additional RCMs, LSMs, HMs based on the exsting structure and configuration namelists

Nice to have:

- Addition of still missing variables (such as cll, `clh`, `clm` incl. `plev_bnds`), and wind speed extremes
- More efficient calculation of more complex diagnostics, such as `cape` and `cin`
- Compilation using GCC and OpenMPI, fix of some compile time warnings, which lead to runtime abort with GCC
- Improved chunking, although a series of tests has shown no improvement when netCDF file chunking is altered to what the netCDF library does as a default
- All temporal aggregations and regriddings as required by the complete CORDEX-CMIP3 standard

Future Todos:

- OpenMP in combination with MPI to speed up calculation further
- Extension to be used with the CORDEX-CMIP6 experiment

## Documentation

This README.md file and comments in the files are only documentation.

## Release

Written and tested by in the CONTRIBUTORS.txt file.

pCMORizer.f90 is released under the MIT License.

For details and restrictions, please read the LICENSE.txt file.

## Some notes on support for TSMP

Generalisation of the features of the tool; instead of two domains, process component models concurrently; some specific adjustments are needed.

Below are some development notes, not the final documentation:
- Start test dev with COSMO EUR-0275 output, most similar to WRF
- The different component models need a somwhat adjusted treatment; WRF is using absolute and not relative times, use of https://code.mpimet.mpg.de/projects/cdi for Fortran to extract time information; various tests needed; eventually it seems the tool can only be used for classic netCDF files, hence a different solution, perhaps a Fortran subroutine or a system call to `date` is needed
- Code scripts and tests to include the new CDI library for the moment
- The Intel/PSMPI 2020 environment seemed to have issues with CDI; hence use 2022 GCC and OpenMPI, this is on the task list anyway; get rid of some compiler warnings from GCC, also beneficial for Intel; also the GCC compiler produces much smaller executable
- Using the GCC/OpenMPI environment for pCMORizer.f90 and CDI, there are issues with the parallel execution, i.e. as soon as going to more than one variable there is an issue, something related with the MPI-SCATTER perhaps, maybe interplay with the slurm parameters
- Right now only variables are ingested which do not need additional processing
- No adjustment of the attributes and meta data so far (nml)
- (Other tools from the netCDF library or, e.g., UDUNITS, also do not solve the time information conversion issue easily.)
- The serial netCDF libs do not work properly with the mpifort compiled code -> netCDF is not working in this case, albeit without error message

## References<a name="ref_refs"></a>

- [CF-conventions](https://cfconventions.org/)
- [CORDEX-CMIP3 downscaling experiment guideleines documents](https://cordex.org/experiment-guidelines/how-to-submit-data-rcms/)
- [Specific reuqirements of the CORDEX FPSCONV initiative](https://www.hymex.org/cordexfps-convection/wiki/doku.php?id=protocol)
- [DKRZ QA checker](https://github.com/IS-ENES-Data/QA-DKRZ)
- [ESGF infrastructure data node for public dissemination](https://esgf-data.dkrz.de/projects/esgf-dkrz/)

## Links<a name="ref_links"></a>

Similar tools and resources:

- [https://github.com/C2SM-RCM/CCLM2CMOR](https://github.com/C2SM-RCM/CCLM2CMOR)
