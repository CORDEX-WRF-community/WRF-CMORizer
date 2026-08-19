# Postprocessing scripts #

The provided postprocessing scripts and codes enable the efficient postprocessing of WRF raw output data, ensuring correct formatting and inclusion of all required variables (e.g. [link](https://docs.google.com/spreadsheets/d/1qUauozwXkq7r1g-L4ALMIkCNINIhhCPx/edit?usp=sharing&ouid=106400114626932444685&rtpof=true&sd=true)) for subsequent analyses.

### Instructions to run the code: ###

1. After cloning the repository, navigate to the pCMORizer folder and:
    - Set the correct compiler (gfortran or intel) and path to your NetCDF-c and NetCDF-Fortran folders in the [Makefile](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/Makefile?ref_type=heads)
   - Adjust the characteristics of your domain and global attributes to fit your domain or create your own template following the [d01](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/runctrl.current.nml_template_d01?ref_type=heads) example
    - Set the correct paths to the WRF raw output files in [run_pCMORizer_in_quasi_parallel.sh](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/run_pCMORizer_in_quasi_parallel.sh?ref_type=heads)
    
2. Run the code:
 - To test the code manually for a single variable:
   ```
  ./run_pCMORizer_in_quasi_parallel.sh <YEAR> <DOMAIN> <VARNAME> <PROJECT> <FREQUENCY>
   ```

   Accepted formats:
    - YEAR: e.g. '2020'
    - DOMAIN: e.g. 'd01' or 'd02'
    - VARNAME: e.g. 'tas'
    - PROJECT: e.g. 'EUROCORDEX', 'I4C', 'STAGE0-URBAN' (any name given to the template - runctrl.current.nml_template_d01_EUROCORDEX)
    - FREQUENCY: 'fx' '1hr' '3hr' '6hr' 'day'

 - To process all requested variables in the selected project, use this command:
   ```
   ./loop_pCMORizer_over_all_variables.sh <VARLIST> <YEAR> <DOMAIN> <PROJECT> <FREQUENCY-OPTIONAL>
   ```
   Accepted formats:
    - VARLIST: a download link for a csv file can be given or csv file itself that has to be located in the running directory
    - YEAR: e.g. '2020'
    - DOMAIN: e.g. 'd01' or 'd02'
    - PROJECT: e.g. 'EUROCORDEX', 'I4C', 'FPS-URBAN' (any name given to the template - runctrl.current.nml_template_d01_EUROCORDEX)
    - FREQUENCY: 'fx' '1hr' '3hr' '6hr' 'day' - optional argument
    
   This script loops over all variables in the given CSV file, reads the required frequency, and if the frequency is "day," it sets the FREQ to "1hr" since the main code does not perform daily aggregation. An additional script is provided for daily aggregation, which reads "1hr" postprocessed data. If the FREQ argument is given, the script subsets the CSV file and extracts only the required frequency for the variables. The script also detects variables that are missing or cannot be postprocessed, providing information on unprocessed variables in a text output file , 'variables_not_processed.txt'.
   
3. Run a script for daily aggreagtion:
   ```
   ./pCMORizer_aggregate_per_day.sh <VARLIST> <PATH_TO_CMORIZED_DATA> <YEAR> <VARNAME(not_obligatory)>  
   ```
   Accepted formats:
    - VARLIST: a download link for a csv file can be given or csv file itself that has to be located in the running directory
    - PATH_TO_CMORIZED_DATA: The full path to the folder with temporal frequencies, e.g. /$PWD/<project_id>/<mip_era>/<activity_id>/<domain_id>/<institution_id>/<driving_source_id>/<driving_experiment_id>/<driving_variant_label>/<source_id>/<version_realization>/ 
     
### Detailed description of the files and their usage:  ###

1. [environment.yaml](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/environment.yaml?ref_type=heads) is a file that you can use to create conda enviroment named "pCMORizer" to run the code with gfortran. To set it run:
  ```
  conda env create -f environment.yaml
  ```
  
2. [Makefile](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/Makefile?ref_type=heads) need to be set up to the corresponding compiler, and is used to compile the main pCMORizer.f90.

3. [pCMORizer.f90](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/pCMORizer.f90?ref_type=heads) is the main Fortran code that reads raw wrfout files, extracts requested variables, and writes and aggregates them in a CF conforming NetCDF file. For the compilation of the code use the [Makefile](https://github.com/impetus4change/wrf/blob/main/CMORizer/Makefile).

4. [CORDEX_CMIP6_variables.csv](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/CORDEX_CMIP6_variables.csv?ref_type=heads) lists all necessary variables and metadata required for CORDEX CMIP6 domains, including the corresponding variable names in the WRF model. The metadata is based on the official CORDEX CMIP6 [variable list](https://cordex.org/experiment-guidelines/cordex-cmip6/data-request-cordex-cmip6-rcms/).

5. [data_request.csv](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/data_request.csv?ref_type=heads) contains the list of requested variables for processing.

6. [generate_vars_namelist.py](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/generate_vars_namelist.py?ref_type=heads) creates a namelist for a desired variable, containing the information necessary for the main code to be executed. To run it, the following command should be used:
   ```
   python generate_vars_namelist.py <variable name(s)>
   ```

7. Template for the core namelist, runctrl.current.nml, for domain [d01](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/runctrl.current.nml_template_d01?ref_type=heads). This template provides information on global attributes of the CMORized file, characteristics of the domain, and the aggregation method to be applied (e.g., from date to date, monthly, yearly).

8. [run_pCMORizer_in_quasi_parallel.sh](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/run_pCMORizer_in_quasi_parallel.sh?ref_type=heads) sends a job to a HPC per variable, per year, per domain, per project, and per frequency.

9. [loop_pCMORizer_over_all_variables.sh](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/loop_pCMORizer_over_all_variables.sh?ref_type=heads) is a small wrapper script that loops over all variables in the chosen CSV table. 

10. [pCMORizer_aggregate_per_day.sh](https://github.com/CORDEX-WRF-community/WRF-CMORizer/blob/main/pCMORizer_aggregate_per_day.sh?ref_type=heads) is a script that loops over all variables from the chosen CSV file with daily frequencies and uses cdo and nco tools to obtain daily mean, maximum, and minimum values. 
