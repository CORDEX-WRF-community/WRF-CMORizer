!===============================================================================
! BOP
!
! NAME:
!   postpro_model_WRF_to_ESGcompliancy.f90
!   See license information at the end of the preamble.
!
! VERSION:
!   v2013-08-09
!   see git log for revision details and history
!
! STATUS:
!   under development
!
! CURRENT / (FORMER) CODE OWNER(S):
!   Klaus GOERGEN | k.goergen@gmx.net | KGo | MIUB/JSC
!   Sebastian KNIST | sebastian.knist@gmx.de | SKn | MIUB/JSC
!
! PURPOSE / DESCRIPTION:
!   This application postprocesses (standard) raw WRF simulation results into
!   CORDEX (and CMIP5) compliant NetCDF files in a dedicated directory tree.
!   See the references below for the specifications used for this program.
!   The code can easily be adjusted in case specifications change. The post-
!   processing is necessary in order to being able to upload, stage and
!   distribute RCM simulation results via the Earth System Grid infrastructure
!   adopted by CORDEX from CMIP5. Files which do not adhere to this are
!   rejected.
!   Other, similar tool-sets are e.g.: (i) CMOR, (ii) XXXX wrf spain [...]
!   This tool has be seen in conjunction with 
!   (1) [...]
!   (2)
!
! CONVENTION:
!   crpgl_ucc_v01
!
! PRG.-LANGUAGE / ENVIRONMENT:
!   - >=F95 compiler
!   - used: gfortran 4.6.3, ifort 11.1
!   - Linux/UNIX OS (-> system calls)
!   - F95 ISO FORTRAN except "SYSTEM" intrinsic function
!
! REQUIREMENTS:
!   - FORTRAN95 compiler
!   - NetCDF F90 library (http://www.unidata.ucar.edu/software/netcdf/)
!     v4.x, used: v4.1.1 and 4.2.1.1 incl. HDF5 -> write NetCDF-4 classic model
!     format
!   - make *nix console application
!   - date *nix console application
!   - uuidgen *nix console application (http://www.uuidgen.com/)
!
! BUILDING:
!   a) command line:
!      gfortran -I/usr/include <prg>.f90 -L/usr/lib -lnetcdff -lnetcdf
!   b) Makefile:
!      make
!
! CATEGORY:
!   PostPro.RCM.WRF.
!
! CALLING SEQUENCE:
!   ./postpro_model_WRF_to_ESGcompliancy > log
!
! CALLED FROM:
!   Standalone command-line tool.
!
! LOCAL VARIABLES:
!   See the variable declarations at the beginning of the code.
!
! INPUTS:
!   - NML files, see the tools main directory, have to be adjusted.
!   - WRF static fields (geo_em*)
!   - WRF outputs (wrfout*)
!   - WRF outputs extremes (wrfxtrm*)
!   No optional inputs, no keyword parameters.
!
! FILES USED:
!   - See inputs.
!
! OUTPUTS:
!   - NetCDF files according to standard specification.
!   - No Optional outputs.
!
! RESTRICTIONS:
!   - No restrictions. No side effects.
!
! BUGS:
!   - None.
!
! PROCEDURE / FEATURES:
!   - The tool can produce all required variables, i.e. output NetCDF files as
!     defined in the CORDEX archive design specifications as available in
!     July/August 2013 in possibly one pass based on standard WRF simulation
!     outputs. No additional processing is needed. Also the reuiqred NetCDF-4
!     classic data model format is written immediately, no later conversion
!     needed.
!   - 3hr data is produced first. This is closest to outputs most groups have
!     anyway. It is possible to produce more 3hr output than needed, i.e. for
!     all variables specified plus more vertical levels. Thereby the tool might
!     also be used for a general data volume reduction after a model run.
!   - The application loops over various averaging intervals, then over
!     variables and then over the existing WRF data files.
!   - Very large WRF outputs and/or even larger model domains than used in
!     CORDEX can be handled as the tool is working on individual output
!     intervals only. The tradeoff is a more I/O overhead but lower RAM
!     requirements.
!   - The FORTRAN code is ISO FORTRAN95-compliant except for the "SYSTEM" calls.
!     Hence the tool should be portable easily.
!   - The tool has minimum requirements in terms of libraries or external tools.
!   - FORTRAN was chosen for its fast execution speed, suitable for the
!     large datasets and the possibility to use the tool in HPC environments as
!     well as on individual workstations with a good performance. Also it is
!     easily possible to parallelize portions of the code via e.g. OpenMP.
!   - By splitting the *vars* namelist, the tool may be run concurrently for
!     different portions of a WRF output dataset.
!   - No code modifications are needed to run the tool and customize it for use
!     with a different WRF experiment. Only the NML files need to be changed.
!   - Different WRF input sources are possible, currently the tool is designed
!     for wrfout and wrfxtrm files.
!   - The tool creates a reference time vector. The time-information contained
!     in each original WRF sim. file is retrieved and according to this
!     information data is sorted into the resulting NetCDF files. This makes the
!     the tool rather robust and flexible, a tradeoff is the possibly longer
!     processing time due to the searches needed for the date and time matching.
!     However these searches are on subsets only and therefore fairly efficient.
!     It also means that the WRF outputs may cover any timespan, daily, monthly,
!     or any overlap of months and/or years and that they may come in filelists
!     even not temporally ordered.
!   - Currently only one input root directory is possible. If data is stored at
!     different locations symbolic links might have to be done beforehand.
!   - The static fields are treated independently by the tool.
!   - Currently the namelists are split into several parts but they may also 
!     be combined. 
!   - The tool is also intended to reduce WRF model output data volume. This
!     means that original raw model outputs are likely to be erased afterwards.
!     Therefore the tool generates slightly more variables than required by the
!     CORDEX data protocol: e.g. CAPE, ....
!
! EXAMPLE:
!   ./postpro_model_WRF_to_ESGcompliancy > log
!
! MODIFICATION / REVISION HISTORY:
!   See either git log or NEWS for details.
!   2013-08-09_KGo v0.1
!
! TODO / PLANNED EXTENSIONS:      ------    outdated
!   - Temporal aggregations, i.e. 6hr, day, mon, seas; all based on orginal
!     outputs
!   - Static fields processing
!   - All variables -> extension of namelist
!   - Variable-dependent processings, e.g. MSLP
!   - time_bnds and lvl vars and processing, ask Grigori about this 00,24 thing
!   - Import of alternative (to NML) ASCII file with long and standard names
!   - Add units to that ASCII file as well
!   - Not yet tested with EUR-11, i.e. large model outputs
!   - Additional variable in runctrl.vars.nml to control additional 3hr outputs
!     to have all vars in that format and have more vertical levels
!   - (OpenMP parallelism for the processing section), via pre-processor flags
!   - (Parallel NetCDF I/O where possible), via pre-processor flags
!
! REFERENCES (some reference tool format):
!   - CORDEX WRF group model identification and naming:
!     https://docs.google.com/spreadsheet/ccc?key=0ArYFyU35McvvdFBqaXdLcERjbFp3U
!     lBZcC1qbm53NFE#gid=0
!   - Standard specification / naming conventions (see NML files)
!     ... CMIP5, CORDEX, txt files
!
! CALLED PROCEDURES:
!   No external calls. -- System calls are needed.
!
! PERFORMANCE:
!   EUR-44: 1min/yr > 3hr/150yr OR 1h/1yr65vars... + averaging, after each run
!
! LICENSE / COPYING:    ------ replace by MIT license, see github
!
!   Copyright (C) 2013 Klaus GOERGEN
!
!   This file is part of postpro_model_WRF_to_ESGcompliancy.
!
!   postpro_model_WRF_to_ESGcompliancy is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 3 of the License, or any later
!   version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program. If not, see <http://www.gnu.org/licenses/>.
!
! EOP
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! passing allocatable arrays between main program and external subroutine

MODULE flhandling
	
	IMPLICIT NONE
	SAVE
	
	CHARACTER (len = 200):: tmpfileFL
	CHARACTER (len = 200), DIMENSION(:), ALLOCATABLE :: fl_wrfout
	CHARACTER (len = 200), DIMENSION(:), ALLOCATABLE :: fl_wrfxtr
	INTEGER :: ft
	
END MODULE flHandling

!-------------------------------------------------------------------------------
! index, Y, M, D, H, depending on the "frequency" of the dataset, i.e. time
! intervals

MODULE RefTimeVecs
	
	IMPLICIT NONE
	SAVE
	
	REAL, DIMENSION(:,:), ALLOCATABLE :: TimeRefArray
	!INTEGER :: Tyear_start, Tyear_end
	
END MODULE RefTimeVecs

!-------------------------------------------------------------------------------
! namelist handling

MODULE NameListHandling
	
	IMPLICIT NONE
	SAVE
	
	INTEGER, PARAMETER :: nvars = 13 
	
	CHARACTER (len = 200) :: Conventions, contact, experiment_id, experiment, &
	driving_experiment, driving_model_id, driving_model_ensemble_member, &
	driving_experiment_name, institution, institute_id, model_id, &
	rcm_version_id, project_id, CORDEX_domain, product, references, forcing, &
	source, title
	
	CHARACTER (len = 200) :: comment, institute_run_id
	
	CHARACTER (LEN = 100), DIMENSION(nvars):: var_wrf, var_cmip, standard_name, &
	long_name, units, filetype, cm1hr, cm3hr, cm6hr, cmDay, cmMon, cmSea, positive
	INTEGER, DIMENSION(nvars):: height, cordexID
	LOGICAL, DIMENSION(nvars):: time1hr, time3hr, time6hr, timeDay, timeMon, timeSea, &
	interpolate
	
	CHARACTER (len = 200) :: DirInputSimResRoot, DirOutputPostProRoot, domain
	
	INTEGER ::  nx, ny, nz, xoffset, yoffset, xfocus, yfocus
	CHARACTER (len = 4) :: ts, te, tstot, tetot
	
	CHARACTER (len = 300) :: PnFnGeo
	
	NAMELIST / globalvars / Conventions, contact, experiment, experiment_id, &
	driving_experiment, driving_model_id, driving_model_ensemble_member, &
	driving_experiment_name, institution, institute_id, model_id, &
	rcm_version_id, project_id, CORDEX_domain, product, references, &
	forcing, source, title
	
	NAMELIST / globalvars_additional / comment, institute_run_id
	
	NAMELIST / vars / var_wrf, var_cmip, standard_name, long_name, units, &
	height, time1hr, time3hr, time6hr, timeDay, timeMon, timeSea, filetype, &
	cm1hr, cm3hr, cm6hr, cmDay, cmMon, cmSea, interpolate, cordexID, positive
	
	NAMELIST / filesystem / DirInputSimResRoot, DirOutputPostProRoot, domain
	
	NAMELIST / model_config / ts, te, nx, ny, nz, xoffset, yoffset, xfocus, &
	yfocus, tstot, tetot
	
	NAMELIST / static_fields / PnFnGeo
	
END MODULE NameListHandling

!===============================================================================

PROGRAM ppWRFCMIP
	
	USE flhandling
	USE RefTimeVecs
	USE NameListHandling
	
	USE netcdf
	
	IMPLICIT NONE
	
	!===============================================================================
	
	INTERFACE
		
		SUBROUTINE generateFilelist
		END SUBROUTINE generateFilelist
		
		SUBROUTINE CreateRefTimeArray( dt, dataIncludesLeapYearDays )
			IMPLICIT NONE
			CHARACTER (LEN = 3), INTENT(IN) :: dt
			INTEGER :: dataIncludesLeapYearDays
		END SUBROUTINE CreateRefTimeArray
		
	END INTERFACE
	
	!===============================================================================
	! filenames
	
	CHARACTER (len = *), PARAMETER :: fnNMLexp = "runctrl.mpiesmlrhist.nml" !"runctrl.erainteval.nml" !"runctrl.access13hist.nml"
	
	!CHARACTER (len = *), PARAMETER :: fnNMLvar = "runctrl.vars.nml" !"runctrl.vars.nml_evp_roff" !"runctrl.vars.nml_water_column" ! "runctrl.vars.nml_vars_on_plevels"  !"runctrl.vars.nml_vars_on_plevels" !"runctrl.vars.nml_pr"
	CHARACTER (len = 100), DIMENSION(:), ALLOCATABLE :: fnNMLvar
	
	!CHARACTER (len = *), PARAMETER :: PathFileNameInTEST = "testWRFin.nc"
	!CHARACTER (len = *), PARAMETER :: PathFileNameOutTEST = "testESGout.nc"
	
	CHARACTER (len = 200) :: pn_out, fn_out, iflWRFin
	
	!-------------------------------------------------------------------------------
	
	! auxilliary vars, just needed during development
	! INTEGER, PARAMETER :: nt = 8
	
	! new NetCDF file
	INTEGER :: ncid, ncidin, ncidin0
	INTEGER :: lon_dimid, lat_dimid, rec_dimid, height_dimid, &
	nb2_dimid
	INTEGER :: varid, x_varid, lon_varid, lat_varid, rlon_varid, rlat_varid, &
	rotated_pole_varid, height_varid, rec_varid, pp_varid, pb_varid, ph_varid, &
	phb_varid, qv_varid, qc_varid, qi_varid, qr_varid, qs_varid, &
	theta_varid, t2_varid, recbnds_varid, rainnc_varid, &
	rainc_varid, snownc_varid, u10_varid, v10_varid, u_varid, v_varid, &
	sfcevp_varid, potevp_varid, sfroff_varid, udroff_varid, acsnom_varid, &
	sinalpha_varid, cosalpha_varid
	
	! input data general query
	INTEGER :: ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in !!!, formatp_in
	INTEGER :: nvar_nml
	! record variable in input data
	INTEGER :: InVarIdRec, InDimLenRec !!!, InVarNdimsRec
	CHARACTER (len = NF90_MAX_NAME) :: InDimNameRec !!!, InVarNameRec
	!!!INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: InDimIdsRec
	CHARACTER (len = 19), DIMENSION(:), ALLOCATABLE :: InVarDataRec
	
	! data
	INTEGER, DIMENSION(:), ALLOCATABLE :: counter_array
	REAL, DIMENSION(:,:), ALLOCATABLE :: temp_data, temp_data_Month_Average, &
	temp_data_Season_Average !For averaging reasons...
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: temp_data_InTime,temp_data_Monthly
	REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: temp_data_Seasonally!To use as temporal storage..
	
	
	REAL, DIMENSION(:,:), ALLOCATABLE :: data_in, psl_in, t2_in, TimeRefArraySelYear, &
	cldfra_inv, u10_in, v10_in, cape, cin, lcl, lfc, prw, clwvi, clivi, &
	sinalpha_in, cosalpha_in, Time_bnds 
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: pp_in, pb_in, ph_in, phb_in, qv_in, qvs, &
	qc_in, qr_in, qi_in, qs_in,  theta_in, t_in, ph_fl, p_in, cldfra_in, &
	u_in, v_in, var3d_in, var_pl, potevp_in, &
	rainnc_in, rainc_in, rad_in, t_p, snownc_in, acsnom_in, GeoInLonLat, &
	sfcevp_in, sfroff_in, udroff_in
	REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: smois_in
	DOUBLE PRECISION :: fct=100.0, trsl=180.0
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GeoInRLat, GeoInRLon
	REAL, DIMENSION(:), ALLOCATABLE ::  pout
	INTEGER (KIND=4):: ix, jy
	REAL :: t_ii, dtHours
	
	REAL, PARAMETER :: cp = 1004 !J kg-1 K-1
	REAL, PARAMETER :: R = 287.05 !J kg-1 K-1
	REAL, PARAMETER :: L = 2501000. !J kg-1
	REAL, PARAMETER :: a = 610.78 ! Pa
	REAL, PARAMETER :: b = 17.27 !
	REAL, PARAMETER :: c = 273.15 !
	REAL, PARAMETER :: d = 35.86 !
	REAL, PARAMETER :: n = L*0.622*a/cp !
	
	
	! time and date handling
	CHARACTER (len = 3), DIMENSION(:), ALLOCATABLE :: frequency
	!!!INTEGER :: WRFfileNyears, WRFfileNmonths
	INTEGER, DIMENSION(:), ALLOCATABLE :: InDateTimeYear, InDateTimeMonth, &
	InDateTimeDay, InDateTimeHour      !, WRFfileIyears, WRFfileImonths
	CHARACTER (LEN=4) :: InDateTimeYearStr, EndInDateTimeYearStr
	!CHARACTER (LEN=2) :: InDateTimeMonthStr
	INTEGER :: InDateTimeYearPrev = 0  !, InDateTimeMonthPrev = 0
	
	CHARACTER (LEN = 100), DIMENSION(nvars) :: cell_methods
	
	!-------------------------------------------------------------------------------
	! statistics
	
	REAL :: stat_mean, slope
	
	!-------------------------------------------------------------------------------
	! general
	!!!!!!!!!!!!!!! total_iters was added 13/6/2017 -- ARIS EDIT !!!!!
	INTEGER :: i, sts, ivar, ifrq, ifl, it, counter, j, np, nl, ii, varnml, &
	obs_interval, record_count, writeIndex, write_iter_count, &
	monthIndex, HoursOfCurrentInterval, isEndOfMonth, isEndOfSeason, &
	currentDayOfMonth, currentMonthOfYear, dayCount, counter_prev_year = 0, &
	FirstYearOfRun, LastYearOfRun, dataIncludesLeapYearDays! !!!j
	! functions
	INTEGER :: CheckForLeapyear
	! other general...
	INTEGER, DIMENSION(1:12) :: NoOfDaysPerMonth ! (31, 29, 31, 30 ... month days!)
	!!!INTEGER :: AllocateStatus, DeAllocateStatus
	LOGICAL :: FileExists   !, comb_flags
	REAL :: cpuTs, cpuTe
	! Command-line arguments for executable...
	INTEGER :: NumberOfArguments, ArgumentsCounter !firstvar_nvar_nml, first_varnml ! number of arguments given, and counter for later loop...
	CHARACTER (len = 20) :: ArgumentName ! name of argument...
	CHARACTER (len = 4) :: dirLen	! length of path name up until last slash ('/') before actual wrfout input data. Used in filesearch pattern
	
	
	
	!-------------------------------------------------------------------------------
	! system calls
	
	CHARACTER (len = *), PARAMETER :: cmdUUID = "uuidgen -t > tmpfileUUID"
	CHARACTER (len = 37) :: trackingID
	
	CHARACTER (len = *), PARAMETER :: cmdDate = "date -u +%Y-%m-%dT%H:%M:%SZ > tmpfileDate"
	CHARACTER (len = 21) :: creationDate, EndDate
	
	
	
	
	!-----------------------------------------
	!Handling command line arguments.
	!Currently supported: frequency, range of nml files, range of variables in selected nml files
	
	!Check if any arguments were given
	NumberOfArguments = command_argument_count()
	
	!Loop over arguments and options
	IF (NumberOfArguments > 0) THEN
	
		DO ArgumentsCounter = 1, NumberOfArguments, 1
		
		CALL get_command_argument(ArgumentsCounter, ArgumentName)
		
			SELECT CASE (adjustl(ArgumentName))
			
				CASE ("1")
				
					ifrq = 1;	!3hr
					
				CASE ("2")
				
					ifrq = 3;	!daily
					
				CASE ("3")
				
					ifrq = 4;	!monthly
					
				CASE ("4")
				
					ifrq = 5;	!seasonal
					
					
					
			END SELECT
			
		END DO
		
	END IF
	
	!===============================================================================
	
	PRINT *, "============================================================"
	PRINT *, "*** NML READING ***"
	PRINT *, fnNMLexp
	!PRINT *, fnNMLvar
	
	OPEN(2,FILE=fnNMLexp)
	READ(UNIT=2,NML=globalvars)
	CLOSE(2)
	
	OPEN(2,FILE=fnNMLexp)
	READ(UNIT=2,NML=globalvars_additional)
	CLOSE(2)
	
	OPEN(2,FILE=fnNMLexp)
	READ(UNIT=2,NML=filesystem)
	CLOSE(2)
	
	OPEN(2,FILE=fnNMLexp)
	READ(UNIT=2,NML=model_config)
	CLOSE(2)
	
	OPEN(2,FILE=fnNMLexp)
	READ(UNIT=2,NML=static_fields)
	CLOSE(2)
	
	!OPEN(2,FILE=fnNMLvar)  !SKn loop over sereval var namelists later
	!READ(UNIT=2,NML=vars)
	!CLOSE(2)
	
	!-------------------------------------------------------------------------------
	! allocate main data input array outside the loops based on nml entries
	! dummy allocation of the ref time array -> has to maintain its values during
	! several looping constructs and is de-allocated before the initial allocation
	
	ALLOCATE( data_in( xfocus, yfocus ), STAT=sts )
	IF (sts /= 0) STOP "*** Not enough memory ***"
	
	ALLOCATE( temp_data(xfocus, yfocus), STAT=sts )
	IF (sts /= 0) STOP "*** Not enough memory ***"
	
	ALLOCATE( temp_data_Month_Average(xfocus, yfocus), STAT = sts)
	IF (sts /= 0) STOP "*** Not enough memory ***"
	
	ALLOCATE( temp_data_Season_Average(xfocus, yfocus), STAT = sts )
	IF (sts /= 0) STOP "*** Not enough memory ***"
	
	ALLOCATE( TimeRefArraySelYear(2,2) )
	ALLOCATE( Time_bnds(2,2) )  !SKn
	!-------------------------------------------------------------------------------
	! get the invariant vars which have to be added all the time
	! lon, lat, rlon, rlat
	! mass grid
	! seperate file
	! subset is double checked with orig geo_em file and previously postprocessed
	! data; match up to the 5th digit
	! geo files match
	
	PRINT *, "============================================================"
	PRINT *, "*** STATIC FIELDS ***"
	PRINT *, TRIM(PnFnGeo)
	
	ALLOCATE( GeoInLonLat(yfocus, xfocus, 2) ) ! F95 order
	PRINT *, SHAPE(GeoInLonLat)
	ALLOCATE( GeoInRLat(yfocus) )
	ALLOCATE( GeoInRLon(xfocus) )
	
	sts = NF90_OPEN(TRIM(PnFnGeo), NF90_NOWRITE, ncidin)
	
	sts = NF90_INQ_VARID(ncidin, "XLONG_M", varid)
	!sts = NF90_INQ_VARID(ncidin, "XLONG", varid)
	sts = NF90_GET_VAR(ncidin, varid, GeoInLonLat(:, :, 1), &
	START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /)) !normal order: x y z t
	
	sts = NF90_INQ_VARID(ncidin, "XLAT_M", varid)
	!sts = NF90_INQ_VARID(ncidin, "XLAT", varid)
	sts = NF90_GET_VAR(ncidin, varid, GeoInLonLat(:, :, 2), &
	START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /))
	
	sts = NF90_INQ_VARID(ncidin, "CLONG", varid)
	!sts = NF90_INQ_VARID(ncidin, "XLONG", varid)
	sts = NF90_GET_VAR(ncidin, varid, GeoInRLon(:), &
	START = (/ xoffset, 1, 1 /), COUNT = (/ xfocus, 1, 1 /))
	
	sts = NF90_INQ_VARID(ncidin, "CLAT", varid)
	!sts = NF90_INQ_VARID(ncidin, "XLAT", varid)
	sts = NF90_GET_VAR(ncidin, varid, GeoInRLat(:), &
	START = (/ 1, yoffset, 1 /), COUNT = (/ 1, yfocus, 1 /))
	
	sts = NF90_CLOSE(ncidin)
	
	! Change rlon (add 180 degrees)
	DO i=1,size(GeoInRLon(:))
		
		IF ( GeoInRLon(i) .lt. 0 ) THEN
			
			ix=ANINT((GeoInRLon(i)+trsl)*fct)
			GeoInRLon(i)=ix/fct
		
		ELSE
			
			ix=ANINT((GeoInRLon(i)-trsl)*fct)
			GeoInRLon(i)=ix/fct
			
		END IF
		
	END DO
	
	DO j=1,size(GeoInRLat(:))
		
		jy=ANINT((GeoInRLat(j)*fct))
		GeoInRLat(j)= jy/fct
		
	END DO
	
	!PRINT *, "rlon = "
	!PRINT *, SHAPE(GeoInRLon)
	!PRINT *, GeoInRLon
	
	!PRINT *, "rlat = "
	!PRINT *, SHAPE(GeoInRLat)
	!PRINT *, GeoInRLat
	
	!PRINT *, "GeoInLonLat = "
	!PRINT *, SHAPE(GeoInLonLat)
	!PRINT *, "lon = ", GeoInLonLat(:,:,1)
	!PRINT *, "lat = ", GeoInLonLat(:,:,2)
	
	!-------------------------------------------------------------------------------
	! main loop over the different variables, namelist controlled
	
	ALLOCATE ( frequency(3) )
	frequency(1) = "3hr"
	frequency(2) = "6hr"
	frequency(3) = "day"
	frequency(4) = "mon"
	frequency(5) = "sem"
	!frequency(6) = "fx"
	!frequency(7) = "1hr"
	
	!DEFAULT VALUES
	obs_interval = 3 !!!! IMPORTANT -- EVERY HOW MANY HOURS DO WE HAVE AN OBSERVATION ???? --  !!!!
	!isEndOfDay = 1 !Not really needed...
	isEndOfMonth = 1
	isEndOfSeason = 1
	
	
	ALLOCATE ( fnNMLvar(7) )
	fnNMLvar(1) = "runctrl.vars.nml"
	fnNMLvar(2) = "runctrl.vars.nml_evp_roff"
	fnNMLvar(3) = "runctrl.vars.nml_water_column"
	fnNMLvar(4) = "runctrl.vars.nml_vars_on_plevels"
	fnNMLvar(5) = "runctrl.vars.nml_pr_mrso"
	fnNMLvar(6) = "runctrl.vars.nml_snow"
	fnNMLvar(7) = "runctrl.vars.nml_radiation"
	!fnNMLvar(8) = "runctrl.vars.nml_cape"
	!fnNMLvar(9) = "runctrl.vars.nml_pr_tas_1hr_test"
	
	
	!DO ifrq = 1, SIZE(frequency), 1
	!DO ifrq = 1, 1, 1
	!ifrq = 4   !Commented from here, for the sake of giving life to it through command-line arguments. 
	
	
	
	!IMPORTANT NOTE:
	!If input data does NOT include last day of a leap year's February (i.e. 29th), then set dataIncludesLeapYearDays to 0,
	!otherwise set it to 1.
	dataIncludesLeapYearDays = 1
	
	
	
	PRINT *, "============================================================"
	PRINT *, "freq = ", frequency(ifrq)
	
	!------------------------------------------------------------------------------
	
	SELECT CASE (frequency(ifrq))
		CASE ('1hr')
			
			cell_methods(:) = cm1hr(:)
			dtHours = 1.
		
		CASE ('3hr')
			
			cell_methods(:) = cm3hr(:)
			dtHours = 3.
		
		CASE ('6hr')
			
			cell_methods(:) = cm6hr(:)
		
		CASE ('day')
			
			cell_methods(:) = cmDay(:)
			dtHours = 24.
		
		CASE ('mon')
			
			cell_methods(:) = cmMon(:)
			dtHours = 720.	
		
		CASE ('sem')
			
			cell_methods(:) = cmSea(:)
			dtHours = 2880.	
		
		CASE DEFAULT
			
			PRINT *, "invalid time interval specified"
			STOP
			
		END SELECT
		
		
		!------------------------------------------------------------------------------
		!get a file list of all wrfout and wrfxtrm files
		!use regex to refine the ls output and the filelist
		!non-std, works for gfortran (fct & subroutine) + ifort
		
		PRINT *, "============================================================"
		PRINT *, "*** FILELIST CREATION ***"
		
		tmpfileFL = "tmpfileFL"
		
		!PRINT *, "filelist search pattern = ", TRIM(DirInputSimResRoot) // "/" // TRIM(domain) // "/" // "*/*wrfout*nc"
		!CALL SYSTEM("ls -1 " // TRIM(DirInputSimResRoot) // "/" // TRIM(domain) // "/*/*wrfout*{" // ts // ".." // te // "}*nc > " // tmpfileFL)
		
		!PRINT *, "filelist search pattern = ", TRIM(DirInputSimResRoot) // "/" // "wrfout*" //TRIM(domain)
		!CALL SYSTEM("ls -1 " // TRIM(DirInputSimResRoot) // "/" // "wrfout*" //TRIM(domain)//"*{" // ts // ".." // te // "}* > " // tmpfileFL)
	

		!dirLen is the length of the input data path up until the last folder before actual wrf input data is located.
		!For example, if wrfout input files are located in path = /mnt/meteo/groups/hindcast_WRF371b/XXXX_WRFv3.7.1_20years/wrfout/
		!then dirLen = '65'.
		!For now, it should be changed manually...
		dirLen = '58'
		
		!Make sure that DirInputSimResRoot is last folder BEFORE ALL YEAR DATA FOLDERS.
		!For example, if path of input wrfout data = /mnt/meteo/groups/hindcast_WRF371b/XXXX_WRFv3.7.1_20years/wrfout/
		!then DirInputSimResRoot = /mnt/meteo/groups/hindcast_WRF371b
		!Then, change filesearch pattern (system call) accordingly to include folders from all years inside DirInputSimResRoot.
		!This is for looking into multiple directories in order to be able to make  computations for consequent years at a single run...
		!i.e.:
		!ls -1 /mnt/meteo/groups/hindcast_WRF371b/{1989-2008}_WRFv3.7.1_20years/wrfout/wrfout_d01_* | sort -r | sort -u -k1.65  > tmpfileFL
		PRINT *, "filelist search pattern = ", TRIM(DirInputSimResRoot) // "/" // "*_WRFv3.7.1_GISS" // "/" // "wrfout_" //TRIM(domain) // "_*"
		CALL SYSTEM("ls -1 " // TRIM(DirInputSimResRoot) // "/" // "*{" // ts // ".." // te // "}*"  // "_WRFv3.7.1_GISS" // "/" // "wrfout_" // TRIM(domain) // "_*" // "| sort -r | sort -u -k1." // dirLen // " >  tmpfileFL")
		
	
		ft = 0
		counter = 0
		!Convert ts and te from string to integers for computational purposes...
		READ(ts(1:4),'(i4)') FirstYearOfRun
		READ(te(1:4),'(i4)') LastYearOfRun
		
		
		CALL generateFilelist
		
		!PRINT *, "filelist search pattern = ", TRIM(DirInputSimResRoot) // "/" // TRIM(domain) // "/" // "*/*wrfxtrm*nc"
		!CALL SYSTEM("ls -1 " // TRIM(DirInputSimResRoot) // "/" // "wrfxtrm*" //TRIM(domain)//"*{" // ts // ".." // te // "}* > " // tmpfileFL)
		!ft = 1
		!CALL generateFilelist
		
		DO i=1,SIZE(fl_wrfout(:)),1
			
			PRINT '(100A)', fl_wrfout(i)
			!  PRINT '(100A)', fl_wrfxtr(i)
			
		END DO
		
		!------------------------------------------------------------------------------
		!creation of the main reference array
		
		PRINT *, "============================================================"
		PRINT *, "*** TIME REFERENCE ARRAY ***"
		
		CALL CreateRefTimeArray( frequency(ifrq), dataIncludesLeapYearDays )
		
		PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
		SHAPE(TimeRefArray)
		PRINT *, "SIZE(TimeRefArray,1)",SIZE(TimeRefArray,1) 
		PRINT *, "SHAPE(TimeRefArray,1)",  SHAPE(TimeRefArray,1)
		
		DO i = 1,605,1
			
			PRINT *, TimeRefArray(i,1:5)
			
		END DO
		
		
		!CHANGES IN THIS LOOP AND THE NEXT ONE FOR ARIS EXPERIMENT
		!-------------------------------------------------------------------------------
		! loop over the different variables
		DO varnml = 1, 7, 1 !loop over different var namelists (not best solution, but one namelist for all vars is to big)
			
		!DO varnml = first_varnml, first_varnml, 1 !choose just specific namelists from list above
			
			OPEN(2,FILE=TRIM(fnNMLvar(varnml)))
			READ(UNIT=2,NML=vars)
			CLOSE(2)
			
			SELECT CASE (frequency(ifrq))
				CASE ('1hr')
					
					cell_methods(:) = cm1hr(:)
				
				CASE ('3hr')
					
					cell_methods(:) = cm3hr(:)
				
				CASE ('6hr')
					
					cell_methods(:) = cm6hr(:)
				
				CASE ('day')
					
					cell_methods(:) = cmDay(:)
				
				CASE ('mon')
					
					cell_methods(:) = cmMon(:)
				
				CASE ('sem')
					
					cell_methods(:) = cmSea(:)
				
				CASE DEFAULT
					
					PRINT *, "invalid time interval specified"
					STOP
					
				END SELECT
				
				
				SELECT CASE (TRIM(fnNMLvar(varnml)))
					CASE ("runctrl.vars.nml")
						
						nvar_nml = 9
					
					CASE ("runctrl.vars.nml_evp_roff")
						
						nvar_nml = 4
					
					CASE ("runctrl.vars.nml_water_column")
						
						nvar_nml = 3
					
					CASE ("runctrl.vars.nml_vars_on_plevels")
						
						nvar_nml = 12
					
					CASE ("runctrl.vars.nml_pr_mrso")
						
						nvar_nml = 4
					
					CASE ("runctrl.vars.nml_snow")
						
						nvar_nml = 6 
					
					CASE ("runctrl.vars.nml_radiation")
						
						nvar_nml = 9
						!CASE ("runctrl.vars.nml_radiation_alternative")
						
						!nvar_nml = 9
						
						!CASE ("runctrl.vars.nml_cape")
						!nvar_nml = 1
						
						!CASE ("runctrl.vars.nml_pr_tas_1hr_test")
						!nvar_nml = 3
						
					END SELECT
					
					print*, "nvar_nml", nvar_nml
					
					DO ivar = 1, nvar_nml, 1
					!DO ivar = firstvar_nvar_nml, firstvar_nvar_nml, 1 !choose just specific variables from namelist
					
						!DO ivar = 1, 1, 1    !choose just specific variables from namelist (look up var entry in individual namelist)
						
						PRINT *,"============================================================"
						PRINT *, "*** ", TRIM(var_cmip(ivar)), " ***"
						
						!First process variables that do not require opening pairs of files (temporal averages and cumulatives.)
						!and then process the rest of them, that require opening pairs of files (temporal averages and cumulatives.)
						
						
						!-------------------------------------------------------------------------------
						! loop over the filelist per variable
						! content of filelist is defined by filename patterns in system call
						
						
						IF (var_cmip(ivar) == "pr" .OR. var_cmip(ivar) == "prc" .OR. &
						var_cmip(ivar) == "prsn" .OR. var_cmip(ivar) == "snm" .OR. &
						var_cmip(ivar) == "evspsbl" .OR. var_cmip(ivar) == "evspsblpot" .OR. &
						var_cmip(ivar) == "mrros" .OR. var_cmip(ivar) == "mrro" .OR. &
						var_cmip(ivar) == "rsds" .OR. var_cmip(ivar) == "rlds" .OR. &
						var_cmip(ivar) == "rsus" .OR. var_cmip(ivar) == "rlus" .OR. &
						var_cmip(ivar) == "rsdt" .OR. var_cmip(ivar) == "rlut" .OR. &
						var_cmip(ivar) == "hfss" .OR. var_cmip(ivar) == "rsut" .OR. &
						var_cmip(ivar) == "hfls" .OR. var_cmip(ivar) == "mrso") THEN
							
							!Perform only the steps necessary to read and process the above variables...
							!and NOT in an iterative manner, since they only require instantaneous readings.
							
							!Because the calculation takes place by reading two files (that correspond to the 
							!temporal instants in time, which signify the requested interval that is desired to
							!be averaged, when the "time" comes, when the interval exceeds the time period spanned
							!by the available data (as set by ts and te variables in the configuration external file),
							!the process will be stopped. So the last average will always be missing. This means:
							!Daily averages --> Last day's average will not be available (since the beginning of
							!the next day is not available). Same goes for month and season... last ones are not calculated.
							
							print *, "Cumulative Variable PROCESSING"
							
							!Designate the time step for the file reading process...
							record_count = dtHours / obs_interval
							!if record_count is NOT equal to 1, then some type of temporal average has been requested...
							
							IF (record_count == 1) THEN !dtHours == obs_interval (maybe 3-hr)
								
								!do for every file....
								DO ifl = 1, SIZE(fl_wrfout), 1 ! operational: loop over complete filelist
									
									print *,' SIZE(fl_wrfout', SIZE(fl_wrfout)
									
									!     DO ifl = 1, SIZE(fl_wrfxtr), 1 ! operational: loop over complete filelist
									!      print *,' SIZE(fl_wrfout', SIZE(fl_wrfxtr)
									
									!DO ifl = 1, 1, 1 ! testing: loop over specific entry in filelist (e.g. just January)
									
									CALL CPU_TIME(cpuTs)
									
									PRINT *, "------------------------------------------------------------"
									
									IF ( filetype(ivar) == "s" ) THEN
										
										iflWRFin = fl_wrfout(ifl)
									
									ELSE IF ( filetype(ivar) == "x" ) THEN
										
										iflWRFin = fl_wrfxtr(ifl)
										
									END IF
									
									PRINT '(100A)', TRIM(iflWRFin)
									
									!-------------------------------------------------------------------------------
									! which timespan is covered by the WRF outputs?
									! assume timespan wrfout = timespan wrfxtrm
									! this determines how many times the tool has to loop over the inputs
									! also check how many years are covered by a single wrfout and wrfxtrm which
									! determines the output file generation
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncid_in)
									sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
									sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
									sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
									NAME = InDimNameRec, LEN = InDimLenRec)
									ALLOCATE(InVarDataRec(InDimLenRec))
									sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
									sts = NF90_CLOSE(ncid_in)
									
									print *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)
									
									ALLOCATE(InDateTimeYear(InDimLenRec))
									ALLOCATE(InDateTimeMonth(InDimLenRec))
									ALLOCATE(InDateTimeDay(InDimLenRec))
									ALLOCATE(InDateTimeHour(InDimLenRec))
									
									! this is the temporal coverage of the WRF input data
									DO i = 1, SIZE(InVarDataRec), 1
										
										READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(i), InDateTimeMonth(i), InDateTimeDay(i), InDateTimeHour(i)
										
									END DO
									
									
									!-------------------------------------------------------------------------------
									! check whether a file is needed at all
									! check for first date
									! for subsequent dates, just check whether it is the same as the previous one
									! if this is not the case, then check whether file exists...
									! if it does not exist (default case): create
									! pathname is always needed
									
									! loop over the individual timesteps in the WRF files...
									! this may take some time but it is robust and also extremely large arrays
									! might be used
									! highly robust code
									
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ARIS EDIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									! Allocate correspondingly the temp_data_InTime variable
									
									
									!	SELECT CASE (frequency(ifrq))
									!		CASE ('3hr')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,1), STAT = sts) !run through the entire list of 3-hour records and use them all for the calculations.
									!		CASE ('day')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,8), STAT = sts) !only run once... the daily data will be processed once..
									!		CASE DEFAULT
									!			PRINT *, "something happened... "
									!			total_iters = InDimLenRec !make sure bad results in default behavior (i.e. what was happening anyway before...)
									!			!STOP
									!	END SELECT
									
									
									ALLOCATE( temp_data_InTime(xfocus,yfocus,InDimLenRec), STAT = sts)
									ALLOCATE( counter_array(InDimLenRec), STAT = sts)
									
									IF ( CheckForLeapyear( InDateTimeYear(1), dataIncludesLeapYearDays ) == 366 ) THEN
										
										NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
									
									ELSE
										
										NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
										
									END IF
									
									IF (isEndOfMonth) THEN
										
										IF (ALLOCATED(temp_data_Monthly)) THEN
											
											DEALLOCATE (temp_data_Monthly)
											
										END IF
										
										ALLOCATE( temp_data_Monthly(xfocus,yfocus, NoOfDaysPerMonth(InDateTimeMonth(1))  ), STAT = sts  )
										isEndOfMonth = 0
										
									END IF
									
									IF (isEndOfSeason) THEN
										
										IF (ALLOCATED(temp_data_Seasonally)) THEN
											
											DEALLOCATE (temp_data_Seasonally)
											
										END IF
										
										ALLOCATE( temp_data_Seasonally(xfocus, yfocus, 3, 31 ), STAT = sts)
										isEndOfSeason = 0
										dayCount = 0
										
									END IF
									
									!!!! Here is the case, when 3-hr intervals have been requested... so no alterations are made,
									!!!! except for removing the irrelevant variables, keeping only the cumulative ones..
									DO it = 1, InDimLenRec, 1
										
										PRINT *, "----------------------------------------"
										PRINT *, "working on date in WRF input file = ", TRIM(InVarDataRec(it))
										PRINT *, "---"
										
										!-------------------------------------------------------------------------------
										! generate path and filename
										! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/EUR-44/CRPGL/
										! ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
										! evspsbl_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_CRPGL-WRFARW331_v1_3hr_
										! 1989010100-1989123121
										
										! output does not yet exist
										! monthy check is basically not even possible, but does not do any harm
										! this is a rerstriction for all those who might have a different file
										! structure
										IF ( InDateTimeYearPrev /= InDateTimeYear(it) ) THEN !.AND. &
											!( InDateTimeMonthPrev /= InDateTimeMonth(it) ) ) THEN
											
											InDateTimeYearPrev = InDateTimeYear(it)
											
											PRINT *, "start of processing or new year encountered -> t ref. vec. and filecheck"
											
											!READ( InDateTimeYear(it), '(4A)' ) InDateTimeYearStr
											!READ( InDateTimeMonth(it), '(2A)' ) InDateTimeMonthStr
											WRITE (InDateTimeYearStr,'(I4.4)') InDateTimeYear(it)
											WRITE (EndInDateTimeYearStr,'(I4.4)') InDateTimeYear(it)+1
											
											
											
											
											!WRITE (InDateTimeMonthStr,'(I2.2)') InDateTimeMonth(it)
											!PRINT *, InDateTimeYearStr, InDateTimeMonthStr
											
											PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
											SHAPE(TimeRefArray)
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												EndDate=EndInDateTimeYearStr//"010100"
											
											ELSE
												
												EndDate=InDateTimeYearStr//"123121"
												
											END IF
											
											
											pn_out = TRIM(product)                       // "/" // &
											TRIM(CORDEX_domain)                 // "/" // &
											TRIM(institute_id)                  // "/" // &
											TRIM(driving_model_id)              // "/" // &
											TRIM(driving_experiment_name)       // "/" // &
											TRIM(driving_model_ensemble_member) // "/" // &
											TRIM(model_id)                      // "/" // &
											TRIM(rcm_version_id)                // "/" // &
											TRIM(frequency(ifrq))               // "/" // &
											TRIM(var_cmip(ivar))
											
											fn_out = TRIM(var_cmip(ivar))                // "_" // &
											TRIM(CORDEX_domain)                 // "_" // &
											TRIM(driving_model_id)              // "_" // &
											TRIM(driving_experiment_name)       // "_" // &
											TRIM(driving_model_ensemble_member) // "_" // &
											TRIM(model_id)                      // "_" // &
											TRIM(rcm_version_id)                // "_" // &
											TRIM(frequency(ifrq))               // "_" // &
											!                   InDateTimeYearStr//"010100-"//InDateTimeYearStr//"123121" // &
											InDateTimeYearStr//"010100-"//TRIM(EndDate)//  &
											".nc"
											
											PRINT *, "pn_out = ", TRIM(pn_out)
											PRINT *, "fn_out = ", TRIM(fn_out)
											
											!-------------------------------------------------------------------------------
											! extract the time info from the ref array which fits the respective year
											! ...as there is no "WHERE" the way I need it in F95, use loops
											! this is needed whenever a new netcdf file is to be used, also if
											! this file exists already
											!        READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
											
											PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
											SHAPE(TimeRefArray)
											
											DEALLOCATE( TimeRefArraySelYear )
											
											print *, "DEALLOCATE( TimeRefArraySelYear )" 
											
											DEALLOCATE( Time_bnds ) !SKn 
											
											print *, "DEALLOCATE( Time_bnds )"
											
											counter = 0
											PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
											PRINT *,'SHAPE(TimeRefArray, 1)',  SHAPE(TimeRefArray, 1)
											PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
											PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
											PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)
											DO i = 1, SIZE(TimeRefArray, 1), 1
												
												IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
													
													counter = counter + 1
													
												END IF
												
											END DO
											PRINT *, "Counter: ", counter
											! holds data of exactly 1 year
											PRINT *, "timesteps in the time ref. subset = ", counter
											ALLOCATE( TimeRefArraySelYear( counter, 5 ) ) ! index, y, m, d, h
											
											print *, "ALLOCATE( TimeRefArraySelYear( counter, 5 ) )"
											
											print *, "TimeRefArraySelYear( counter, 1 )", TimeRefArraySelYear( counter, 1 )
											print *, "TimeRefArraySelYear( counter, 2 )", TimeRefArraySelYear( counter, 2 )
											print *, "TimeRefArraySelYear( counter, 3 )", TimeRefArraySelYear( counter, 3 )
											print *, "TimeRefArraySelYear( counter, 4 )", TimeRefArraySelYear( counter, 4 )
											print *, "TimeRefArraySelYear( counter, 5 )", TimeRefArraySelYear( counter, 5 )
											
											! find the matching elements of the respecitve year and copy them
											
											ALLOCATE( Time_bnds( 2, counter ) )
											
											print *, "ALLOCATE( Time_bnds( 2, counter ) )", Time_bnds( 2, counter )
											
											counter = 0
											DO i = 1, SIZE(TimeRefArray, 1), 1
												
												IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
													
													counter = counter + 1
													TimeRefArraySelYear(counter,1:5) = TimeRefArray(i,1:5)
													
												END IF
												
											END DO
											
											!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySelYear(:,:) )
											
											!-------------------------------------------------------------------------------
											! check for existance of the file and generate file if needed
											
											INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )
											
											! could exist already from a previous run of tool and due to multiple
											! months in a file (i.e. 1 WRF output file may cover several months)
											IF ( FileExists ) THEN
												
												PRINT *, "path and file exist, continue filling"
											
											ELSE
												
												PRINT *, "path and file do not yet exist, create path and NetCDF file first"
												PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)
												
												CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )
												
												!-------------------------------------------------------------------------------
												! non-std, works for gfortran (fct & subroutine) + ifort
												! comment lines in the NetCDF file global attribute definition
												! turn standard checking in Makefile off
												! trackingID = "xxxxxxxx-xxxx-Mxxx-Nxxx-xxxxxxxxxxxx"
												! creationDate = "YYYY-MM-DD-THH:MM:SSZ"
												
												CALL SYSTEM(cmdUUID)
												OPEN(1,FILE="tmpfileUUID",STATUS='old')
												READ(1,*) trackingID
												CLOSE(1)
												PRINT *, "uuidgen externally generated trackingID = ", trackingID
												
												CALL SYSTEM(cmdDate)
												OPEN(1,FILE="tmpfileDate",STATUS='old')
												READ(1,*) creationDate
												CLOSE(1)
												PRINT *, "date externally generated creation date = ", creationDate
												
												!-------------------------------------------------------------------------------
												! create NetCDF file
												! NF90_CLASSIC_MODEL = NetCDF4_classic
												! NF90_HDF5 = NetCDF4 based on HDF5
												! NF90_CLOBBER = old NetCDF
												! sts = NF90_CREATE(PathFileNameOutTEST, NF90_HDF5, ncid)
												
												!comb_flags = IOR(NF90_HDF5, NF90_CLASSIC_MODEL)
												!https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/NF90_005fCREATE.html
												!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid)   !not sure whether this is the right data format spec. I guess it may be right using compression but not the other fancy stuff from NetCDF4
												!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncid)   !if anything, then use this here
												sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_HDF5, ncid)
												
												
												! always included
												sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
												sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid)
												IF ( height(ivar) /= -999 ) THEN
													
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
														
														sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
													
													ELSE
														
														sts = NF90_DEF_DIM(ncid, "plev", 1, height_dimid)
														
													END IF
													
												END IF
												sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
												
												! for mean variables
												IF ( cell_methods(ivar) == "mean" ) THEN
													
													sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
													
												END IF
												
												
												! always included
												sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lon_varid)
												sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
												sts = nf90_put_att(ncid, lon_varid, "long_name", "longitude")
												sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
												
												! always included
												sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lat_varid)
												sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
												sts = nf90_put_att(ncid, lat_varid, "long_name", "latitude")
												sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
												
												! always included
												sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid)
												sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
												sts = nf90_put_att(ncid, rlon_varid, "long_name", "longitude in rotated pole grid")
												sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
												sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
												
												! always included
												sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid)
												sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
												sts = nf90_put_att(ncid, rlat_varid, "long_name", "latitude in rotated pole grid")
												sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
												sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")
												
												! always included
												! restriction to one domain only
												sts = nf90_def_var(ncid, "rotated_pole", NF90_DOUBLE, rotated_pole_varid)
												sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
												sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", 39.25)
												sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", -162.0)
												
												! depends whether height is set in the nml
												IF ( height(ivar) /= -999 ) THEN
													
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
														
														sts = nf90_def_var(ncid, "height", NF90_DOUBLE, height_varid)
														sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
														sts = nf90_put_att(ncid, height_varid, "long_name", "height")
														sts = nf90_put_att(ncid, height_varid, "units", "m")
														sts = nf90_put_att(ncid, height_varid, "positive", "up")
														sts = nf90_put_att(ncid, height_varid, "axis", "Z")
													
													ELSE
														
														sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, height_varid)
														sts = nf90_put_att(ncid, height_varid, "standard_name","air_pressure")
														sts = nf90_put_att(ncid, height_varid, "long_name", "pressure")
														sts = nf90_put_att(ncid, height_varid, "units", "Pa")
														sts = nf90_put_att(ncid, height_varid, "positive", "down")
														sts = nf90_put_att(ncid, height_varid, "axis", "Z")                 
														
													END IF
													
												END IF
												
												!missing: lvl, depends
												
												! always included
												sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid)
												sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
												sts = nf90_put_att(ncid, rec_varid, "long_name", "time")
												sts = nf90_put_att(ncid, rec_varid, "units", "days since 1949-12-01T00:00:00Z")
												sts = nf90_put_att(ncid, rec_varid, "calendar", "standard")
												sts = nf90_put_att(ncid, rec_varid, "axis", "T")
												
												! for mean variables
												IF ( cell_methods(ivar) == "mean" ) THEN
													! Add time bounds to time varaible
													
													sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")	
													sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid)
													sts = nf90_put_att(ncid, recbnds_varid, "standard_name", "time")
													sts = nf90_put_att(ncid, recbnds_varid, "long_name", "time")
													sts = nf90_put_att(ncid, recbnds_varid, "units", "days since 1949-12-01T00:00:00Z")
													sts = nf90_put_att(ncid, recbnds_varid, "calendar", "standard")
													!sts = nf90_put_att(ncid, recbnds_varid, "axis", "T")
													
													print *,'rec_varid', rec_varid
													print *,'recbnds_varid', recbnds_varid
													
												END IF
												
												! always included
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact", contact)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "creation_date", creationDate)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment", experiment)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment_id", experiment_id)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment", driving_experiment)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_id", driving_model_id)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_ensemble_member", driving_model_ensemble_member)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment_name", driving_experiment_name)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "frequency", frequency(ifrq))
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution", institution)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_id", institute_id)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "model_id", model_id)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "forcing",forcing)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact",contact)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source",source)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "rcm_version_id", rcm_version_id)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "CORDEX_domain", CORDEX_domain)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id", trackingID)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "title",title)
												sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_run_id", institute_run_id)
												
												! always included
												IF ( height(ivar) /= -999 ) THEN
													
													sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, height_dimid, rec_dimid /), x_varid)
												
												ELSE
													
													sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid)
													
												END IF
												
												sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
												sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
												sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
												
												IF ( positive(ivar) /= '-999' ) THEN
													
													sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
													
												END IF
												
												sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
												
												IF ( height(ivar) /= -999 ) THEN
													
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
														
														sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat height")
													
													ELSE
														
														sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat plev")
														
													END IF
												
												ELSE
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat")
													
												END IF
												
												sts = nf90_put_att(ncid, x_varid, "grid_mapping", "rotated_pole")
												
												!IF ( height(ivar) == 850 ) THEN
												sts = nf90_put_att(ncid, x_varid, "missing_value", 1.e20)
												sts = nf90_put_att(ncid, x_varid, "_FillValue", 1.e20)
												!END IF
												
												sts = NF90_ENDDEF(ncid)
												
												IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
													
													NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
												
												ELSE
													
													NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
													
												END IF
												
												IF (dtHours == 720) THEN
													
													HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
												
												ELSE IF (dtHours == 2880) THEN
													
													IF (TimeRefArraySelYear(it,3) == 12 .OR. &
													TimeRefArraySelYear(it,3) == 1 .OR. &
													TimeRefArraySelYear(it,3) == 2) THEN
														
														HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
														NoOfDaysPerMonth(12)) * 24 !It's Winter..
														!adding forthcoming december, but, still.. it's 31 days anyway
													
													ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
													TimeRefArraySelYear(it,3) == 4 .OR. &
													TimeRefArraySelYear(it,3) == 5) THEN
														
														HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
														NoOfDaysPerMonth(5)) * 24 !It's Spring..
													
													ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
													TimeRefArraySelYear(it,3) == 7 .OR. &
													TimeRefArraySelYear(it,3) == 8) THEN
														
														HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
														NoOfDaysPerMonth(8)) * 24 !It's Summer..
													
													ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
														!TimeRefArraySelYear(it,3) == 10 .OR. &
														!TimeRefArraySelYear(it,3) == 11) THEN
														
														HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
														NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
														
													END IF
												
												ELSE
													
													HoursOfCurrentInterval = dtHours
													
												END IF
												
												
												! add time, whole year from above
												IF ( cell_methods(ivar) == "point" ) THEN
													
													print*, 'cell_methods:', cell_methods(ivar)
													sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
													
													print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
													print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
													print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
													print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
													print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
													
												END IF
												
												IF ( cell_methods(ivar) == "mean" ) THEN
													
													print*, 'cell_methods:', cell_methods(ivar)
													sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.) )
													
													print *, 'sts NF90_PUT_VAR time', sts
													
													print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
													print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
													print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
													print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
													print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
													
													Time_bnds(1,:) = TimeRefArraySelYear(:,1)
													Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
													
													print*, 'recbnds_varid', recbnds_varid
													
													sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
													
												END IF
												!print *,'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												
												
												sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
												START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
												sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
												START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
												sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
												sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
												
												! add time_bnds, calc here
												
												! add height from NML
												IF ( height(ivar) /= -999 ) THEN
													
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
														
														sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
													
													ELSE
														
														sts = NF90_PUT_VAR(ncid, height_varid, height(ivar)*100. )
														
													END IF
													
												END IF
												
												sts = NF90_CLOSE(ncid)
												
											END IF ! file exists y/n
											
										END IF ! checking with previous year and month
										
										!-------------------------------------------------------------------------------
										! match timestep of WRFin with the subset of the ref time vec which belong to
										! the NetCDF file of the year currently open to receive data
										! NOT SURE WHETHER THIS IS NEEDED AT ALL, THIS IS THE WRONG DIRECTION ???
										! see whether the current time of the timestep fits anywhere in the
										
										! ARIS CHANGE --------------- This small excerpt below is used to locate where
										! in the final nc file is the data to be written (i.e. in position "counter").
										! However, due to the nature of the TimeRefArraySelYear variable, the original
										! WRF data and the aforementioned variable do not have the same frequency and
										! instants in time as records. This incongruence completely destroys the purpose
										! of computing averages. As a result, this excerpt here must ONLY be executed if
										! the WRF input data and the requested output frequency coincide (i.e. no temporal)
										! averaging takes place.
										
										IF (record_count.eq.1) THEN
											
											PRINT *, "reading WRF sim. res. = ", TRIM(InVarDataRec(it)), it
											!PRINT *, SIZE(TimeRefArraySelYear,1)
											!PRINT *, SHAPE(TimeRefArraySelYear)
											!PRINT *, "current transferred input time: ", InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
											
											counter = 0
											DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
												
												counter = counter + 1
												
												!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
												
												IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(it)  ) .AND. &
												( TimeRefArraySelYear(i,3) == InDateTimeMonth(it) ) .AND. &
												( TimeRefArraySelYear(i,4) == InDateTimeDay(it)   ) .AND. &
												( TimeRefArraySelYear(i,5) == InDateTimeHour(it)  ) ) THEN
													
													EXIT
													
												END IF
												
											END DO
											
											PRINT *, "index where in the NC file the WRF data is sorted in = ", &
											counter
!SUPER IMPORTANT ADDITION TOO
											counter_array(it) = counter
										END IF
										! If this is not run here, then it must be run directly after the DO it = ... loop
										! is finished (somewhere in line 1900+)...
										
										
										!-------------------------------------------------------------------------------
										! read orig WRF outpoutses
										! there is always a corresponding time-slot in the NC file
										! extracted time from above
										! "it" controls it all: timestep in the individual WRF file
										
										sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)
										
										!!!!!!!!!!!!! ONLY READ THE VARIABLES THAT PROVIDE CUMULATIVE WRF OUTPUTS.......!!!
										IF (var_cmip(ivar) == "pr") THEN 
											
											ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
											ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
											sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
											
											sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )   !read two timesteps to calculate 3hr sum
											
											sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
										
										ELSE IF (var_cmip(ivar) == "prc") THEN
											
											ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
											
											sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
										
										ELSE IF (var_cmip(ivar) == "prsn") THEN
											
											ALLOCATE( snownc_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
											
											sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
										
										ELSE IF (var_cmip(ivar) == "snm") THEN
											
											ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
											
											sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
											
										
										ELSE IF (var_cmip(ivar) == "evspsbl") THEN
											
											ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
											
											print*, "LETS SEE IF SFCEVP variable was found:", sts
											
											sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
										ELSE IF (var_cmip(ivar) == "evspsblpot") THEN
											
											ALLOCATE( potevp_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
											
											sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
										
										ELSE IF (var_cmip(ivar) == "mrros") THEN
											
											ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
											
											sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
										ELSE IF (var_cmip(ivar) == "mrro") THEN
											
											ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
											ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
											sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
											
											sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
											sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
											
										
										ELSE IF ((var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")  &
										.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
										.or. (var_cmip(ivar) == "rlut")                                 &
										.or. (var_cmip(ivar) == "rsdt") .or. (var_cmip(ivar) == "rsut") &
										.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
											
											ALLOCATE( rad_in ( xfocus, yfocus, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
											
											sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,:), &
											START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
											
											
											print*, var_cmip(ivar), rad_in(50,50,1), rad_in(50,50,2)
											print*, 'difference in J m-2', (rad_in(50,50,2) - rad_in(50,50,1))
											print*, 'in mean W m-2', (rad_in(50,50,2) - rad_in(50,50,1))/ (obs_interval*3600.)
											
											! alternative: since accumulated values as read above get so large in 
											! long term simulations that their differences loose accuracy, use 
											! instantaneous values instead and calculate means
											
										
										ELSE IF (var_cmip(ivar) == "mrso") THEN
											
											ALLOCATE( smois_in( xfocus, yfocus, 4, 2 ), STAT=sts )
											
											sts = NF90_INQ_VARID(ncidin, "SMOIS", varid)
											
											sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,:), &
											START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, 4, 2 /) )
											
											
										END IF
										
										sts = NF90_CLOSE(ncidin)
										
										!-------------------------------------------------------------------------------
										! some analysis of the data
										
										print *,"shape of array" , SHAPE(temp_data)
										print *,"size of array" , SIZE(temp_data)
										!
										!       stat_mean = SUM(temp_data(:,:,5))/(MAX(1,SIZE(temp_data(:,:,5))))
										!       PRINT *, stat_mean
										stat_mean = SUM(temp_data(:,:))/SIZE(temp_data(:,:))
										PRINT *, stat_mean
										
										!-------------------------------------------------------------------------------
										! processing
										
										
										!       ***pr***
										IF (var_cmip(ivar) == "pr") THEN 
											
											temp_data(:,:) = ((rainnc_in(:,:,2) + rainc_in(:,:,2)) - (rainnc_in(:,:,1) + rainc_in(:,:,1)))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
											!ATTENTION: implement adjustable time intervals that the differences are devided by
											
										END IF
										
										
										!       ***prc***
										IF (var_cmip(ivar) == "prc") THEN
											
											temp_data(:,:) = (rainc_in(:,:,2) - rainc_in(:,:,1))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
											
										END IF
										
										!       ***prsn***
										IF (var_cmip(ivar) == "prsn") THEN
											
											temp_data(:,:) = (snownc_in(:,:,2) - snownc_in(:,:,1))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
											
										END IF
										
										!       ***snm***
										IF (var_cmip(ivar) == "snm") THEN
											!print*, "STARTSRTTASTARST"
											!print*, "acsnom_in(:,:,2):", acsnom_in(:,:,2)
											!print*, "acsnom_in(:,:,1):", acsnom_in(:,:,1)
											
											temp_data(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1))/(obs_interval*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
											
											!print*, "temp_data(:,:)", temp_data(:,:)
											!print*, "OVEROVEROVEROVOVER"
										END IF
										
										!       ***evspsbl***
										IF (var_cmip(ivar) == "evspsbl") THEN
											
											print*, "sfcevp_in(:,:,2):", sfcevp_in(:,:,2)
											print*, "sfcevp_in(:,:,1):", sfcevp_in(:,:,1)
											temp_data(:,:) = (sfcevp_in(:,:,2) - sfcevp_in(:,:,1))/(obs_interval*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
											print*, "temp_data(:,:):", temp_data(:,:)
										END IF
										
										!       ***evspblpot**
										IF (var_cmip(ivar) == "evspsblpot") THEN
											
											temp_data(:,:) = (potevp_in(:,:,2) - (potevp_in(:,:,1)))/L   !unit [W m-2]/[J kg-1] -> [kg m-2 s-1]
											
											! THERE IS STH WRONG WITH THE UNITS: WRF's POTEVP is accumulated and declared to be in W m-2. 
											! It doesen't make sense to accumulate in W m-2, but even if assume it as J m-2 or derive kg m-2 
											! by using latent heat of vaporization you never get values that have a reasonable magnitude...
											
											
										END IF
										
										!       ***mrros***
										IF (var_cmip(ivar) == "mrros") THEN
											
											temp_data(:,:) = (sfroff_in(:,:,2) - sfroff_in(:,:,1))/(obs_interval*3600.)       !unit [mm/3hr] to [kg m-2 s-1]
											
										END IF
										
										!       ***mrro***
										IF (var_cmip(ivar) == "mrro") THEN
											
											temp_data(:,:) = ((sfroff_in(:,:,2) - sfroff_in(:,:,1)) + (udroff_in(:,:,2) - udroff_in(:,:,1)))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
											
										END IF
										
										!!!!!!! This if has probably forgotten to include 3 variables (rlut, rsdt, rsut).
										!!!!!! They are read from the files earlier altogether, but have not been included here.
										!!!!!!!! Probably a mistake... add them some time...
										!       ***rsds, rlds, rsus, rlus***
										IF ( (var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")      &
										.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
										.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
											
											IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation") THEN
												
												temp_data(:,:) = (rad_in(:,:,2) - rad_in(:,:,1)) /(obs_interval*3600.)       ! take difference of accumulated values
											
											ELSE IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation_alternative") THEN
												
												temp_data(:,:) = (rad_in(:,:,2) + rad_in(:,:,1)) / 2.              ! take mean of instantaneous values
												
											END IF 
											
										END IF
										
										
										!       ***mrso***
										IF (var_cmip(ivar) == "mrso") THEN
											
											temp_data(:,:) = (((smois_in(:,:,1,1)*0.1 + smois_in(:,:,2,1)*0.3 + smois_in(:,:,3,1)*0.6 + smois_in(:,:,4,1)*1.0 ) + &
											(smois_in(:,:,1,2)*0.1 + smois_in(:,:,2,2)*0.3 + smois_in(:,:,3,2)*0.6 + smois_in(:,:,4,2)*1.0 ))/2.)*1000. 
											
										END IF
										
										
										
										!Up to this point, temp_data has played the role of data_in... for 3-hr intervals (strictly!)
										
										temp_data_InTime(:,:,it) = temp_data(:,:)
										
									END DO !DO it = 1, ...
									!variable "it" will be InDimLenRec+1 at this point.
									
									!Locate the position in the file, to write the data
									
									!								counter = 0
									!								DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
									!									
									!									counter = counter + 1
									!									!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
									!									
									!									IF (dtHours.eq.24) THEN !we are averaging per day.. so compare DAYS
									!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
									!										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) .AND. &
									!										( TimeRefArraySelYear(i,4) == InDateTimeDay(1)   ) ) THEN
									!											
									!											EXIT
									!											
									!										END IF
									!									
									!									ELSE IF (dtHours.eq.720) THEN !we are averaging per month... so compare Months
									!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
									!										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) ) THEN
									!											
									!											EXIT
									!											
									!										END IF
									!									ELSE IF (dtHours.eq.2880) THEN !we are averaging over seasons.. so ..................
									!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
									!										( TimeRefArraySelYear(i,3) > InDateTimeMonth(1)  ) ) THEN
									!											
									!											EXIT
									!											
									!										END IF
									!									END IF
									!										
									!								END DO
									!								
									!								PRINT *, "index where in the NC file the WRF data is sorted in = ", &
									!								counter
									!!!!! THIS WAS ALL COMMENTED OUT BECAUSE this is a 3-hr processing (record_count == 1).. no average requested.
									!!!!! The procedure to find where in the output file the data will be written took place earlier! (counter)
									
									
									DO writeIndex = 1,InDimLenRec,1
										!-------------------------------------------------------------------------------
										! write data to NetCDF file
										
										counter = counter_array(writeIndex) !CHANGED
										data_in = temp_data_InTime(:,:,writeIndex)
										!!!!!!! ARIS CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!
										
										print *,'write data to NetCDF file'
										print *,'fn_out',fn_out
										sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid )
										IF (sts/=0) EXIT
										print *, 'NF90_OPEN',  sts
										sts = NF90_INQ_VARID(ncid, TRIM(var_cmip(ivar)), x_varid)
										print *, 'NF90_INQ_VARID', ncid
										print *, 'var_cmip(ivar)', var_cmip(ivar)
										print *, 'x_varid', x_varid
										print *, 'counter', counter
										print *, 'xfocus', xfocus
										print *, 'yfocus', yfocus
										IF ( height(ivar) /= -999 ) THEN
											
											sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
											START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1, 1 /) )
										
										ELSE
											
											sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
											START=(/ 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1 /) )
											
										END IF
										
										print *,'NF90_PUT_VAR', sts
										print *, 'ncid', ncid
										print *, 'x_varid', x_varid
										print *, 'data_in(50:52,50:52)', data_in(50:52,50:52)
										sts = NF90_CLOSE(ncid)
										
										print *, pn_out//"/"//fn_out
										print *, TRIM(var_cmip(ivar)), xfocus, yfocus, counter, ncid, x_varid
										
										!-------------------------------------------------------------------------------
										
									END DO !writeIndex ...
									
									
									
									!-------------------------------------------------------------------------------
									! next WRF file contains different number of output intervals
									
									DEALLOCATE(InVarDataRec)
									
									DEALLOCATE(InDateTimeYear)
									DEALLOCATE(InDateTimeMonth)
									DEALLOCATE(InDateTimeDay)
									DEALLOCATE(InDateTimeHour)
									
									!-------------------------------------------------------------------------------
									
									PRINT *, "------------------------------------------------------------"
									CALL CPU_TIME(cpuTe)
									PRINT '("CPU timing for 1 WRF file (e.g. 1 month worth of data) = ",F6.3," sec")',cpuTe-cpuTs
									
									
									
								END DO !ifl - specific WRF input file, filelist loop
								
								
								
								
								
								
							
							ELSE IF (dtHours == 24) THEN !daily averages
								
								!print *, "counter", counter
								!do for every file....
								DO ifl = 1, SIZE(fl_wrfout)-1, 1 ! operational: loop over complete filelist
									
									print *,' SIZE(fl_wrfout', SIZE(fl_wrfout)
									
									!     DO ifl = 1, SIZE(fl_wrfxtr), 1 ! operational: loop over complete filelist
									!      print *,' SIZE(fl_wrfout', SIZE(fl_wrfxtr)
									
									!DO ifl = 1, 1, 1 ! testing: loop over specific entry in filelist (e.g. just January)
									
									CALL CPU_TIME(cpuTs)
									
									PRINT *, "------------------------------------------------------------"
									
									IF ( filetype(ivar) == "s" ) THEN
										
										iflWRFin = fl_wrfout(ifl)
									
									ELSE IF ( filetype(ivar) == "x" ) THEN
										
										iflWRFin = fl_wrfxtr(ifl)
										
									END IF
									
									PRINT '(100A)', TRIM(iflWRFin)
									
									!-------------------------------------------------------------------------------
									! which timespan is covered by the WRF outputs?
									! assume timespan wrfout = timespan wrfxtrm
									! this determines how many times the tool has to loop over the inputs
									! also check how many years are covered by a single wrfout and wrfxtrm which
									! determines the output file generation
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncid_in)
									sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
									sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
									sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
									NAME = InDimNameRec, LEN = InDimLenRec)
									!print *, "counter", it
									ALLOCATE(InVarDataRec(InDimLenRec))
									
									sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
									sts = NF90_CLOSE(ncid_in)
									
									print *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)
									
									ALLOCATE(InDateTimeYear(InDimLenRec))
									ALLOCATE(InDateTimeMonth(InDimLenRec))
									ALLOCATE(InDateTimeDay(InDimLenRec))
									ALLOCATE(InDateTimeHour(InDimLenRec))
									
									! this is the temporal coverage of the WRF input data
									DO i = 1, SIZE(InVarDataRec), 1
										
										READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(i), InDateTimeMonth(i), InDateTimeDay(i), InDateTimeHour(i)
										
									END DO
									
									
									!-------------------------------------------------------------------------------
									! check whether a file is needed at all
									! check for first date
									! for subsequent dates, just check whether it is the same as the previous one
									! if this is not the case, then check whether file exists...
									! if it does not exist (default case): create
									! pathname is always needed
									
									! loop over the individual timesteps in the WRF files...
									! this may take some time but it is robust and also extremely large arrays
									! might be used
									! highly robust code
									
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ARIS EDIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									! Allocate correspondingly the temp_data_InTime variable
									
									
									!	SELECT CASE (frequency(ifrq))
									!		CASE ('3hr')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,1), STAT = sts) !run through the entire list of 3-hour records and use them all for the calculations.
									!		CASE ('day')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,8), STAT = sts) !only run once... the daily data will be processed once..
									!		CASE DEFAULT
									!			PRINT *, "something happened... "
									!			total_iters = InDimLenRec !make sure bad results in default behavior (i.e. what was happening anyway before...)
									!			!STOP
									!	END SELECT
									
									
									
									IF ( CheckForLeapyear( InDateTimeYear(1), dataIncludesLeapYearDays ) == 366 ) THEN
										
										NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
									
									ELSE
										
										NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
										
									END IF
									
									
									!!!! Here is the case, when 3-hr intervals have been requested... so no alterations are made,
									!!!! except for removing the irrelevant variables, keeping only the cumulative ones..
									!								DO it = 1, InDimLenRec, 1
									it = 1
									PRINT *, "----------------------------------------"
									PRINT *, "working on date in WRF input file = ", TRIM(InVarDataRec(it))
									PRINT *, "---"
									
									!-------------------------------------------------------------------------------
									! generate path and filename
									! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/EUR-44/CRPGL/
									! ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
									! evspsbl_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_CRPGL-WRFARW331_v1_3hr_
									! 1989010100-1989123121
									
									! output does not yet exist
									! monthy check is basically not even possible, but does not do any harm
									! this is a rerstriction for all those who might have a different file
									! structure
									
									!print *, "Let's see: InDateTimeYearPrev", InDateTimeYearPrev
									!print *, "Let's see: InDateTimeYear(1)", InDateTimeYear(it)
									IF ( InDateTimeYearPrev /= InDateTimeYear(it) ) THEN !.AND. &
										!( InDateTimeMonthPrev /= InDateTimeMonth(it) ) ) THEN
										
										InDateTimeYearPrev = InDateTimeYear(it)
										
										PRINT *, "start of processing or new year encountered -> t ref. vec. and filecheck"
										
										!READ( InDateTimeYear(it), '(4A)' ) InDateTimeYearStr
										!READ( InDateTimeMonth(it), '(2A)' ) InDateTimeMonthStr
										WRITE (InDateTimeYearStr,'(I4.4)') InDateTimeYear(it)
										WRITE (EndInDateTimeYearStr,'(I4.4)') InDateTimeYear(it)+1
										
										
										!WRITE (InDateTimeMonthStr,'(I2.2)') InDateTimeMonth(it)
										!PRINT *, InDateTimeYearStr, InDateTimeMonthStr
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										IF ( cell_methods(ivar) == "mean" ) THEN
											
											EndDate=EndInDateTimeYearStr//"010100"
										
										ELSE
											
											EndDate=InDateTimeYearStr//"123121"
											
										END IF
										
										
										pn_out = TRIM(product)                       // "/" // &
										TRIM(CORDEX_domain)                 // "/" // &
										TRIM(institute_id)                  // "/" // &
										TRIM(driving_model_id)              // "/" // &
										TRIM(driving_experiment_name)       // "/" // &
										TRIM(driving_model_ensemble_member) // "/" // &
										TRIM(model_id)                      // "/" // &
										TRIM(rcm_version_id)                // "/" // &
										TRIM(frequency(ifrq))               // "/" // &
										TRIM(var_cmip(ivar))
										
										fn_out = TRIM(var_cmip(ivar))                // "_" // &
										TRIM(CORDEX_domain)                 // "_" // &
										TRIM(driving_model_id)              // "_" // &
										TRIM(driving_experiment_name)       // "_" // &
										TRIM(driving_model_ensemble_member) // "_" // &
										TRIM(model_id)                      // "_" // &
										TRIM(rcm_version_id)                // "_" // &
										TRIM(frequency(ifrq))               // "_" // &
										!                   InDateTimeYearStr//"010100-"//InDateTimeYearStr//"123121" // &
										InDateTimeYearStr//"010100-"//TRIM(EndDate)//  &
										".nc"
										
										PRINT *, "pn_out = ", TRIM(pn_out)
										PRINT *, "fn_out = ", TRIM(fn_out)
										
										!-------------------------------------------------------------------------------
										! extract the time info from the ref array which fits the respective year
										! ...as there is no "WHERE" the way I need it in F95, use loops
										! this is needed whenever a new netcdf file is to be used, also if
										! this file exists already
										!        READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										DEALLOCATE( TimeRefArraySelYear )
										
										print *, "DEALLOCATE( TimeRefArraySelYear )" 
										
										DEALLOCATE( Time_bnds ) !SKn 
										
										print *, "DEALLOCATE( Time_bnds )"
										
										counter = 0
										PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
										PRINT *,'SHAPE(TimeRefArray, 1)',  SHAPE(TimeRefArray, 1)
										PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
										PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
										PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												
											END IF
											
										END DO
										
										PRINT *, "Counter: ", counter
										! holds data of exactly 1 year
										PRINT *, "timesteps in the time ref. subset = ", counter
										ALLOCATE( TimeRefArraySelYear( counter, 5 ) ) ! index, y, m, d, h
										
										print *, "ALLOCATE( TimeRefArraySelYear( counter, 5 ) )"
										
										print *, "TimeRefArraySelYear( counter, 1 )", TimeRefArraySelYear( counter, 1 )
										print *, "TimeRefArraySelYear( counter, 2 )", TimeRefArraySelYear( counter, 2 )
										print *, "TimeRefArraySelYear( counter, 3 )", TimeRefArraySelYear( counter, 3 )
										print *, "TimeRefArraySelYear( counter, 4 )", TimeRefArraySelYear( counter, 4 )
										print *, "TimeRefArraySelYear( counter, 5 )", TimeRefArraySelYear( counter, 5 )
										
										! find the matching elements of the respecitve year and copy them
										
										ALLOCATE( Time_bnds( 2, counter ) )
										
										print *, "ALLOCATE( Time_bnds( 2, counter ) )", Time_bnds( 2, counter )
										
										counter = 0
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												TimeRefArraySelYear(counter,1:5) = TimeRefArray(i,1:5)
												
											END IF
											
										END DO
										
										!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySelYear(:,:) )
										
										!-------------------------------------------------------------------------------
										! check for existance of the file and generate file if needed
										
										INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )
										
										! could exist already from a previous run of tool and due to multiple
										! months in a file (i.e. 1 WRF output file may cover several months)
										IF ( FileExists ) THEN
											
											PRINT *, "path and file exist, continue filling"
										
										ELSE
											
											PRINT *, "path and file do not yet exist, create path and NetCDF file first"
											PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)
											
											CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )
											
											!-------------------------------------------------------------------------------
											! non-std, works for gfortran (fct & subroutine) + ifort
											! comment lines in the NetCDF file global attribute definition
											! turn standard checking in Makefile off
											! trackingID = "xxxxxxxx-xxxx-Mxxx-Nxxx-xxxxxxxxxxxx"
											! creationDate = "YYYY-MM-DD-THH:MM:SSZ"
											
											CALL SYSTEM(cmdUUID)
											OPEN(1,FILE="tmpfileUUID",STATUS='old')
											READ(1,*) trackingID
											CLOSE(1)
											PRINT *, "uuidgen externally generated trackingID = ", trackingID
											
											CALL SYSTEM(cmdDate)
											OPEN(1,FILE="tmpfileDate",STATUS='old')
											READ(1,*) creationDate
											CLOSE(1)
											PRINT *, "date externally generated creation date = ", creationDate
											
											!-------------------------------------------------------------------------------
											! create NetCDF file
											! NF90_CLASSIC_MODEL = NetCDF4_classic
											! NF90_HDF5 = NetCDF4 based on HDF5
											! NF90_CLOBBER = old NetCDF
											! sts = NF90_CREATE(PathFileNameOutTEST, NF90_HDF5, ncid)
											
											!comb_flags = IOR(NF90_HDF5, NF90_CLASSIC_MODEL)
											!https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/NF90_005fCREATE.html
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid)   !not sure whether this is the right data format spec. I guess it may be right using compression but not the other fancy stuff from NetCDF4
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncid)   !if anything, then use this here
											sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_HDF5, ncid)
											
											
											! always included
											sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
											sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid)
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
												
												ELSE
													
													sts = NF90_DEF_DIM(ncid, "plev", 1, height_dimid)
													
												END IF
												
											END IF
											sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
												
											END IF
											
											
											! always included
											sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lon_varid)
											sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "long_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
											
											! always included
											sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lat_varid)
											sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "long_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
											
											! always included
											sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid)
											sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
											sts = nf90_put_att(ncid, rlon_varid, "long_name", "longitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
											
											! always included
											sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid)
											sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
											sts = nf90_put_att(ncid, rlat_varid, "long_name", "latitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")
											
											! always included
											! restriction to one domain only
											sts = nf90_def_var(ncid, "rotated_pole", NF90_DOUBLE, rotated_pole_varid)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", 39.25)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", -162.0)
											
											! depends whether height is set in the nml
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_def_var(ncid, "height", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
													sts = nf90_put_att(ncid, height_varid, "long_name", "height")
													sts = nf90_put_att(ncid, height_varid, "units", "m")
													sts = nf90_put_att(ncid, height_varid, "positive", "up")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")
												
												ELSE
													
													sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name","air_pressure")
													sts = nf90_put_att(ncid, height_varid, "long_name", "pressure")
													sts = nf90_put_att(ncid, height_varid, "units", "Pa")
													sts = nf90_put_att(ncid, height_varid, "positive", "down")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")                 
													
												END IF
												
											END IF
											
											!missing: lvl, depends
											
											! always included
											sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid)
											sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "long_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "units", "days since 1949-12-01T00:00:00Z")
											sts = nf90_put_att(ncid, rec_varid, "calendar", "standard")
											sts = nf90_put_att(ncid, rec_varid, "axis", "T")
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												! Add time bounds to time varaible
												
												sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")	
												sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid)
												sts = nf90_put_att(ncid, recbnds_varid, "standard_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "long_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "units", "days since 1949-12-01T00:00:00Z")
												sts = nf90_put_att(ncid, recbnds_varid, "calendar", "standard")
												!sts = nf90_put_att(ncid, recbnds_varid, "axis", "T")
												
												print *,'rec_varid', rec_varid
												print *,'recbnds_varid', recbnds_varid
												
											END IF
											
											! always included
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact", contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "creation_date", creationDate)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment", experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment_id", experiment_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment", driving_experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_id", driving_model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_ensemble_member", driving_model_ensemble_member)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment_name", driving_experiment_name)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "frequency", frequency(ifrq))
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution", institution)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_id", institute_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "model_id", model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "forcing",forcing)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact",contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source",source)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "rcm_version_id", rcm_version_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "CORDEX_domain", CORDEX_domain)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id", trackingID)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "title",title)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_run_id", institute_run_id)
											
											! always included
											IF ( height(ivar) /= -999 ) THEN
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, height_dimid, rec_dimid /), x_varid)
											
											ELSE
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid)
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
											
											IF ( positive(ivar) /= '-999' ) THEN
												
												sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
											
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat height")
												
												ELSE
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat plev")
													
												END IF
											
											ELSE
												
												sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat")
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "grid_mapping", "rotated_pole")
											
											!IF ( height(ivar) == 850 ) THEN
											sts = nf90_put_att(ncid, x_varid, "missing_value", 1.e20)
											sts = nf90_put_att(ncid, x_varid, "_FillValue", 1.e20)
											!END IF
											
											sts = NF90_ENDDEF(ncid)
											
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF
											!print *,'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
											
											
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											
											! add time_bnds, calc here
											
											! add height from NML
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
												
												ELSE
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar)*100. )
													
												END IF
												
											END IF
											
											sts = NF90_CLOSE(ncid)
											
										END IF ! file exists y/n
										
									END IF ! checking with previous year and month
									
									!-------------------------------------------------------------------------------
									! match timestep of WRFin with the subset of the ref time vec which belong to
									! the NetCDF file of the year currently open to receive data
									! NOT SURE WHETHER THIS IS NEEDED AT ALL, THIS IS THE WRONG DIRECTION ???
									! see whether the current time of the timestep fits anywhere in the
									
									
									!-------------------------------------------------------------------------------
									! read orig WRF outpoutses
									! there is always a corresponding time-slot in the NC file
									! extracted time from above
									! "it" controls it all: timestep in the individual WRF file
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)
									
									if (ifl /= SIZE(fl_wrfout,1)) THEN
										sts = NF90_OPEN(fl_wrfout(ifl+1), NF90_NOWRITE, ncidin0)
									END IF
									
									
									!!!!!!!!!!!!! ONLY READ THE VARIABLES THAT PROVIDE CUMULATIVE WRF OUTPUTS.......!!!
									IF (var_cmip(ivar) == "pr") THEN 
										
										ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										!get first day
										sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus, 1 /) )   !read two timesteps to calculate 3hr sum
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin0, rainnc_varid, rainnc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )   !read last timestep of previous wrfout file
										
										sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prc") THEN
										
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prsn") THEN
										
										ALLOCATE( snownc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin0, snownc_varid, snownc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "snm") THEN
										
										ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin0, acsnom_varid, acsnom_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "evspsbl") THEN
										
										ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfcevp_varid, sfcevp_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "evspsblpot") THEN
										
										ALLOCATE( potevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin0, potevp_varid, potevp_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrros") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrro") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin0, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin0, udroff_varid, udroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF ((var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")  &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "rlut")                                 &
									.or. (var_cmip(ivar) == "rsdt") .or. (var_cmip(ivar) == "rsut") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										ALLOCATE( rad_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin0, varid, rad_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										
										print*, var_cmip(ivar), rad_in(50,50,1), rad_in(50,50,2)
										print*, 'difference in J m-2', (rad_in(50,50,2) - rad_in(50,50,1))
										print*, 'in mean W m-2', (rad_in(50,50,2) - rad_in(50,50,1))/ (obs_interval*3600.)
										
										! alternative: since accumulated values as read above get so large in 
										! long term simulations that their differences loose accuracy, use 
										! instantaneous values instead and calculate means
										
									
									ELSE IF (var_cmip(ivar) == "mrso") THEN
										
										ALLOCATE( smois_in( xfocus, yfocus, 4, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,1), &
										START = (/ xoffset, yoffset, 1, 1 /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin0, varid, smois_in(:,:,:,2), &
										START = (/ xoffset, yoffset, 1, 1 /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
										
										
										!!!!!!!!!!!!!!! THIS ELSE IS NEVER GOING TO BE EXECUTED... IT IS JUST COMMENTED OUT
										!									ELSE 
										!										
										!										ALLOCATE( temp_data ( xfocus, yfocus ), STAT=sts )
										!										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										!										!print *, "VALUE FOR NAME: ", TRIM(var_wrf(ivar)), "+ VARID: ", varid
										!										
										!										
										!										sts = NF90_GET_VAR(ncidin, varid, temp_data(:,:), &
										!										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										!										!			data_in(:,:) = SUM(temp_data, 3)/8.
										!										!			print *, temp_data(50,50,1), temp_data(50,50,2), temp_data(50,50,3), temp_data(50,50,4)
										!										!			print *, temp_data(50,50,5), temp_data(50,50,6), temp_data(50,50,7), temp_data(50,50,8)
										!										!			print *, "Mean of data_in(50,50):", data_in(50,50)
									END IF
									
									sts = NF90_CLOSE(ncidin)
									sts = NF90_CLOSE(ncidin0) !close both files
									
									
									!-------------------------------------------------------------------------------
									! some analysis of the data
									
									print *,"shape of array" , SHAPE(temp_data)
									print *,"size of array" , SIZE(temp_data)
									!
									!       stat_mean = SUM(temp_data(:,:,5))/(MAX(1,SIZE(temp_data(:,:,5))))
									!       PRINT *, stat_mean
									stat_mean = SUM(temp_data(:,:))/SIZE(temp_data(:,:))
									PRINT *, stat_mean
									
									!-------------------------------------------------------------------------------
									! processing
									
									
									!       ***pr***
									IF (var_cmip(ivar) == "pr") THEN 
										
										temp_data(:,:) = ((rainnc_in(:,:,2) + rainc_in(:,:,2)) - (rainnc_in(:,:,1) + rainc_in(:,:,1)))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										!ATTENTION: implement adjustable time intervals that the differences are devided by
										
									END IF
									
									
									!       ***prc***
									IF (var_cmip(ivar) == "prc") THEN
										
										temp_data(:,:) = (rainc_in(:,:,2) - rainc_in(:,:,1))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***prsn***
									IF (var_cmip(ivar) == "prsn") THEN
										
										temp_data(:,:) = (snownc_in(:,:,2) - snownc_in(:,:,1))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***snm***
									IF (var_cmip(ivar) == "snm") THEN
										
										temp_data(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1))/(dtHours*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspsbl***
									IF (var_cmip(ivar) == "evspsbl") THEN
										
										temp_data(:,:) = (sfcevp_in(:,:,2) - sfcevp_in(:,:,1))/(dtHours*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspblpot**
									IF (var_cmip(ivar) == "evspsblpot") THEN
										
										temp_data(:,:) = (potevp_in(:,:,2) - (potevp_in(:,:,1)))/L   !unit [W m-2]/[J kg-1] -> [kg m-2 s-1]
										
										! THERE IS STH WRONG WITH THE UNITS: WRF's POTEVP is accumulated and declared to be in W m-2. 
										! It doesen't make sense to accumulate in W m-2, but even if assume it as J m-2 or derive kg m-2 
										! by using latent heat of vaporization you never get values that have a reasonable magnitude...
										
										
									END IF
									
									!       ***mrros***
									IF (var_cmip(ivar) == "mrros") THEN
										
										temp_data(:,:) = (sfroff_in(:,:,2) - sfroff_in(:,:,1))/(dtHours*3600.)       !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***mrro***
									IF (var_cmip(ivar) == "mrro") THEN
										
										temp_data(:,:) = ((sfroff_in(:,:,2) - sfroff_in(:,:,1)) + (udroff_in(:,:,2) - udroff_in(:,:,1)))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!!!!!!! This if has probably forgotten to include 3 variables (rlut, rsdt, rsut).
									!!!!!! They are read from the files earlier altogether, but have not been included here.
									!!!!!!!! Probably a mistake... add them some time...
									!       ***rsds, rlds, rsus, rlus***
									IF ( (var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")      &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) - rad_in(:,:,1)) /(dtHours*3600.)       ! take difference of accumulated values
										
										ELSE IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation_alternative") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) + rad_in(:,:,1)) / 2.              ! take mean of instantaneous values
											
										END IF 
										
									END IF
									
									
									!       ***mrso***
									IF (var_cmip(ivar) == "mrso") THEN
										
										temp_data(:,:) = (((smois_in(:,:,1,1)*0.1 + smois_in(:,:,2,1)*0.3 + smois_in(:,:,3,1)*0.6 + smois_in(:,:,4,1)*1.0 ) + &
										(smois_in(:,:,1,2)*0.1 + smois_in(:,:,2,2)*0.3 + smois_in(:,:,3,2)*0.6 + smois_in(:,:,4,2)*1.0 ))/2.)*1000. 
										
									END IF
									
									
									
									!Up to this point, temp_data has played the role of data_in... for 3-hr intervals (strictly!)
									
									!									temp_data_InTime(:,:,it) = temp_data(:,:)
									
									!								END DO !DO it = 1, ...
									!variable "it" will be InDimLenRec+1 at this point.
									
									!Locate the position in the file, to write the data

									counter = 0
									DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
										!									
										counter = counter + 1

										!									!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
										!									
										!									IF (dtHours.eq.24) THEN !we are averaging per day.. so compare DAYS
										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) .AND. &
										( TimeRefArraySelYear(i,4) == InDateTimeDay(1)   ) ) THEN

											EXIT
											
										END IF
										!									
										!									ELSE IF (dtHours.eq.720) THEN !we are averaging per month... so compare Months
										!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										!										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) ) THEN
										!											
										!											EXIT
										!											
										!										END IF
										!									ELSE IF (dtHours.eq.2880) THEN !we are averaging over seasons.. so ..................
										!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										!										( TimeRefArraySelYear(i,3) > InDateTimeMonth(1)  ) ) THEN
										!											
										!											EXIT
										!											
										!										END IF
										!									END IF
										!										
									END DO
									!								
									PRINT *, "index where in the NC file the WRF data is sorted in = ", &
									counter
									
									!now the writing process takes place only once!
									!									DO writeIndex = 1,InDimLenRec,1
									!-------------------------------------------------------------------------------
									! write data to NetCDF file
									
									data_in = temp_data 
									!!!!!!! ARIS CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!
									
									print *,'write data to NetCDF file'
									print *,'fn_out',fn_out
									sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid )
									IF (sts/=0) EXIT
									print *, 'NF90_OPEN',  sts
									sts = NF90_INQ_VARID(ncid, TRIM(var_cmip(ivar)), x_varid)
									print *, 'NF90_INQ_VARID', ncid
									print *, 'var_cmip(ivar)', var_cmip(ivar)
									print *, 'x_varid', x_varid
									print *, 'counter', counter
									print *, 'xfocus', xfocus
									print *, 'yfocus', yfocus
									IF ( height(ivar) /= -999 ) THEN
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1, 1 /) )
									
									ELSE
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									END IF
									
									print *,'NF90_PUT_VAR', sts
									print *, 'ncid', ncid
									print *, 'x_varid', x_varid
									print *, 'data_in(50:52,50:52)', data_in(50:52,50:52)
									sts = NF90_CLOSE(ncid)
									
									print *, pn_out//"/"//fn_out
									print *, TRIM(var_cmip(ivar)), xfocus, yfocus, counter, ncid, x_varid
									
									!-------------------------------------------------------------------------------
									
									!									END DO !writeIndex ...
									!single write process... no iteration here..									
									
									
									!-------------------------------------------------------------------------------
									! next WRF file contains different number of output intervals
									print *, "deallocation...."
									DEALLOCATE(InVarDataRec)
									print *, "in process....."
									
									DEALLOCATE(InDateTimeYear)
									DEALLOCATE(InDateTimeMonth)
									DEALLOCATE(InDateTimeDay)
									DEALLOCATE(InDateTimeHour)
									
									!-------------------------------------------------------------------------------
									
									PRINT *, "------------------------------------------------------------"
									CALL CPU_TIME(cpuTe)
									PRINT '("CPU timing for 1 WRF file (e.g. 1 month worth of data) = ",F6.3," sec")',cpuTe-cpuTs
									
									
									
								END DO !ifl - specific WRF input file, filelist loop
								
								
								
								
							
							ELSE IF (dtHours == 720) THEN !monthly averages
								
								!do for every OTHER file....
								ifl = 1
								DO WHILE (ifl < (SIZE(fl_wrfout) - 31))
									!							DO ifl = 1, SIZE(fl_wrfout)-1, 1 ! operational: loop over complete filelist
									
									print *,' SIZE(fl_wrfout', SIZE(fl_wrfout)
									
									!     DO ifl = 1, SIZE(fl_wrfxtr), 1 ! operational: loop over complete filelist
									!      print *,' SIZE(fl_wrfout', SIZE(fl_wrfxtr)
									
									!DO ifl = 1, 1, 1 ! testing: loop over specific entry in filelist (e.g. just January)
									
									CALL CPU_TIME(cpuTs)
									
									PRINT *, "------------------------------------------------------------"
									
									IF ( filetype(ivar) == "s" ) THEN
										
										iflWRFin = fl_wrfout(ifl)
									
									ELSE IF ( filetype(ivar) == "x" ) THEN
										
										iflWRFin = fl_wrfxtr(ifl)
										
									END IF
									
									PRINT '(100A)', TRIM(iflWRFin)
									
									!-------------------------------------------------------------------------------
									! which timespan is covered by the WRF outputs?
									! assume timespan wrfout = timespan wrfxtrm
									! this determines how many times the tool has to loop over the inputs
									! also check how many years are covered by a single wrfout and wrfxtrm which
									! determines the output file generation
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncid_in)
									sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
									sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
									sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
									NAME = InDimNameRec, LEN = InDimLenRec)
									ALLOCATE(InVarDataRec(InDimLenRec))
									sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
									sts = NF90_CLOSE(ncid_in)
									
									print *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)
									
									ALLOCATE(InDateTimeYear(InDimLenRec))
									ALLOCATE(InDateTimeMonth(InDimLenRec))
									ALLOCATE(InDateTimeDay(InDimLenRec))
									ALLOCATE(InDateTimeHour(InDimLenRec))
									
									! this is the temporal coverage of the WRF input data
									DO i = 1, SIZE(InVarDataRec), 1
										
										READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(i), InDateTimeMonth(i), InDateTimeDay(i), InDateTimeHour(i)
										
									END DO
									
									
									!-------------------------------------------------------------------------------
									! check whether a file is needed at all
									! check for first date
									! for subsequent dates, just check whether it is the same as the previous one
									! if this is not the case, then check whether file exists...
									! if it does not exist (default case): create
									! pathname is always needed
									
									! loop over the individual timesteps in the WRF files...
									! this may take some time but it is robust and also extremely large arrays
									! might be used
									! highly robust code
									
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ARIS EDIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									! Allocate correspondingly the temp_data_InTime variable
									
									
									!	SELECT CASE (frequency(ifrq))
									!		CASE ('3hr')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,1), STAT = sts) !run through the entire list of 3-hour records and use them all for the calculations.
									!		CASE ('day')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,8), STAT = sts) !only run once... the daily data will be processed once..
									!		CASE DEFAULT
									!			PRINT *, "something happened... "
									!			total_iters = InDimLenRec !make sure bad results in default behavior (i.e. what was happening anyway before...)
									!			!STOP
									!	END SELECT
									
									
									
									IF ( CheckForLeapyear( InDateTimeYear(1), dataIncludesLeapYearDays ) == 366 ) THEN
										
										NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
									
									ELSE
										
										NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
										
									END IF
									
									
									!!!! Here is the case, when 3-hr intervals have been requested... so no alterations are made,
									!!!! except for removing the irrelevant variables, keeping only the cumulative ones..
									!								DO it = 1, InDimLenRec, 1
									it = 1
									PRINT *, "----------------------------------------"
									PRINT *, "working on date in WRF input file = ", TRIM(InVarDataRec(it))
									PRINT *, "---"
									
									!-------------------------------------------------------------------------------
									! generate path and filename
									! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/EUR-44/CRPGL/
									! ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
									! evspsbl_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_CRPGL-WRFARW331_v1_3hr_
									! 1989010100-1989123121
									
									! output does not yet exist
									! monthy check is basically not even possible, but does not do any harm
									! this is a rerstriction for all those who might have a different file
									! structure
									IF ( ( InDateTimeYearPrev /= InDateTimeYear(it) ) .AND. ( InDateTimeYear(it) /= LastYearOfRun + 1 ) ) THEN !.AND. &
										!( InDateTimeMonthPrev /= InDateTimeMonth(it) ) ) THEN
										
										InDateTimeYearPrev = InDateTimeYear(it)
										counter_prev_year = counter
										
										PRINT *, "start of processing or new year encountered -> t ref. vec. and filecheck"
										
										!READ( InDateTimeYear(it), '(4A)' ) InDateTimeYearStr
										!READ( InDateTimeMonth(it), '(2A)' ) InDateTimeMonthStr
										
										!Change naming rules for producing just one output file when running for multiple years
										!WRITE (EndInDateTimeYearStr,'(I4.4)') InDateTimeYear(it)+1
										WRITE (InDateTimeYearStr,'(I4)') FirstYearOfRun
										!WRITE (InDateTimeYearStr,'(I4.4)') InDateTimeYear(it)		
										WRITE (EndInDateTimeYearStr,'(I4)') LastYearOfRun
										
									
										!WRITE (InDateTimeMonthStr,'(I2.2)') InDateTimeMonth(it)
										!PRINT *, InDateTimeYearStr, InDateTimeMonthStr
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										IF ( cell_methods(ivar) == "mean" ) THEN
											
											WRITE (EndInDateTimeYearStr,'(I4)') (LastYearOfRun + 1)
											EndDate=EndInDateTimeYearStr//"010100"
											
										
										ELSE
											
											EndDate=EndInDateTimeYearStr//"123121"
											
										END IF
										
										
										pn_out = TRIM(product)                       // "/" // &
										TRIM(CORDEX_domain)                 // "/" // &
										TRIM(institute_id)                  // "/" // &
										TRIM(driving_model_id)              // "/" // &
										TRIM(driving_experiment_name)       // "/" // &
										TRIM(driving_model_ensemble_member) // "/" // &
										TRIM(model_id)                      // "/" // &
										TRIM(rcm_version_id)                // "/" // &
										TRIM(frequency(ifrq))               // "/" // &
										TRIM(var_cmip(ivar))
										
										fn_out = TRIM(var_cmip(ivar))                // "_" // &
										TRIM(CORDEX_domain)                 // "_" // &
										TRIM(driving_model_id)              // "_" // &
										TRIM(driving_experiment_name)       // "_" // &
										TRIM(driving_model_ensemble_member) // "_" // &
										TRIM(model_id)                      // "_" // &
										TRIM(rcm_version_id)                // "_" // &
										TRIM(frequency(ifrq))               // "_" // &
										!                   InDateTimeYearStr//"010100-"//InDateTimeYearStr//"123121" // &
										InDateTimeYearStr//"010100-"//TRIM(EndDate)//  &
										".nc"
										
										PRINT *, "pn_out = ", TRIM(pn_out)
										PRINT *, "fn_out = ", TRIM(fn_out)
										
										!-------------------------------------------------------------------------------
										! extract the time info from the ref array which fits the respective year
										! ...as there is no "WHERE" the way I need it in F95, use loops
										! this is needed whenever a new netcdf file is to be used, also if
										! this file exists already
										!        READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										DEALLOCATE( TimeRefArraySelYear )
										
										print *, "DEALLOCATE( TimeRefArraySelYear )" 
										
										DEALLOCATE( Time_bnds ) !SKn 
										
										print *, "DEALLOCATE( Time_bnds )"
										
										counter = 0
										PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
										PRINT *,'SHAPE(TimeRefArray, 1)',  SHAPE(TimeRefArray, 1)
										PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
										PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
										PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												
											END IF
											
										END DO
										PRINT *, "Counter: ", counter
										! holds data of exactly 1 year
										PRINT *, "timesteps in the time ref. subset = ", counter
										ALLOCATE( TimeRefArraySelYear( counter, 5 ) ) ! index, y, m, d, h
										
										print *, "ALLOCATE( TimeRefArraySelYear( counter, 5 ) )"
										
										print *, "TimeRefArraySelYear( counter, 1 )", TimeRefArraySelYear( counter, 1 )
										print *, "TimeRefArraySelYear( counter, 2 )", TimeRefArraySelYear( counter, 2 )
										print *, "TimeRefArraySelYear( counter, 3 )", TimeRefArraySelYear( counter, 3 )
										print *, "TimeRefArraySelYear( counter, 4 )", TimeRefArraySelYear( counter, 4 )
										print *, "TimeRefArraySelYear( counter, 5 )", TimeRefArraySelYear( counter, 5 )
										
										! find the matching elements of the respecitve year and copy them
										
										ALLOCATE( Time_bnds( 2, counter ) )
										
										print *, "ALLOCATE( Time_bnds( 2, counter ) )", Time_bnds( 2, counter )
										
										counter = 0
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												TimeRefArraySelYear(counter,1:5) = TimeRefArray(i,1:5)
												
											END IF
											
										END DO
										
										!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySelYear(:,:) )
										
										!-------------------------------------------------------------------------------
										! check for existance of the file and generate file if needed
										
										INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )
										
										! could exist already from a previous run of tool and due to multiple
										! months in a file (i.e. 1 WRF output file may cover several months)
										IF ( FileExists ) THEN
											
											PRINT *, "path and file exist, continue filling"
											
											!!!!Change indices for TimeRefArraySelYear (time) and Time_bnds (time bounds) records
											!!!!so that they are written in proper position in same file right after the last year's data
											sts = NF90_OPEN(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid)
											
											!Not sure if all of this is necessary though...
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.), START = (/ counter_prev_year + 1, 1 /) , COUNT = (/ SIZE(TimeRefArraySelYear(:,1), 1) /) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, counter_prev_year+1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF											
											
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											
										
										ELSE
											
											PRINT *, "path and file do not yet exist, create path and NetCDF file first"
											PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)
											
											CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )
											
											!-------------------------------------------------------------------------------
											! non-std, works for gfortran (fct & subroutine) + ifort
											! comment lines in the NetCDF file global attribute definition
											! turn standard checking in Makefile off
											! trackingID = "xxxxxxxx-xxxx-Mxxx-Nxxx-xxxxxxxxxxxx"
											! creationDate = "YYYY-MM-DD-THH:MM:SSZ"
											
											CALL SYSTEM(cmdUUID)
											OPEN(1,FILE="tmpfileUUID",STATUS='old')
											READ(1,*) trackingID
											CLOSE(1)
											PRINT *, "uuidgen externally generated trackingID = ", trackingID
											
											CALL SYSTEM(cmdDate)
											OPEN(1,FILE="tmpfileDate",STATUS='old')
											READ(1,*) creationDate
											CLOSE(1)
											PRINT *, "date externally generated creation date = ", creationDate
											
											!-------------------------------------------------------------------------------
											! create NetCDF file
											! NF90_CLASSIC_MODEL = NetCDF4_classic
											! NF90_HDF5 = NetCDF4 based on HDF5
											! NF90_CLOBBER = old NetCDF
											! sts = NF90_CREATE(PathFileNameOutTEST, NF90_HDF5, ncid)
											
											!comb_flags = IOR(NF90_HDF5, NF90_CLASSIC_MODEL)
											!https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/NF90_005fCREATE.html
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid)   !not sure whether this is the right data format spec. I guess it may be right using compression but not the other fancy stuff from NetCDF4
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncid)   !if anything, then use this here
											sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_HDF5, ncid)
											
											
											! always included
											sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
											sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid)
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
												
												ELSE
													
													sts = NF90_DEF_DIM(ncid, "plev", 1, height_dimid)
													
												END IF
												
											END IF
											sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
												
											END IF
											
											
											! always included
											sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lon_varid)
											sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "long_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
											
											! always included
											sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lat_varid)
											sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "long_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
											
											! always included
											sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid)
											sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
											sts = nf90_put_att(ncid, rlon_varid, "long_name", "longitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
											
											! always included
											sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid)
											sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
											sts = nf90_put_att(ncid, rlat_varid, "long_name", "latitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")
											
											! always included
											! restriction to one domain only
											sts = nf90_def_var(ncid, "rotated_pole", NF90_DOUBLE, rotated_pole_varid)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", 39.25)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", -162.0)
											
											! depends whether height is set in the nml
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_def_var(ncid, "height", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
													sts = nf90_put_att(ncid, height_varid, "long_name", "height")
													sts = nf90_put_att(ncid, height_varid, "units", "m")
													sts = nf90_put_att(ncid, height_varid, "positive", "up")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")
												
												ELSE
													
													sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name","air_pressure")
													sts = nf90_put_att(ncid, height_varid, "long_name", "pressure")
													sts = nf90_put_att(ncid, height_varid, "units", "Pa")
													sts = nf90_put_att(ncid, height_varid, "positive", "down")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")                 
													
												END IF
												
											END IF
											
											!missing: lvl, depends
											
											! always included
											sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid)
											sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "long_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "units", "days since 1949-12-01T00:00:00Z")
											sts = nf90_put_att(ncid, rec_varid, "calendar", "standard")
											sts = nf90_put_att(ncid, rec_varid, "axis", "T")
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												! Add time bounds to time varaible
												
												sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")	
												sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid)
												sts = nf90_put_att(ncid, recbnds_varid, "standard_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "long_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "units", "days since 1949-12-01T00:00:00Z")
												sts = nf90_put_att(ncid, recbnds_varid, "calendar", "standard")
												!sts = nf90_put_att(ncid, recbnds_varid, "axis", "T")
												
												print *,'rec_varid', rec_varid
												print *,'recbnds_varid', recbnds_varid
												
											END IF
											
											! always included
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact", contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "creation_date", creationDate)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment", experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment_id", experiment_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment", driving_experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_id", driving_model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_ensemble_member", driving_model_ensemble_member)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment_name", driving_experiment_name)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "frequency", frequency(ifrq))
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution", institution)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_id", institute_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "model_id", model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "forcing",forcing)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact",contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source",source)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "rcm_version_id", rcm_version_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "CORDEX_domain", CORDEX_domain)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id", trackingID)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "title",title)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_run_id", institute_run_id)
											
											! always included
											IF ( height(ivar) /= -999 ) THEN
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, height_dimid, rec_dimid /), x_varid)
											
											ELSE
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid)
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
											
											IF ( positive(ivar) /= '-999' ) THEN
												
												sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
											
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat height")
												
												ELSE
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat plev")
													
												END IF
											
											ELSE
												
												sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat")
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "grid_mapping", "rotated_pole")
											
											!IF ( height(ivar) == 850 ) THEN
											sts = nf90_put_att(ncid, x_varid, "missing_value", 1.e20)
											sts = nf90_put_att(ncid, x_varid, "_FillValue", 1.e20)
											!END IF
											
											sts = NF90_ENDDEF(ncid)
											
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF
											!print *,'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
											
											
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											
											! add time_bnds, calc here
											
											! add height from NML
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
												
												ELSE
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar)*100. )
													
												END IF
												
											END IF
											
											sts = NF90_CLOSE(ncid)
											
										END IF ! file exists y/n
										
									END IF ! checking with previous year and month
									
									!-------------------------------------------------------------------------------
									! match timestep of WRFin with the subset of the ref time vec which belong to
									! the NetCDF file of the year currently open to receive data
									! NOT SURE WHETHER THIS IS NEEDED AT ALL, THIS IS THE WRONG DIRECTION ???
									! see whether the current time of the timestep fits anywhere in the
									
									
									!-------------------------------------------------------------------------------
									! read orig WRF outpoutses
									! there is always a corresponding time-slot in the NC file
									! extracted time from above
									! "it" controls it all: timestep in the individual WRF file
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)
									
									!< and NOT <= is used, because we want to avoid the last december..
									
									if (ifl < (SIZE(fl_wrfout,1)-31) ) THEN
										sts = NF90_OPEN(fl_wrfout(ifl+NoOfDaysPerMonth(InDateTimeMonth(1))), NF90_NOWRITE, ncidin0)
									END IF
									!InDateTimeMonth is the month in which the first 3-hr interval refers to.
									!All intervals refer to the same month anyway, so the first is used.
									!According to whether it is a leapyear or not, the appropriate number of days are addded...
									!in order to refer to the first day of the next month..
									
									!!!!!!!!!!!!! ONLY READ THE VARIABLES THAT PROVIDE CUMULATIVE WRF OUTPUTS.......!!!
									IF (var_cmip(ivar) == "pr") THEN 
										
										ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										!get first day
										sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus, 1 /) )   !read two timesteps to calculate 3hr sum
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin0, rainnc_varid, rainnc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )   !read last timestep of previous wrfout file
										
										sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prc") THEN
										
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prsn") THEN
										
										ALLOCATE( snownc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin0, snownc_varid, snownc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "snm") THEN
										
										ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin0, acsnom_varid, acsnom_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "evspsbl") THEN
										
										ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfcevp_varid, sfcevp_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "evspsblpot") THEN
										
										ALLOCATE( potevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin0, potevp_varid, potevp_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrros") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrro") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin0, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin0, udroff_varid, udroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF ((var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")  &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "rlut")                                 &
									.or. (var_cmip(ivar) == "rsdt") .or. (var_cmip(ivar) == "rsut") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										ALLOCATE( rad_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin0, varid, rad_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										
										print*, var_cmip(ivar), rad_in(50,50,1), rad_in(50,50,2)
										print*, 'difference in J m-2', (rad_in(50,50,2) - rad_in(50,50,1))
										print*, 'in mean W m-2', (rad_in(50,50,2) - rad_in(50,50,1))/ (obs_interval*3600.)
										
										! alternative: since accumulated values as read above get so large in 
										! long term simulations that their differences loose accuracy, use 
										! instantaneous values instead and calculate means
										
									
									ELSE IF (var_cmip(ivar) == "mrso") THEN
										
										ALLOCATE( smois_in( xfocus, yfocus, 4, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,1), &
										START = (/ xoffset, yoffset, 1, 1 /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin0, varid, smois_in(:,:,:,2), &
										START = (/ xoffset, yoffset, 1, 1 /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
										
										
										!!!!!!!!!!!!!!! THIS ELSE IS NEVER GOING TO BE EXECUTED... IT IS JUST COMMENTED OUT
										!									ELSE 
										!										
										!										ALLOCATE( temp_data ( xfocus, yfocus ), STAT=sts )
										!										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										!										!print *, "VALUE FOR NAME: ", TRIM(var_wrf(ivar)), "+ VARID: ", varid
										!										
										!										
										!										sts = NF90_GET_VAR(ncidin, varid, temp_data(:,:), &
										!										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										!										!			data_in(:,:) = SUM(temp_data, 3)/8.
										!										!			print *, temp_data(50,50,1), temp_data(50,50,2), temp_data(50,50,3), temp_data(50,50,4)
										!										!			print *, temp_data(50,50,5), temp_data(50,50,6), temp_data(50,50,7), temp_data(50,50,8)
										!										!			print *, "Mean of data_in(50,50):", data_in(50,50)
									END IF
									
									sts = NF90_CLOSE(ncidin)
									sts = NF90_CLOSE(ncidin0) !close both files
									
									
									!-------------------------------------------------------------------------------
									! some analysis of the data
									
									print *,"shape of array" , SHAPE(temp_data)
									print *,"size of array" , SIZE(temp_data)
									!
									!       stat_mean = SUM(temp_data(:,:,5))/(MAX(1,SIZE(temp_data(:,:,5))))
									!       PRINT *, stat_mean
									stat_mean = SUM(temp_data(:,:))/SIZE(temp_data(:,:))
									PRINT *, stat_mean
									
									!-------------------------------------------------------------------------------
									! processing
									
									
									!       ***pr***
									IF (var_cmip(ivar) == "pr") THEN 
										
										temp_data(:,:) = ((rainnc_in(:,:,2) + rainc_in(:,:,2)) - (rainnc_in(:,:,1) + rainc_in(:,:,1)))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										!ATTENTION: implement adjustable time intervals that the differences are devided by
										
									END IF
									
									
									!       ***prc***
									IF (var_cmip(ivar) == "prc") THEN
										
										temp_data(:,:) = (rainc_in(:,:,2) - rainc_in(:,:,1))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***prsn***
									IF (var_cmip(ivar) == "prsn") THEN
										
										temp_data(:,:) = (snownc_in(:,:,2) - snownc_in(:,:,1))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***snm***
									IF (var_cmip(ivar) == "snm") THEN
										
										temp_data(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1))/(dtHours*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspsbl***
									IF (var_cmip(ivar) == "evspsbl") THEN
										
										temp_data(:,:) = (sfcevp_in(:,:,2) - sfcevp_in(:,:,1))/(dtHours*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspblpot**
									IF (var_cmip(ivar) == "evspsblpot") THEN
										
										temp_data(:,:) = (potevp_in(:,:,2) - (potevp_in(:,:,1)))/L   !unit [W m-2]/[J kg-1] -> [kg m-2 s-1]
										
										! THERE IS STH WRONG WITH THE UNITS: WRF's POTEVP is accumulated and declared to be in W m-2. 
										! It doesen't make sense to accumulate in W m-2, but even if assume it as J m-2 or derive kg m-2 
										! by using latent heat of vaporization you never get values that have a reasonable magnitude...
										
										
									END IF
									
									!       ***mrros***
									IF (var_cmip(ivar) == "mrros") THEN
										
										temp_data(:,:) = (sfroff_in(:,:,2) - sfroff_in(:,:,1))/(dtHours*3600.)       !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***mrro***
									IF (var_cmip(ivar) == "mrro") THEN
										
										temp_data(:,:) = ((sfroff_in(:,:,2) - sfroff_in(:,:,1)) + (udroff_in(:,:,2) - udroff_in(:,:,1)))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!!!!!!! This if has probably forgotten to include 3 variables (rlut, rsdt, rsut).
									!!!!!! They are read from the files earlier altogether, but have not been included here.
									!!!!!!!! Probably a mistake... add them some time...
									!       ***rsds, rlds, rsus, rlus***
									IF ( (var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")      &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) - rad_in(:,:,1)) /(dtHours*3600.)       ! take difference of accumulated values
										
										ELSE IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation_alternative") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) + rad_in(:,:,1)) / 2.              ! take mean of instantaneous values
											
										END IF 
										
									END IF
									
									
									!       ***mrso***
									IF (var_cmip(ivar) == "mrso") THEN
										
										temp_data(:,:) = (((smois_in(:,:,1,1)*0.1 + smois_in(:,:,2,1)*0.3 + smois_in(:,:,3,1)*0.6 + smois_in(:,:,4,1)*1.0 ) + &
										(smois_in(:,:,1,2)*0.1 + smois_in(:,:,2,2)*0.3 + smois_in(:,:,3,2)*0.6 + smois_in(:,:,4,2)*1.0 ))/2.)*1000. 
										
									END IF
									
									
									
									!Up to this point, temp_data has played the role of data_in... for 3-hr intervals (strictly!)
									
									!									temp_data_InTime(:,:,it) = temp_data(:,:)
									
									!								END DO !DO it = 1, ...
									!variable "it" will be InDimLenRec+1 at this point.
									
									!Locate the position in the file, to write the data
									
									counter = 0
									DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
										!									
										counter = counter + 1
										!									!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
										!									
										!									IF (dtHours.eq.24) THEN !we are averaging per day.. so compare DAYS
										!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										!										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) .AND. &
										!										( TimeRefArraySelYear(i,4) == InDateTimeDay(1)   ) ) THEN
										!											
										!											EXIT
										!											
										!										END IF
										!									
										!								ELSE IF (dtHours.eq.720) THEN !we are averaging per month... so compare Months
										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) ) THEN
											
											EXIT
											
										END IF
										!									ELSE IF (dtHours.eq.2880) THEN !we are averaging over seasons.. so ..................
										!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										!										( TimeRefArraySelYear(i,3) > InDateTimeMonth(1)  ) ) THEN
										!											
										!											EXIT
										!											
										!										END IF
										!									END IF
										!										
									END DO
									
									counter = counter + counter_prev_year
									PRINT *, "index where in the NC file the WRF data is sorted in = ", &
									counter
									
									!now the writing process takes place only once!
									!									DO writeIndex = 1,InDimLenRec,1
									!-------------------------------------------------------------------------------
									! write data to NetCDF file
									
									data_in = temp_data

									!!!!!!! ARIS CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!
									
									print *,'write data to NetCDF file'
									print *,'fn_out',fn_out
									sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid )
									IF (sts/=0) EXIT
									print *, 'NF90_OPEN',  sts
									sts = NF90_INQ_VARID(ncid, TRIM(var_cmip(ivar)), x_varid)
									print *, 'NF90_INQ_VARID', ncid
									print *, 'var_cmip(ivar)', var_cmip(ivar)
									print *, 'x_varid', x_varid
									print *, 'counter', counter
									print *, 'xfocus', xfocus
									print *, 'yfocus', yfocus
									IF ( height(ivar) /= -999 ) THEN
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1, 1 /) )
									
									ELSE
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									END IF
									
									print *,'NF90_PUT_VAR', sts
									print *, 'ncid', ncid
									print *, 'x_varid', x_varid
									print *, 'data_in(50:52,50:52)', data_in(50:52,50:52)
									sts = NF90_CLOSE(ncid)
									
									print *, pn_out//"/"//fn_out
									print *, TRIM(var_cmip(ivar)), xfocus, yfocus, counter, ncid, x_varid
									
									!-------------------------------------------------------------------------------
									
									!									END DO !writeIndex ...
									!single write process... no iteration here..									
									
									ifl = ifl + NoOfDaysPerMonth(InDateTimeMonth(1))

									
									!-------------------------------------------------------------------------------
									! next WRF file contains different number of output intervals
									
									DEALLOCATE(InVarDataRec)
									
									DEALLOCATE(InDateTimeYear)
									DEALLOCATE(InDateTimeMonth)
									DEALLOCATE(InDateTimeDay)
									DEALLOCATE(InDateTimeHour)
									
									!-------------------------------------------------------------------------------
									
									PRINT *, "------------------------------------------------------------"
									CALL CPU_TIME(cpuTe)
									PRINT '("CPU timing for 1 WRF file (e.g. 1 month worth of data) = ",F6.3," sec")',cpuTe-cpuTs
									
									
									
								END DO !ifl - specific WRF input file, filelist loop
								
								
								
								
							
							ELSE !dtHours == 2880
								
								
								
								
								!do for EVEN FEWER files.... (that is, 3...)
								ifl = 1
								DO WHILE (ifl < (SIZE(fl_wrfout) - 31))
									!							DO ifl = 1, SIZE(fl_wrfout)-1, 1 ! operational: loop over complete filelist
									
									print *,' SIZE(fl_wrfout', SIZE(fl_wrfout)
									
									!     DO ifl = 1, SIZE(fl_wrfxtr), 1 ! operational: loop over complete filelist
									!      print *,' SIZE(fl_wrfout', SIZE(fl_wrfxtr)
									
									!DO ifl = 1, 1, 1 ! testing: loop over specific entry in filelist (e.g. just January)
									
									CALL CPU_TIME(cpuTs)
									
									PRINT *, "------------------------------------------------------------"
									
									IF ( filetype(ivar) == "s" ) THEN
										
										iflWRFin = fl_wrfout(ifl)
									
									ELSE IF ( filetype(ivar) == "x" ) THEN
										
										iflWRFin = fl_wrfxtr(ifl)
										
									END IF
									
									PRINT '(100A)', TRIM(iflWRFin)
									
									!-------------------------------------------------------------------------------
									! which timespan is covered by the WRF outputs?
									! assume timespan wrfout = timespan wrfxtrm
									! this determines how many times the tool has to loop over the inputs
									! also check how many years are covered by a single wrfout and wrfxtrm which
									! determines the output file generation
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncid_in)
									sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
									sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
									sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
									NAME = InDimNameRec, LEN = InDimLenRec)
									ALLOCATE(InVarDataRec(InDimLenRec))
									sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
									sts = NF90_CLOSE(ncid_in)
									
									print *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)
									
									ALLOCATE(InDateTimeYear(InDimLenRec))
									ALLOCATE(InDateTimeMonth(InDimLenRec))
									ALLOCATE(InDateTimeDay(InDimLenRec))
									ALLOCATE(InDateTimeHour(InDimLenRec))
									
									! this is the temporal coverage of the WRF input data
									DO i = 1, SIZE(InVarDataRec), 1
										
										READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(i), InDateTimeMonth(i), InDateTimeDay(i), InDateTimeHour(i)
										
									END DO
									
									
									!-------------------------------------------------------------------------------
									! check whether a file is needed at all
									! check for first date
									! for subsequent dates, just check whether it is the same as the previous one
									! if this is not the case, then check whether file exists...
									! if it does not exist (default case): create
									! pathname is always needed
									
									! loop over the individual timesteps in the WRF files...
									! this may take some time but it is robust and also extremely large arrays
									! might be used
									! highly robust code
									
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ARIS EDIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									! Allocate correspondingly the temp_data_InTime variable
									
									
									!	SELECT CASE (frequency(ifrq))
									!		CASE ('3hr')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,1), STAT = sts) !run through the entire list of 3-hour records and use them all for the calculations.
									!		CASE ('day')
									!			ALLOCATE( temp_data_InTime(xfocus,yfocus,8), STAT = sts) !only run once... the daily data will be processed once..
									!		CASE DEFAULT
									!			PRINT *, "something happened... "
									!			total_iters = InDimLenRec !make sure bad results in default behavior (i.e. what was happening anyway before...)
									!			!STOP
									!	END SELECT
									
									
									
									IF ( CheckForLeapyear( InDateTimeYear(1), dataIncludesLeapYearDays ) == 366 ) THEN
										
										NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
									
									ELSE
										
										NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
										
									END IF
									
									!NOTE: THIS 'IF' CONDITION BELOW WILL IGNORE FIRST 2 MONTHS IF INPUT DATA
									!DOES NOT INCLUDE DECEMBER FROM THE YEAR PRIOR TO THE STARTING YEAR OF THIS RUN
									!(WHICH IS TYPICALLY THE CASE). THEREFORE CALCULATIONS WILL BEGIN FROM MARCH,
									!AND OUTPUT DATA FROM FIRST SEASON (WINTER [DEC/JAN/FEB]) WILL BE FILLED IN WITH
									!MISSING VALUES.

									
									!if the execution is on the beginning of a year, then skip two months
									!because they are an incomplete season..
									IF (ifl < 58) THEN
										ifl = ifl + NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2)
										!start all over, after the already allocated arrays above
										!are deallocated to avoid segmentation faults..
										DEALLOCATE(InVarDataRec)
										
										DEALLOCATE(InDateTimeYear)
										DEALLOCATE(InDateTimeMonth)
										DEALLOCATE(InDateTimeDay)
										DEALLOCATE(InDateTimeHour)
										
										CYCLE !restart to next iteration
									END IF
									
									!!!! Here is the case, when 3-hr intervals have been requested... so no alterations are made,
									!!!! except for removing the irrelevant variables, keeping only the cumulative ones..
									!								DO it = 1, InDimLenRec, 1
									it = 1
									PRINT *, "----------------------------------------"
									PRINT *, "working on date in WRF input file = ", TRIM(InVarDataRec(it))
									PRINT *, "---"
									
									!-------------------------------------------------------------------------------
									! generate path and filename
									! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/EUR-44/CRPGL/
									! ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
									! evspsbl_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_CRPGL-WRFARW331_v1_3hr_
									! 1989010100-1989123121
									
									! output does not yet exist
									! monthy check is basically not even possible, but does not do any harm
									! this is a rerstriction for all those who might have a different file
									! structure
									IF ( ( InDateTimeYearPrev /= InDateTimeYear(it) ) .AND. ( InDateTimeYear(it) /= LastYearOfRun + 1 ) ) THEN !.AND. &
										!( InDateTimeMonthPrev /= InDateTimeMonth(it) ) ) THEN
										
										InDateTimeYearPrev = InDateTimeYear(it)
										
										!This shouldn't happen similarly to daily or monthly averaging since
										!it needs to happen on every change of season (and not on change of year)
										!counter_prev_year = counter
										
										
										PRINT *, "start of processing or new year encountered -> t ref. vec. and filecheck"
										
										!Change naming rules for producing just one output file when running for multiple years
										!WRITE (EndInDateTimeYearStr,'(I4.4)') InDateTimeYear(it)+1
										WRITE (InDateTimeYearStr,'(I4)') FirstYearOfRun
										!WRITE (InDateTimeYearStr,'(I4.4)') InDateTimeYear(it)		
										WRITE (EndInDateTimeYearStr,'(I4)') LastYearOfRun
										
									
										!WRITE (InDateTimeMonthStr,'(I2.2)') InDateTimeMonth(it)
										!PRINT *, InDateTimeYearStr, InDateTimeMonthStr
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										IF ( cell_methods(ivar) == "mean" ) THEN
											
											WRITE (EndInDateTimeYearStr,'(I4)') (LastYearOfRun + 1)
											EndDate=EndInDateTimeYearStr//"010100"
											
										
										ELSE
											
											EndDate=EndInDateTimeYearStr//"123121"
											
										END IF
										
										
										pn_out = TRIM(product)                       // "/" // &
										TRIM(CORDEX_domain)                 // "/" // &
										TRIM(institute_id)                  // "/" // &
										TRIM(driving_model_id)              // "/" // &
										TRIM(driving_experiment_name)       // "/" // &
										TRIM(driving_model_ensemble_member) // "/" // &
										TRIM(model_id)                      // "/" // &
										TRIM(rcm_version_id)                // "/" // &
										TRIM(frequency(ifrq))               // "/" // &
										TRIM(var_cmip(ivar))
										
										fn_out = TRIM(var_cmip(ivar))                // "_" // &
										TRIM(CORDEX_domain)                 // "_" // &
										TRIM(driving_model_id)              // "_" // &
										TRIM(driving_experiment_name)       // "_" // &
										TRIM(driving_model_ensemble_member) // "_" // &
										TRIM(model_id)                      // "_" // &
										TRIM(rcm_version_id)                // "_" // &
										TRIM(frequency(ifrq))               // "_" // &
										!                   InDateTimeYearStr//"010100-"//InDateTimeYearStr//"123121" // &
										InDateTimeYearStr//"010100-"//TRIM(EndDate)//  &
										".nc"
										
										PRINT *, "pn_out = ", TRIM(pn_out)
										PRINT *, "fn_out = ", TRIM(fn_out)
										
										!-------------------------------------------------------------------------------
										! extract the time info from the ref array which fits the respective year
										! ...as there is no "WHERE" the way I need it in F95, use loops
										! this is needed whenever a new netcdf file is to be used, also if
										! this file exists already
										!        READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										DEALLOCATE( TimeRefArraySelYear )
										
										print *, "DEALLOCATE( TimeRefArraySelYear )" 
										
										DEALLOCATE( Time_bnds ) !SKn 
										
										print *, "DEALLOCATE( Time_bnds )"
										
										counter = 0
										PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
										PRINT *,'SHAPE(TimeRefArray, 1)',  SHAPE(TimeRefArray, 1)
										PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
										PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
										PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												
											END IF
											
										END DO
										PRINT *, "Counter: ", counter
										! holds data of exactly 1 year
										PRINT *, "timesteps in the time ref. subset = ", counter
										ALLOCATE( TimeRefArraySelYear( counter, 5 ) ) ! index, y, m, d, h
										
										print *, "ALLOCATE( TimeRefArraySelYear( counter, 5 ) )"
										
										print *, "TimeRefArraySelYear( counter, 1 )", TimeRefArraySelYear( counter, 1 )
										print *, "TimeRefArraySelYear( counter, 2 )", TimeRefArraySelYear( counter, 2 )
										print *, "TimeRefArraySelYear( counter, 3 )", TimeRefArraySelYear( counter, 3 )
										print *, "TimeRefArraySelYear( counter, 4 )", TimeRefArraySelYear( counter, 4 )
										print *, "TimeRefArraySelYear( counter, 5 )", TimeRefArraySelYear( counter, 5 )
										
										! find the matching elements of the respecitve year and copy them
										
										ALLOCATE( Time_bnds( 2, counter ) )
										
										print *, "ALLOCATE( Time_bnds( 2, counter ) )", Time_bnds( 2, counter )
										
										counter = 0
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												TimeRefArraySelYear(counter,1:5) = TimeRefArray(i,1:5)
												
											END IF
											
										END DO
										
										!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySelYear(:,:) )
										
										!-------------------------------------------------------------------------------
										! check for existance of the file and generate file if needed
										
										INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )
										
										! could exist already from a previous run of tool and due to multiple
										! months in a file (i.e. 1 WRF output file may cover several months)
										IF ( FileExists ) THEN
											
											PRINT *, "path and file exist, continue filling"
											
											!!!!Change indices for TimeRefArraySelYear (time) and Time_bnds (time bounds) records
											!!!!so that they are written in proper position in same file right after the last year's data
											sts = NF90_OPEN(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid)
											
											!Not sure if all of this is necessary though...
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.), START = (/ counter_prev_year + 1, 1 /) , COUNT = (/ SIZE(TimeRefArraySelYear(:,1), 1) /) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, counter_prev_year+1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF
											
											!!!!!!!!!!!!!!!!!!!!!!!!
											!!!I'm pretty sure that this is not necessary here...
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											!!!!!!!!!!!!!!!!!!!!!!!!
											
											sts = NF90_CLOSE(ncid)
											
											!!!!End adjustments regarding index w.r.t. Time and Time_bnds when running for consecutive years
																				
										ELSE
											
											PRINT *, "path and file do not yet exist, create path and NetCDF file first"
											PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)
											
											CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )
											
											!-------------------------------------------------------------------------------
											! non-std, works for gfortran (fct & subroutine) + ifort
											! comment lines in the NetCDF file global attribute definition
											! turn standard checking in Makefile off
											! trackingID = "xxxxxxxx-xxxx-Mxxx-Nxxx-xxxxxxxxxxxx"
											! creationDate = "YYYY-MM-DD-THH:MM:SSZ"
											
											CALL SYSTEM(cmdUUID)
											OPEN(1,FILE="tmpfileUUID",STATUS='old')
											READ(1,*) trackingID
											CLOSE(1)
											PRINT *, "uuidgen externally generated trackingID = ", trackingID
											
											CALL SYSTEM(cmdDate)
											OPEN(1,FILE="tmpfileDate",STATUS='old')
											READ(1,*) creationDate
											CLOSE(1)
											PRINT *, "date externally generated creation date = ", creationDate
											
											!-------------------------------------------------------------------------------
											! create NetCDF file
											! NF90_CLASSIC_MODEL = NetCDF4_classic
											! NF90_HDF5 = NetCDF4 based on HDF5
											! NF90_CLOBBER = old NetCDF
											! sts = NF90_CREATE(PathFileNameOutTEST, NF90_HDF5, ncid)
											
											!comb_flags = IOR(NF90_HDF5, NF90_CLASSIC_MODEL)
											!https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/NF90_005fCREATE.html
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid)   !not sure whether this is the right data format spec. I guess it may be right using compression but not the other fancy stuff from NetCDF4
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncid)   !if anything, then use this here
											sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_HDF5, ncid)
											
											
											! always included
											sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
											sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid)
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
												
												ELSE
													
													sts = NF90_DEF_DIM(ncid, "plev", 1, height_dimid)
													
												END IF
												
											END IF
											sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
												
											END IF
											
											
											! always included
											sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lon_varid)
											sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "long_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
											
											! always included
											sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lat_varid)
											sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "long_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
											
											! always included
											sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid)
											sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
											sts = nf90_put_att(ncid, rlon_varid, "long_name", "longitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
											
											! always included
											sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid)
											sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
											sts = nf90_put_att(ncid, rlat_varid, "long_name", "latitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")
											
											! always included
											! restriction to one domain only
											sts = nf90_def_var(ncid, "rotated_pole", NF90_DOUBLE, rotated_pole_varid)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", 39.25)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", -162.0)
											
											! depends whether height is set in the nml
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_def_var(ncid, "height", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
													sts = nf90_put_att(ncid, height_varid, "long_name", "height")
													sts = nf90_put_att(ncid, height_varid, "units", "m")
													sts = nf90_put_att(ncid, height_varid, "positive", "up")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")
												
												ELSE
													
													sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name","air_pressure")
													sts = nf90_put_att(ncid, height_varid, "long_name", "pressure")
													sts = nf90_put_att(ncid, height_varid, "units", "Pa")
													sts = nf90_put_att(ncid, height_varid, "positive", "down")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")                 
													
												END IF
												
											END IF
											
											!missing: lvl, depends
											
											! always included
											sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid)
											sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "long_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "units", "days since 1949-12-01T00:00:00Z")
											sts = nf90_put_att(ncid, rec_varid, "calendar", "standard")
											sts = nf90_put_att(ncid, rec_varid, "axis", "T")
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												! Add time bounds to time varaible
												
												sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")	
												sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid)
												sts = nf90_put_att(ncid, recbnds_varid, "standard_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "long_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "units", "days since 1949-12-01T00:00:00Z")
												sts = nf90_put_att(ncid, recbnds_varid, "calendar", "standard")
												!sts = nf90_put_att(ncid, recbnds_varid, "axis", "T")
												
												print *,'rec_varid', rec_varid
												print *,'recbnds_varid', recbnds_varid
												
											END IF
											
											! always included
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact", contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "creation_date", creationDate)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment", experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment_id", experiment_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment", driving_experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_id", driving_model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_ensemble_member", driving_model_ensemble_member)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment_name", driving_experiment_name)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "frequency", frequency(ifrq))
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution", institution)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_id", institute_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "model_id", model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "forcing",forcing)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact",contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source",source)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "rcm_version_id", rcm_version_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "CORDEX_domain", CORDEX_domain)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id", trackingID)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "title",title)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_run_id", institute_run_id)
											
											! always included
											IF ( height(ivar) /= -999 ) THEN
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, height_dimid, rec_dimid /), x_varid)
											
											ELSE
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid)
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
											
											IF ( positive(ivar) /= '-999' ) THEN
												
												sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
											
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat height")
												
												ELSE
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat plev")
													
												END IF
											
											ELSE
												
												sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat")
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "grid_mapping", "rotated_pole")
											
											!IF ( height(ivar) == 850 ) THEN
											sts = nf90_put_att(ncid, x_varid, "missing_value", 1.e20)
											sts = nf90_put_att(ncid, x_varid, "_FillValue", 1.e20)
											!END IF
											
											sts = NF90_ENDDEF(ncid)
											
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF
											!print *,'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
											
											
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											
											! add time_bnds, calc here
											
											! add height from NML
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
												
												ELSE
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar)*100. )
													
												END IF
												
											END IF
											
											sts = NF90_CLOSE(ncid)
											
										END IF ! file exists y/n
										
									END IF ! checking with previous year and month
									
									!-------------------------------------------------------------------------------
									! match timestep of WRFin with the subset of the ref time vec which belong to
									! the NetCDF file of the year currently open to receive data
									! NOT SURE WHETHER THIS IS NEEDED AT ALL, THIS IS THE WRONG DIRECTION ???
									! see whether the current time of the timestep fits anywhere in the
									
									
									!-------------------------------------------------------------------------------
									! read orig WRF outpoutses
									! there is always a corresponding time-slot in the NC file
									! extracted time from above
									! "it" controls it all: timestep in the individual WRF file
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)
									
									!< and NOT <= is used, because we want to avoid the last december..
									
									if (ifl < (SIZE(fl_wrfout,1)-31) ) THEN
										sts = NF90_OPEN(fl_wrfout(ifl + NoOfDaysPerMonth(InDateTimeMonth(1)) + &
										NoOfDaysPerMonth(InDateTimeMonth(1)+1) + &
										NoOfDaysPerMonth(InDateTimeMonth(1)+2) ), &
										NF90_NOWRITE, ncidin0)
									END IF
									!InDateTimeMonth is the month in which the first 3-hr interval refers to.
									!All intervals refer to the same month anyway, so the first is used.
									!According to whether it is a leapyear or not, the appropriate number of days are addded...
									!in order to refer to the first day of the next month..
									
									!!!!!!!!!!!!! ONLY READ THE VARIABLES THAT PROVIDE CUMULATIVE WRF OUTPUTS.......!!!
									IF (var_cmip(ivar) == "pr") THEN 
										
										ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										!get first day
										sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus, 1 /) )   !read two timesteps to calculate 3hr sum
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin0, rainnc_varid, rainnc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )   !read last timestep of previous wrfout file
										
										sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prc") THEN
										
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prsn") THEN
										
										ALLOCATE( snownc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin0, snownc_varid, snownc_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "snm") THEN
										
										ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin0, acsnom_varid, acsnom_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), &
										COUNT = (/ xfocus, yfocus,1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "evspsbl") THEN
										
										ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfcevp_varid, sfcevp_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "evspsblpot") THEN
										
										ALLOCATE( potevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin0, potevp_varid, potevp_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrros") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrro") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin0, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin0, udroff_varid, udroff_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF ((var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")  &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "rlut")                                 &
									.or. (var_cmip(ivar) == "rsdt") .or. (var_cmip(ivar) == "rsut") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										ALLOCATE( rad_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,1), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin0, varid, rad_in(:,:,2), &
										START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										
										print*, var_cmip(ivar), rad_in(50,50,1), rad_in(50,50,2)
										print*, 'difference in J m-2', (rad_in(50,50,2) - rad_in(50,50,1))
										print*, 'in mean W m-2', (rad_in(50,50,2) - rad_in(50,50,1))/ (obs_interval*3600.)
										
										! alternative: since accumulated values as read above get so large in 
										! long term simulations that their differences loose accuracy, use 
										! instantaneous values instead and calculate means
										
									
									ELSE IF (var_cmip(ivar) == "mrso") THEN
										
										ALLOCATE( smois_in( xfocus, yfocus, 4, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,1), &
										START = (/ xoffset, yoffset, 1, 1 /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
										
										!get next day
										sts = NF90_INQ_VARID(ncidin0, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin0, varid, smois_in(:,:,:,2), &
										START = (/ xoffset, yoffset, 1, 1 /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
										
										
										!!!!!!!!!!!!!!! THIS ELSE IS NEVER GOING TO BE EXECUTED... IT IS JUST COMMENTED OUT
										!									ELSE 
										!										
										!										ALLOCATE( temp_data ( xfocus, yfocus ), STAT=sts )
										!										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										!										!print *, "VALUE FOR NAME: ", TRIM(var_wrf(ivar)), "+ VARID: ", varid
										!										
										!										
										!										sts = NF90_GET_VAR(ncidin, varid, temp_data(:,:), &
										!										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										!										!			data_in(:,:) = SUM(temp_data, 3)/8.
										!										!			print *, temp_data(50,50,1), temp_data(50,50,2), temp_data(50,50,3), temp_data(50,50,4)
										!										!			print *, temp_data(50,50,5), temp_data(50,50,6), temp_data(50,50,7), temp_data(50,50,8)
										!										!			print *, "Mean of data_in(50,50):", data_in(50,50)
									END IF
									
									sts = NF90_CLOSE(ncidin)
									sts = NF90_CLOSE(ncidin0) !close both files
									
									
									!-------------------------------------------------------------------------------
									! some analysis of the data
									
									print *,"shape of array" , SHAPE(temp_data)
									print *,"size of array" , SIZE(temp_data)
									!
									!       stat_mean = SUM(temp_data(:,:,5))/(MAX(1,SIZE(temp_data(:,:,5))))
									!       PRINT *, stat_mean
									stat_mean = SUM(temp_data(:,:))/SIZE(temp_data(:,:))
									PRINT *, stat_mean
									
									!-------------------------------------------------------------------------------
									! processing
									
									
									!       ***pr***
									IF (var_cmip(ivar) == "pr") THEN 
										
										temp_data(:,:) = ((rainnc_in(:,:,2) + rainc_in(:,:,2)) - (rainnc_in(:,:,1) + rainc_in(:,:,1)))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										!ATTENTION: implement adjustable time intervals that the differences are devided by
										
									END IF
									
									
									!       ***prc***
									IF (var_cmip(ivar) == "prc") THEN
										
										temp_data(:,:) = (rainc_in(:,:,2) - rainc_in(:,:,1))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***prsn***
									IF (var_cmip(ivar) == "prsn") THEN
										
										temp_data(:,:) = (snownc_in(:,:,2) - snownc_in(:,:,1))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***snm***
									IF (var_cmip(ivar) == "snm") THEN
										
										temp_data(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1))/(dtHours*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspsbl***
									IF (var_cmip(ivar) == "evspsbl") THEN
										
										temp_data(:,:) = (sfcevp_in(:,:,2) - sfcevp_in(:,:,1))/(dtHours*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspblpot**
									IF (var_cmip(ivar) == "evspsblpot") THEN
										
										temp_data(:,:) = (potevp_in(:,:,2) - (potevp_in(:,:,1)))/L   !unit [W m-2]/[J kg-1] -> [kg m-2 s-1]
										
										! THERE IS STH WRONG WITH THE UNITS: WRF's POTEVP is accumulated and declared to be in W m-2. 
										! It doesen't make sense to accumulate in W m-2, but even if assume it as J m-2 or derive kg m-2 
										! by using latent heat of vaporization you never get values that have a reasonable magnitude...
										
										
									END IF
									
									!       ***mrros***
									IF (var_cmip(ivar) == "mrros") THEN
										
										temp_data(:,:) = (sfroff_in(:,:,2) - sfroff_in(:,:,1))/(dtHours*3600.)       !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***mrro***
									IF (var_cmip(ivar) == "mrro") THEN
										
										temp_data(:,:) = ((sfroff_in(:,:,2) - sfroff_in(:,:,1)) + (udroff_in(:,:,2) - udroff_in(:,:,1)))/(dtHours*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!!!!!!! This if has probably forgotten to include 3 variables (rlut, rsdt, rsut).
									!!!!!! They are read from the files earlier altogether, but have not been included here.
									!!!!!!!! Probably a mistake... add them some time...
									!       ***rsds, rlds, rsus, rlus***
									IF ( (var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")      &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) - rad_in(:,:,1)) /(dtHours*3600.)       ! take difference of accumulated values
										
										ELSE IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation_alternative") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) + rad_in(:,:,1)) / 2.              ! take mean of instantaneous values
											
										END IF 
										
									END IF
									
									
									!       ***mrso***
									IF (var_cmip(ivar) == "mrso") THEN
										
										temp_data(:,:) = (((smois_in(:,:,1,1)*0.1 + smois_in(:,:,2,1)*0.3 + smois_in(:,:,3,1)*0.6 + smois_in(:,:,4,1)*1.0 ) + &
										(smois_in(:,:,1,2)*0.1 + smois_in(:,:,2,2)*0.3 + smois_in(:,:,3,2)*0.6 + smois_in(:,:,4,2)*1.0 ))/2.)*1000. 
										
									END IF
									
									
									
									!Up to this point, temp_data has played the role of data_in... for 3-hr intervals (strictly!)
									
									!									temp_data_InTime(:,:,it) = temp_data(:,:)
									
									!								END DO !DO it = 1, ...
									!variable "it" will be InDimLenRec+1 at this point.
									
									!Locate the position in the file, to write the data
									
									counter = 0
									DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
										!									
										counter = counter + 1
										!									!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
										!									
										!									IF (dtHours.eq.24) THEN !we are averaging per day.. so compare DAYS
										!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										!										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) .AND. &
										!										( TimeRefArraySelYear(i,4) == InDateTimeDay(1)   ) ) THEN
										!											
										!											EXIT
										!											
										!										END IF
										!									
										!									ELSE IF (dtHours.eq.720) THEN !we are averaging per month... so compare Months
										!										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										!										( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) ) THEN
										!											
										!											EXIT
										!											
										!										END IF
										!									ELSE IF (dtHours.eq.2880) THEN !we are averaging over seasons.. so ..................
										IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
										( TimeRefArraySelYear(i,3) > InDateTimeMonth(1)  ) ) THEN
											
												!Setting counter manually to zero at this point, because when computing seasonal avg. then
												!counter_prev_year should come into play in end of last season of year,
												!and not in end of year (in contrast to daily/monthly avg.).
												!Otherwise counter will be added twice right after this loop.
												IF (currentMonthOfYear == 12) THEN
																							
													counter = 1
													EXIT

												END IF
											
										END IF
										!									END IF
										!										
									END DO
									
									counter = counter + counter_prev_year									
									PRINT *, "index where in the NC file the WRF data is sorted in = ", &
									counter
									
									!Adjusting counter here specifically, because when computing seasonal avg. then
									!previous counter should be stored and come into play in end of last season of year,
									!and not in end of year (in contrast to daily/monthly avg.).
									IF (currentDayOfMonth == NoOfDaysPerMonth(currentMonthOfYear) .AND. currentMonthOfYear == 11) THEN 
										
										counter_prev_year = counter
									
									END IF	
									
									!now the writing process takes place only once!
									!									DO writeIndex = 1,InDimLenRec,1
									!-------------------------------------------------------------------------------
									! write data to NetCDF file
									
									data_in = temp_data 
									!!!!!!! ARIS CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!
									
									print *,'write data to NetCDF file'
									print *,'fn_out',fn_out
									sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid )
									IF (sts/=0) EXIT
									print *, 'NF90_OPEN',  sts
									sts = NF90_INQ_VARID(ncid, TRIM(var_cmip(ivar)), x_varid)
									print *, 'NF90_INQ_VARID', ncid
									print *, 'var_cmip(ivar)', var_cmip(ivar)
									print *, 'x_varid', x_varid
									print *, 'counter', counter
									print *, 'xfocus', xfocus
									print *, 'yfocus', yfocus
									IF ( height(ivar) /= -999 ) THEN
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1, 1 /) )
									
									ELSE
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									END IF
									
									print *,'NF90_PUT_VAR', sts
									print *, 'ncid', ncid
									print *, 'x_varid', x_varid
									print *, 'data_in(50:52,50:52)', data_in(50:52,50:52)
									sts = NF90_CLOSE(ncid)
									
									print *, pn_out//"/"//fn_out
									print *, TRIM(var_cmip(ivar)), xfocus, yfocus, counter, ncid, x_varid
									
									!-------------------------------------------------------------------------------
									
									!									END DO !writeIndex ...
									!single write process... no iteration here..									
									
									!Jump the next 3 months (i.e. the season)
									ifl = ifl + NoOfDaysPerMonth(InDateTimeMonth(1)) + &
									NoOfDaysPerMonth(InDateTimeMonth(1)+1) + &
									NoOfDaysPerMonth(InDateTimeMonth(1)+2)
									
									!-------------------------------------------------------------------------------
									! next WRF file contains different number of output intervals
									
									DEALLOCATE(InVarDataRec)
									
									DEALLOCATE(InDateTimeYear)
									DEALLOCATE(InDateTimeMonth)
									DEALLOCATE(InDateTimeDay)
									DEALLOCATE(InDateTimeHour)
									
									!-------------------------------------------------------------------------------
									
									PRINT *, "------------------------------------------------------------"
									CALL CPU_TIME(cpuTe)
									PRINT '("CPU timing for 1 WRF file (e.g. 1 month worth of data) = ",F6.3," sec")',cpuTe-cpuTs
									
									
									
								END DO !ifl - specific WRF input file, filelist loop
								
								
								
								
							END IF
							
						
						ELSE !variable does not require multiple files for averaging...
							!so just read daily values and calculate daily averages and then
							!whatever type of averages are desired...
							
							
							!THOSE THAT DO NOT REQUIRE CUMULATIVE CALCULATIONS - OPENING PAIRS OF FILES, ETC....
							
							DO ifl = 1, SIZE(fl_wrfout), 1 ! operational: loop over complete filelist
								
								print *,' SIZE(fl_wrfout', SIZE(fl_wrfout)
								
								!     DO ifl = 1, SIZE(fl_wrfxtr), 1 ! operational: loop over complete filelist
								!      print *,' SIZE(fl_wrfout', SIZE(fl_wrfxtr)
								
								!DO ifl = 1, 1, 1 ! testing: loop over specific entry in filelist (e.g. just January)
								
								CALL CPU_TIME(cpuTs)
								
								PRINT *, "------------------------------------------------------------"
								
								IF ( filetype(ivar) == "s" ) THEN
									
									iflWRFin = fl_wrfout(ifl)
								
								ELSE IF ( filetype(ivar) == "x" ) THEN
									
									iflWRFin = fl_wrfxtr(ifl)
									
								END IF
								
								PRINT '(100A)', TRIM(iflWRFin)
								
								!-------------------------------------------------------------------------------
								! which timespan is covered by the WRF outputs?
								! assume timespan wrfout = timespan wrfxtrm
								! this determines how many times the tool has to loop over the inputs
								! also check how many years are covered by a single wrfout and wrfxtrm which
								! determines the output file generation
								
								sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncid_in)
								sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
								sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
								sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
								NAME = InDimNameRec, LEN = InDimLenRec)
								ALLOCATE(InVarDataRec(InDimLenRec))
								sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
								sts = NF90_CLOSE(ncid_in)
								
								print *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)
								
								ALLOCATE(InDateTimeYear(InDimLenRec))
								ALLOCATE(InDateTimeMonth(InDimLenRec))
								ALLOCATE(InDateTimeDay(InDimLenRec))
								ALLOCATE(InDateTimeHour(InDimLenRec))
								
								! this is the temporal coverage of the WRF input data
								DO i = 1, SIZE(InVarDataRec), 1
									
									READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(i), InDateTimeMonth(i), InDateTimeDay(i), InDateTimeHour(i)
									
								END DO
								
								
								!-------------------------------------------------------------------------------
								! check whether a file is needed at all
								! check for first date
								! for subsequent dates, just check whether it is the same as the previous one
								! if this is not the case, then check whether file exists...
								! if it does not exist (default case): create
								! pathname is always needed
								
								! loop over the individual timesteps in the WRF files...
								! this may take some time but it is robust and also extremely large arrays
								! might be used
								! highly robust code
								
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ARIS EDIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								! Allocate correspondingly the temp_data_InTime variable
								
								
								!	SELECT CASE (frequency(ifrq))
								!		CASE ('3hr')
								!			ALLOCATE( temp_data_InTime(xfocus,yfocus,1), STAT = sts) !run through the entire list of 3-hour records and use them all for the calculations.
								!		CASE ('day')
								!			ALLOCATE( temp_data_InTime(xfocus,yfocus,8), STAT = sts) !only run once... the daily data will be processed once..
								!		CASE DEFAULT
								!			PRINT *, "something happened... "
								!			total_iters = InDimLenRec !make sure bad results in default behavior (i.e. what was happening anyway before...)
								!			!STOP
								!	END SELECT
								
								record_count = dtHours / obs_interval
								!if record_count is NOT equal to 1, then some type of temporal average has been requested...
								
								ALLOCATE( temp_data_InTime(xfocus,yfocus,InDimLenRec), STAT = sts)
							
								ALLOCATE( counter_array(InDimLenRec), STAT = sts)
								
								IF ( CheckForLeapyear( InDateTimeYear(1), dataIncludesLeapYearDays ) == 366 ) THEN
									
									NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
								
								ELSE
									
									NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
									
								END IF
								
								
								IF (isEndOfMonth) THEN
									
									IF (ALLOCATED(temp_data_Monthly)) THEN
									
										!print*, "ALLOCATED MUCH???", temp_data_Monthly
										
										DEALLOCATE (temp_data_Monthly)
										
									END IF
									
									ALLOCATE( temp_data_Monthly(xfocus,yfocus, NoOfDaysPerMonth(InDateTimeMonth(1))  ), STAT = sts  )
									isEndOfMonth = 0
									
								END IF
								
								IF (isEndOfSeason) THEN
									
									IF (ALLOCATED(temp_data_Seasonally)) THEN
										
										DEALLOCATE (temp_data_Seasonally)
										
									END IF
									
									ALLOCATE( temp_data_Seasonally(xfocus, yfocus, 3, 31 ), STAT = sts)
									isEndOfSeason = 0
									dayCount = 0
									
								END IF
								
								DO it = 1, InDimLenRec, 1
									
									PRINT *, "----------------------------------------"
									PRINT *, "working on date in WRF input file = ", TRIM(InVarDataRec(it))
									PRINT *, "---"
									
									!-------------------------------------------------------------------------------
									! generate path and filename
									! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/EUR-44/CRPGL/
									! ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
									! evspsbl_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_CRPGL-WRFARW331_v1_3hr_
									! 1989010100-1989123121
									
									! output does not yet exist
									! monthy check is basically not even possible, but does not do any harm
									! this is a rerstriction for all those who might have a different file
									! structure
	
									
									IF ( ( InDateTimeYearPrev /= InDateTimeYear(it) ) .AND. ( InDateTimeYear(it) /= LastYearOfRun + 1 ) ) THEN !.AND. &
										!( InDateTimeMonthPrev /= InDateTimeMonth(it) ) ) THEN
										

										InDateTimeYearPrev = InDateTimeYear(it)
										
										!When averaging data seasonally, then the counter should be stored and changed
										!on every change of season, and not on every change of year
										IF (dtHours /= 2880) THEN
											counter_prev_year = counter
										END IF
										
										PRINT *, "start of processing or new year encountered -> t ref. vec. and filecheck"
										
										!READ( InDateTimeYear(it), '(4A)' ) InDateTimeYearStr
										!READ( InDateTimeMonth(it), '(2A)' ) InDateTimeMonthStr
										
										!Change naming rules for producing just one output file when running for multiple years
										!WRITE (EndInDateTimeYearStr,'(I4.4)') InDateTimeYear(it)+1
										WRITE (InDateTimeYearStr,'(I4)') FirstYearOfRun
										!WRITE (InDateTimeYearStr,'(I4.4)') InDateTimeYear(it)		
										WRITE (EndInDateTimeYearStr,'(I4)') LastYearOfRun
										
									
										!WRITE (InDateTimeMonthStr,'(I2.2)') InDateTimeMonth(it)
										!PRINT *, InDateTimeYearStr, InDateTimeMonthStr
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										IF ( cell_methods(ivar) == "mean" ) THEN
											
											WRITE (EndInDateTimeYearStr,'(I4)') (LastYearOfRun + 1)
											EndDate=EndInDateTimeYearStr//"010100"
											
										
										ELSE
											
											EndDate=EndInDateTimeYearStr//"123121"
											
										END IF
										
										
										pn_out = TRIM(product)                       // "/" // &
										TRIM(CORDEX_domain)                 // "/" // &
										TRIM(institute_id)                  // "/" // &
										TRIM(driving_model_id)              // "/" // &
										TRIM(driving_experiment_name)       // "/" // &
										TRIM(driving_model_ensemble_member) // "/" // &
										TRIM(model_id)                      // "/" // &
										TRIM(rcm_version_id)                // "/" // &
										TRIM(frequency(ifrq))               // "/" // &
										TRIM(var_cmip(ivar))
										
										fn_out = TRIM(var_cmip(ivar))                // "_" // &
										TRIM(CORDEX_domain)                 // "_" // &
										TRIM(driving_model_id)              // "_" // &
										TRIM(driving_experiment_name)       // "_" // &
										TRIM(driving_model_ensemble_member) // "_" // &
										TRIM(model_id)                      // "_" // &
										TRIM(rcm_version_id)                // "_" // &
										TRIM(frequency(ifrq))               // "_" // &
										!                   InDateTimeYearStr//"010100-"//InDateTimeYearStr//"123121" // &
										InDateTimeYearStr//"010100-"//TRIM(EndDate)//  &
										".nc"
										
										PRINT *, "pn_out = ", TRIM(pn_out)
										PRINT *, "fn_out = ", TRIM(fn_out)
										
										!-------------------------------------------------------------------------------
										! extract the time info from the ref array which fits the respective year
										! ...as there is no "WHERE" the way I need it in F95, use loops
										! this is needed whenever a new netcdf file is to be used, also if
										! this file exists already
										!        READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2)' ) InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
										
										PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
										SHAPE(TimeRefArray)
										
										DEALLOCATE( TimeRefArraySelYear )
										
										print *, "DEALLOCATE( TimeRefArraySelYear )" 
										
										DEALLOCATE( Time_bnds ) !SKn 
										
										print *, "DEALLOCATE( Time_bnds )"
										
										counter = 0
										PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
										PRINT *,'SHAPE(TimeRefArray, 1)',  SHAPE(TimeRefArray, 1)
										PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
										PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
										PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												
											END IF
											
										END DO
										PRINT *, "Counter: ", counter
										! holds data of exactly 1 year
										PRINT *, "timesteps in the time ref. subset = ", counter
										ALLOCATE( TimeRefArraySelYear( counter, 5 ) ) ! index, y, m, d, h
										
										print *, "ALLOCATE( TimeRefArraySelYear( counter, 5 ) )"
										
										print *, "TimeRefArraySelYear( counter, 1 )", TimeRefArraySelYear( counter, 1 )
										print *, "TimeRefArraySelYear( counter, 2 )", TimeRefArraySelYear( counter, 2 )
										print *, "TimeRefArraySelYear( counter, 3 )", TimeRefArraySelYear( counter, 3 )
										print *, "TimeRefArraySelYear( counter, 4 )", TimeRefArraySelYear( counter, 4 )
										print *, "TimeRefArraySelYear( counter, 5 )", TimeRefArraySelYear( counter, 5 )
										
										! find the matching elements of the respecitve year and copy them
										
										ALLOCATE( Time_bnds( 2, counter ) )
										
										print *, "ALLOCATE( Time_bnds( 2, counter ) )", Time_bnds( 2, counter )
										
										counter = 0
										DO i = 1, SIZE(TimeRefArray, 1), 1
											
											IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
												
												counter = counter + 1
												TimeRefArraySelYear(counter,1:5) = TimeRefArray(i,1:5)
												
											END IF
											
										END DO
										
										!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySelYear(:,:) )
										
										!-------------------------------------------------------------------------------
										! check for existance of the file and generate file if needed
										
										INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )
										
										! could exist already from a previous run of tool and due to multiple
										! months in a file (i.e. 1 WRF output file may cover several months)
										IF ( FileExists ) THEN
											
											PRINT *, "path and file exist, continue filling"
																						
											!!!!Change indices for TimeRefArraySelYear (time) and Time_bnds (time bounds) records
											!!!!so that they are written in proper position in same file right after the last year's data
											sts = NF90_OPEN(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid)
											
											!Not sure if all of this is necessary though...
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.), START = (/ counter_prev_year + 1, 1 /) , COUNT = (/ SIZE(TimeRefArraySelYear(:,1), 1) /) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, counter_prev_year+1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF
											!print *,'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
											
											
											!!!!!!!!!!!!!!!!!!!!!!!!
											!!!I'm pretty sure that this is not necessary here...
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											!!!!!!!!!!!!!!!!!!!!!!!!
											
											sts = NF90_CLOSE(ncid)
											
											!!!!End adjustments regarding index w.r.t. Time and Time_bnds when running for consecutive years
										
										ELSE
											
											PRINT *, "path and file do not yet exist, create path and NetCDF file first"
											PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)
											
											CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )
											
											!-------------------------------------------------------------------------------
											! non-std, works for gfortran (fct & subroutine) + ifort
											! comment lines in the NetCDF file global attribute definition
											! turn standard checking in Makefile off
											! trackingID = "xxxxxxxx-xxxx-Mxxx-Nxxx-xxxxxxxxxxxx"
											! creationDate = "YYYY-MM-DD-THH:MM:SSZ"
											
											CALL SYSTEM(cmdUUID)
											OPEN(1,FILE="tmpfileUUID",STATUS='old')
											READ(1,*) trackingID
											CLOSE(1)
											PRINT *, "uuidgen externally generated trackingID = ", trackingID
											
											CALL SYSTEM(cmdDate)
											OPEN(1,FILE="tmpfileDate",STATUS='old')
											READ(1,*) creationDate
											CLOSE(1)
											PRINT *, "date externally generated creation date = ", creationDate
											
											!-------------------------------------------------------------------------------
											! create NetCDF file
											! NF90_CLASSIC_MODEL = NetCDF4_classic
											! NF90_HDF5 = NetCDF4 based on HDF5
											! NF90_CLOBBER = old NetCDF
											! sts = NF90_CREATE(PathFileNameOutTEST, NF90_HDF5, ncid)
											
											!comb_flags = IOR(NF90_HDF5, NF90_CLASSIC_MODEL)
											!https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/NF90_005fCREATE.html
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid)   !not sure whether this is the right data format spec. I guess it may be right using compression but not the other fancy stuff from NetCDF4
											!sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncid)   !if anything, then use this here
											sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_HDF5, ncid)
											
											
											! always included
											sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
											sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid)
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
												
												ELSE
													
													sts = NF90_DEF_DIM(ncid, "plev", 1, height_dimid)
													
												END IF
												
											END IF
											sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
												
											END IF
											
											
											! always included
											sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lon_varid)
											sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "long_name", "longitude")
											sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
											
											! always included
											sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lat_varid)
											sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "long_name", "latitude")
											sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
											
											! always included
											sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid)
											sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
											sts = nf90_put_att(ncid, rlon_varid, "long_name", "longitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
											
											! always included
											sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid)
											sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
											sts = nf90_put_att(ncid, rlat_varid, "long_name", "latitude in rotated pole grid")
											sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
											sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")
											
											! always included
											! restriction to one domain only
											sts = nf90_def_var(ncid, "rotated_pole", NF90_DOUBLE, rotated_pole_varid)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", 39.25)
											sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", -162.0)
											
											! depends whether height is set in the nml
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_def_var(ncid, "height", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
													sts = nf90_put_att(ncid, height_varid, "long_name", "height")
													sts = nf90_put_att(ncid, height_varid, "units", "m")
													sts = nf90_put_att(ncid, height_varid, "positive", "up")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")
												
												ELSE
													
													sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, height_varid)
													sts = nf90_put_att(ncid, height_varid, "standard_name","air_pressure")
													sts = nf90_put_att(ncid, height_varid, "long_name", "pressure")
													sts = nf90_put_att(ncid, height_varid, "units", "Pa")
													sts = nf90_put_att(ncid, height_varid, "positive", "down")
													sts = nf90_put_att(ncid, height_varid, "axis", "Z")                 
													
												END IF
												
											END IF
											
											!missing: lvl, depends
											
											! always included
											sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid)
											sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "long_name", "time")
											sts = nf90_put_att(ncid, rec_varid, "units", "days since 1949-12-01T00:00:00Z")
											sts = nf90_put_att(ncid, rec_varid, "calendar", "standard")
											sts = nf90_put_att(ncid, rec_varid, "axis", "T")
											
											! for mean variables
											IF ( cell_methods(ivar) == "mean" ) THEN
												! Add time bounds to time varaible
												
												sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")	
												sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid)
												sts = nf90_put_att(ncid, recbnds_varid, "standard_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "long_name", "time")
												sts = nf90_put_att(ncid, recbnds_varid, "units", "days since 1949-12-01T00:00:00Z")
												sts = nf90_put_att(ncid, recbnds_varid, "calendar", "standard")
												!sts = nf90_put_att(ncid, recbnds_varid, "axis", "T")
												
												print *,'rec_varid', rec_varid
												print *,'recbnds_varid', recbnds_varid
												
											END IF
											
											! always included
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact", contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "creation_date", creationDate)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment", experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "experiment_id", experiment_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment", driving_experiment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_id", driving_model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_model_ensemble_member", driving_model_ensemble_member)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment_name", driving_experiment_name)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "frequency", frequency(ifrq))
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution", institution)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_id", institute_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "model_id", model_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "forcing",forcing)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact",contact)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source",source)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "rcm_version_id", rcm_version_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "CORDEX_domain", CORDEX_domain)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id", trackingID)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "title",title)
											sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_run_id", institute_run_id)
											
											! always included
											IF ( height(ivar) /= -999 ) THEN
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, height_dimid, rec_dimid /), x_varid)
											
											ELSE
												
												sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid)
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
											sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
											
											IF ( positive(ivar) /= '-999' ) THEN
												
												sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
											
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat height")
												
												ELSE
													
													sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat plev")
													
												END IF
											
											ELSE
												
												sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat")
												
											END IF
											
											sts = nf90_put_att(ncid, x_varid, "grid_mapping", "rotated_pole")
											
											!IF ( height(ivar) == 850 ) THEN
											sts = nf90_put_att(ncid, x_varid, "missing_value", 1.e20)
											sts = nf90_put_att(ncid, x_varid, "_FillValue", 1.e20)
											!END IF
											
											sts = NF90_ENDDEF(ncid)
											
											IF ( CheckForLeapyear( INT(TimeRefArraySelYear(it,2)), dataIncludesLeapYearDays ) == 366 ) THEN
												
												NoOfDaysPerMonth = (/31,29,31,30,31,30,31,31,30,31,30,31/)
											
											ELSE
												
												NoOfDaysPerMonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
												
											END IF
											
											IF (dtHours == 720) THEN
												
												HoursOfCurrentInterval = NoOfDaysPerMonth(TimeRefArraySelYear(it,3)) * 24
											
											ELSE IF (dtHours == 2880) THEN
												
												IF (TimeRefArraySelYear(it,3) == 12 .OR. &
												TimeRefArraySelYear(it,3) == 1 .OR. &
												TimeRefArraySelYear(it,3) == 2) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(1) + NoOfDaysPerMonth(2) + & 
													NoOfDaysPerMonth(12)) * 24 !It's Winter..
													!adding forthcoming december, but, still.. it's 31 days anyway
												
												ELSE IF (TimeRefArraySelYear(it,3) == 3 .OR. &
												TimeRefArraySelYear(it,3) == 4 .OR. &
												TimeRefArraySelYear(it,3) == 5) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(3) + NoOfDaysPerMonth(4) + & 
													NoOfDaysPerMonth(5)) * 24 !It's Spring..
												
												ELSE IF (TimeRefArraySelYear(it,3) == 6 .OR. &
												TimeRefArraySelYear(it,3) == 7 .OR. &
												TimeRefArraySelYear(it,3) == 8) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(6) + NoOfDaysPerMonth(7) + & 
													NoOfDaysPerMonth(8)) * 24 !It's Summer..
												
												ELSE !IF (TimeRefArraySelYear(it,3) == 9 .OR. &
													!TimeRefArraySelYear(it,3) == 10 .OR. &
													!TimeRefArraySelYear(it,3) == 11) THEN
													
													HoursOfCurrentInterval = (NoOfDaysPerMonth(9) + NoOfDaysPerMonth(10) + & 
													NoOfDaysPerMonth(11)) * 24 !Es ist Herbst..
													
												END IF
											
											ELSE
												
												HoursOfCurrentInterval = dtHours
												
											END IF
											
											
											! add time, whole year from above
											IF ( cell_methods(ivar) == "point" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySelYear(:,1) )
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
											END IF
											
											IF ( cell_methods(ivar) == "mean" ) THEN
												
												print*, 'cell_methods:', cell_methods(ivar)
												sts = NF90_PUT_VAR(ncid, rec_varid, (TimeRefArraySelYear(:,1)+(HoursOfCurrentInterval/2.)/24.) )
												
												print *, 'sts NF90_PUT_VAR time', sts
												
												print*, 'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
												print*, 'TimeRefArraySelYear(:,2)', TimeRefArraySelYear(:,2)
												print*, 'TimeRefArraySelYear(:,3)', TimeRefArraySelYear(:,3)
												print*, 'TimeRefArraySelYear(:,4)', TimeRefArraySelYear(:,4)
												print*, 'TimeRefArraySelYear(:,5)', TimeRefArraySelYear(:,5)
												
												Time_bnds(1,:) = TimeRefArraySelYear(:,1)
												Time_bnds(2,:) = TimeRefArraySelYear(:,1)+HoursOfCurrentInterval/24.
												
												print*, 'recbnds_varid', recbnds_varid
												
												sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
												
											END IF
											!print *,'TimeRefArraySelYear(:,1)', TimeRefArraySelYear(:,1)
											
											
											sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
											START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
											START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )
											sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
											sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
											
											! add time_bnds, calc here
											
											! add height from NML
											IF ( height(ivar) /= -999 ) THEN
												
													IF (  ( height(ivar) /= 200) .and. ( height(ivar) /= 850 ) .and. (height(ivar) /= 500 )) THEN !.and. (height(ivar) /= 100 )  .and. (height(ivar) /= 300 ) .and. (height(ivar) /= 400 ) .and. (height(ivar) /= 600 ) .and. (height(ivar) /= 700 ) .and. (height(ivar) /= 800 ) .and. (height(ivar) /= 900 )) THEN
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
												
												ELSE
													
													sts = NF90_PUT_VAR(ncid, height_varid, height(ivar)*100. )
													
												END IF
												
											END IF
											
											sts = NF90_CLOSE(ncid)
											
										END IF ! file exists y/n
										
									END IF ! checking with previous year and month
									
									!-------------------------------------------------------------------------------
									! match timestep of WRFin with the subset of the ref time vec which belong to
									! the NetCDF file of the year currently open to receive data
									! NOT SURE WHETHER THIS IS NEEDED AT ALL, THIS IS THE WRONG DIRECTION ???
									! see whether the current time of the timestep fits anywhere in the
									
									! ARIS CHANGE --------------- This small excerpt below is used to locate where
									! in the final nc file is the data to be written (i.e. in position "counter").
									! However, due to the nature of the TimeRefArraySelYear variable, the original
									! WRF data and the aforementioned variable do not have the same frequency and
									! instants in time as records. This incongruence completely destroys the purpose
									! of computing averages. As a result, this excerpt here must ONLY be executed if
									! the WRF input data and the requested output frequency coincide (i.e. no temporal)
									! averaging takes place.
									
									IF (record_count.eq.1) THEN
										
										PRINT *, "reading WRF sim. res. = ", TRIM(InVarDataRec(it)), it
										!PRINT *, SIZE(TimeRefArraySelYear,1)
										!PRINT *, SHAPE(TimeRefArraySelYear)
										!PRINT *, "current transferred input time: ", InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
										
										counter = 0
										DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
											
											counter = counter + 1
											
											!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
											
											IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(it)  ) .AND. &
											( TimeRefArraySelYear(i,3) == InDateTimeMonth(it) ) .AND. &
											( TimeRefArraySelYear(i,4) == InDateTimeDay(it)   ) .AND. &
											( TimeRefArraySelYear(i,5) == InDateTimeHour(it)  ) ) THEN
												
												EXIT
												
											END IF
											
										END DO
										
										PRINT *, "index where in the NC file the WRF data is sorted in = ", &
										counter
										!SUPER IMPORTANT ADDITION										
										counter_array(it) = counter
									END IF
									! If this is not run here, then it must be run directly after the DO it = ... loop
									! is finished (somewhere in line 1900+)...
									
									
									!-------------------------------------------------------------------------------
									! read orig WRF outpouts
									! there is always a corresponding time-slot in the NC file
									! extracted time from above
									! "it" controls it all: timestep in the individual WRF file
									
									sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)
									
									IF ( (var_cmip(ivar) == "psl") .or. (height(ivar) == 850) &
									.or.(height(ivar) == 500) .or. (height(ivar) == 200) &
									.or. (var_cmip(ivar) == "prw") .or. (var_cmip(ivar) == "clwvi") &
									.or. (var_cmip(ivar) == "clivi") & .or. (var_cmip(ivar) == "cape")) THEN
									! .or.(height(ivar) == 100) .or. (height(ivar) == 300) &
									! .or.(height(ivar) == 400) .or. (height(ivar) == 600) &
									! .or.(height(ivar) == 700) .or. (height(ivar) == 800) &
									! .or.(height(ivar) == 900)) THEN
										
										! SKn: It is not necessary to read all 3D variables for every single output variable.
										!      Here it is done to have a more compact structure, but it could be separated 
										!      in multiple if-blocks for every variable.
										
										ALLOCATE( pp_in( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( pb_in( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( ph_in( xfocus, yfocus, nz+1 ), STAT=sts )
										ALLOCATE( phb_in( xfocus, yfocus, nz+1 ), STAT=sts )
										ALLOCATE( theta_in( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( qv_in( xfocus, yfocus, nz  ), STAT=sts )
										ALLOCATE( qc_in( xfocus, yfocus, nz  ), STAT=sts )
										ALLOCATE( qi_in( xfocus, yfocus, nz  ), STAT=sts )
										ALLOCATE( qr_in( xfocus, yfocus, nz  ), STAT=sts )
										ALLOCATE( qs_in( xfocus, yfocus, nz  ), STAT=sts )
										
										ALLOCATE( t_in( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( ph_fl( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( u_in( xfocus+1, yfocus, nz ), STAT=sts )
										ALLOCATE( v_in( xfocus, yfocus+1, nz ), STAT=sts )
										ALLOCATE( var3d_in( xfocus, yfocus, nz ), STAT=sts )
										
										ALLOCATE( psl_in ( xfocus, yfocus ), STAT=sts )
										ALLOCATE( t2_in ( xfocus, yfocus ), STAT=sts )          
										
										
										ALLOCATE( t_p( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( qvs( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( cape( xfocus, yfocus ), STAT=sts )
										ALLOCATE( cin( xfocus, yfocus ), STAT=sts )
										ALLOCATE( lcl( xfocus, yfocus ), STAT=sts )
										ALLOCATE( lfc( xfocus, yfocus ), STAT=sts )
										
										ALLOCATE( prw( xfocus, yfocus ), STAT=sts )
										ALLOCATE( clwvi( xfocus, yfocus ), STAT=sts )
										ALLOCATE( clivi( xfocus, yfocus ), STAT=sts )
										
										ALLOCATE( p_in( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( pp_in( xfocus, yfocus, nz ), STAT=sts )
										ALLOCATE( var_pl( xfocus, yfocus, 3 ), STAT=sts )
										!Change dimension due to experimental runs of 'ta' variable
										!ALLOCATE( var_pl( xfocus, yfocus, 9 ), STAT=sts )
										
										ALLOCATE( pout( 3 ), STAT=sts )
										!Change dimension due to experimental runs of 'ta' variable
										!ALLOCATE( pout( 9 ), STAT=sts ) 
										
										ALLOCATE( sinalpha_in( xfocus, yfocus ), STAT=sts )
										ALLOCATE( cosalpha_in( xfocus, yfocus ), STAT=sts )
										
										
										print *,'read 3D vars'
										
										sts = NF90_INQ_VARID(ncidin, "P", pp_varid)
										sts = NF90_INQ_VARID(ncidin, "PB", pb_varid)
										sts = NF90_INQ_VARID(ncidin, "PH", ph_varid)
										sts = NF90_INQ_VARID(ncidin, "PHB", phb_varid)
										sts = NF90_INQ_VARID(ncidin, "T", theta_varid)
										sts = NF90_INQ_VARID(ncidin, "QVAPOR", qv_varid)
										sts = NF90_INQ_VARID(ncidin, "QCLOUD", qc_varid)
										sts = NF90_INQ_VARID(ncidin, "QICE", qi_varid)
										sts = NF90_INQ_VARID(ncidin, "QRAIN", qr_varid)
										sts = NF90_INQ_VARID(ncidin, "QSNOW", qs_varid)
										sts = NF90_INQ_VARID(ncidin, "U", u_varid)
										sts = NF90_INQ_VARID(ncidin, "V", v_varid)
										sts = NF90_INQ_VARID(ncidin, "SINALPHA", sinalpha_varid)
										sts = NF90_INQ_VARID(ncidin, "COSALPHA", cosalpha_varid)
										
										
										sts = NF90_INQ_VARID(ncidin, "T2", t2_varid)
										
										sts = NF90_GET_VAR(ncidin, t2_varid, t2_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, pp_varid, pp_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )          
										
										sts = NF90_GET_VAR(ncidin, pb_varid, pb_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, ph_varid, ph_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz+1, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, phb_varid, phb_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz+1, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, theta_varid, theta_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, qv_varid, qv_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, qc_varid, qc_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, qi_varid, qi_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, qr_varid, qr_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, qs_varid, qs_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )          
										
										sts = NF90_GET_VAR(ncidin, u_varid, u_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus+1, yfocus, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, v_varid, v_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus+1, nz, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, sinalpha_varid, sinalpha_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, cosalpha_varid, cosalpha_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										
										
										
										
									
									ELSE IF (var_cmip(ivar) == "clt") THEN
										
										ALLOCATE( cldfra_in( xfocus, yfocus, nz ), STAT=sts )     
										
										sts = NF90_INQ_VARID(ncidin, "CLDFRA", varid)
										
										sts = NF90_GET_VAR(ncidin, varid, cldfra_in(:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
										
									
									ELSE IF (var_cmip(ivar) == "pr") THEN 
										
										ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )   !read two timesteps to calculate 3hr sum
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prc") THEN
										
										ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
										
										sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
									
									ELSE IF (var_cmip(ivar) == "prsn") THEN
										
										ALLOCATE( snownc_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
										
										sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
									
									ELSE IF (var_cmip(ivar) == "snm") THEN
										
										ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
										
										sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
										print*, sts
									
									ELSE IF (var_cmip(ivar) == "evspsbl") THEN
										
										ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
										
										sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
									
									ELSE IF (var_cmip(ivar) == "evspsblpot") THEN
										
										ALLOCATE( potevp_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
										
										sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
									
									ELSE IF (var_cmip(ivar) == "mrros") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
									
									ELSE IF (var_cmip(ivar) == "mrro") THEN
										
										ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
										sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
										
										sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
										sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
										
									
									ELSE IF ((var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")  &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "rlut")                                 &
									.or. (var_cmip(ivar) == "rsdt") .or. (var_cmip(ivar) == "rsut") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										ALLOCATE( rad_in ( xfocus, yfocus, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										
										sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
										
										
										print*, var_cmip(ivar), rad_in(50,50,1), rad_in(50,50,2)
										print*, 'difference in J m-2', (rad_in(50,50,2) - rad_in(50,50,1))
										print*, 'in mean W m-2', (rad_in(50,50,2) - rad_in(50,50,1))/ (obs_interval*3600.)
										
										! alternative: since accumulated values as read above get so large in 
										! long term simulations that their differences loose accuracy, use 
										! instantaneous values instead and calculate means
										
										
										
									
									ELSE IF (var_cmip(ivar) == "mrso") THEN
										
										ALLOCATE( smois_in( xfocus, yfocus, 4, 2 ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "SMOIS", varid)
										
										sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,:), &
										START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, 4, 2 /) )
										
									
									ELSE IF (var_cmip(ivar) == "sfcWind") THEN
										
										ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts ) 
										ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
										
										sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)
										
										sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									
									ELSE IF ((var_cmip(ivar) == "uas") .or. (var_cmip(ivar) == "vas")) THEN
										
										ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts )
										ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts )
										ALLOCATE( sinalpha_in( xfocus, yfocus ), STAT=sts )
										ALLOCATE( cosalpha_in( xfocus, yfocus ), STAT=sts )
										
										sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
										sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)
										sts = NF90_INQ_VARID(ncidin, "SINALPHA", sinalpha_varid)
										sts = NF90_INQ_VARID(ncidin, "COSALPHA", cosalpha_varid)
										
										sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, sinalpha_varid, sinalpha_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										
										sts = NF90_GET_VAR(ncidin, cosalpha_varid, cosalpha_in(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )          
										
									
									ELSE
										
										ALLOCATE( temp_data ( xfocus, yfocus ), STAT=sts )
										sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
										!print *, "VALUE FOR NAME: ", TRIM(var_wrf(ivar)), "+ VARID: ", varid
										
										
										sts = NF90_GET_VAR(ncidin, varid, temp_data(:,:), &
										START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
										!			data_in(:,:) = SUM(temp_data, 3)/8.
										!			print *, temp_data(50,50,1), temp_data(50,50,2), temp_data(50,50,3), temp_data(50,50,4)
										!			print *, temp_data(50,50,5), temp_data(50,50,6), temp_data(50,50,7), temp_data(50,50,8)
										!			print *, "Mean of data_in(50,50):", data_in(50,50)
									END IF
									
									sts = NF90_CLOSE(ncidin)
									
									!-------------------------------------------------------------------------------
									! some analysis of the data
									
									print *,"shape of array" , SHAPE(temp_data)
									print *,"size of array" , SIZE(temp_data)
									!
									!       stat_mean = SUM(temp_data(:,:,5))/(MAX(1,SIZE(temp_data(:,:,5))))
									!       PRINT *, stat_mean
									stat_mean = SUM(temp_data(:,:))/SIZE(temp_data(:,:))
									PRINT *, stat_mean
									
									!-------------------------------------------------------------------------------
									! processing
									
									!       ***psl***   ***vars on pressure levels***     
									
									IF ( (var_cmip(ivar) == "psl") .or. (height(ivar) == 850) &
									.or.(height(ivar) == 500) .or. (height(ivar) == 200) ) THEN !&
									! .or.(height(ivar) == 100) .or. (height(ivar) == 300) &
									! .or.(height(ivar) == 400) .or. (height(ivar) == 600) &
									! .or.(height(ivar) == 700) .or. (height(ivar) == 800) &
									! .or.(height(ivar) == 900) ) THEN
										
										DO nl = 1,nz
											!print*, "ph_in(:,:,", nl, "):", ph_in(:,:,nl)
											!print*, "phb_in(:,:,", nl, "):", phb_in(:,:,nl)
											
											ph_fl(:,:,nl) = ((ph_in(:,:,nl)+phb_in(:,:,nl))+(ph_in(:,:,nl+1)+phb_in(:,:,nl+1)))/2./9.81
											!print*, "ph_fl(:,:,", nl, "):", ph_fl(:,:,nl)
										END DO
										
										!print*, "pp_in:", pp_in
										!print*, "pb_in:", pb_in
										!print*, "theta_in:", theta_in(:,:,:)
										t_in(:,:,:) = (theta_in(:,:,:)+300.)*((pp_in(:,:,:)+pb_in(:,:,:))/100000.)**(287./1004.)
										!print*, "t_in:", t_in
										
										IF  (var_cmip(ivar) == "psl") THEN
											
											psl_in(:,:) = (pp_in(:,:,1)+pb_in(:,:,1))*((t_in(:,:,1)*(1.+0.61*qv_in(:,:,1))+0.0065*ph_fl(:,:,1))/(t_in(:,:,1)*(1+0.61*qv_in(:,:,1))))**(9.81/(287.*0.0065))
											
											temp_data(:,:)=psl_in(:,:)
										
										ELSE
											
											pout = (/ 85000.,50000.,20000. /)
											!Change for experimenting with 'ta' vars
											!pout = (/ 50000.,20000.,10000.,30000.,40000.,60000.,70000.,80000.,90000. /)
											

											p_in = pp_in+pb_in
											!print*, "p_in:", p_in
											
											
											!DO np = 1,3     !SKn: could loop over heigts per variable or calculate t850, t500, t200 as individual variables 
											
											IF (height(ivar) == 850) THEN
												
												np = 1
											
											ELSE IF (height(ivar) == 500) THEN
												
												np = 2
											
											ELSE IF (height(ivar) == 200) THEN
												
												np = 3
												
											! ELSE IF (height(ivar) == 100) THEN
												
												! np = 3
												
											! ELSE IF (height(ivar) == 300) THEN
												
												! np = 4
												
											! ELSE IF (height(ivar) == 400) THEN
												
												! np = 5

											! ELSE IF (height(ivar) == 600) THEN
												
												! np = 6

											! ELSE IF (height(ivar) == 700) THEN
												
												! np = 7

											! ELSE IF (height(ivar) == 800) THEN
												
												! np = 8

											! ELSE IF (height(ivar) == 900) THEN
												
												! np = 9	
												
											END IF
											
											!print *,'np', np      
											
											IF ( (var_cmip(ivar) == "ta850") .or. (var_cmip(ivar) == "ta500") .or. (var_cmip(ivar) == "ta200") ) THEN !&
													! .or. (var_cmip(ivar) == "ta100") .or. (var_cmip(ivar) == "ta300") .or. (var_cmip(ivar) == "ta400") &
													! .or. (var_cmip(ivar) == "ta600") .or. (var_cmip(ivar) == "ta700") .or. (var_cmip(ivar) == "ta800") &
													! .or. (var_cmip(ivar) == "ta900")	) THEN
												var3d_in(:,:,:) = t_in(:,:,:)
											!print*, 'var3d_in(1,1,1)', var3d_in(1,1,1), var_cmip(ivar)
											ELSE IF ( (var_cmip(ivar) == "hus850") ) THEN
												var3d_in(:,:,:) = qv_in(:,:,:)
											!print*, 'var3d_in(50,50,10)', var3d_in(50,50,10), var_cmip(ivar)
											ELSE IF ( (var_cmip(ivar) == "ua850") .or. (var_cmip(ivar) == "ua500") .or. (var_cmip(ivar) == "ua200") ) THEN
												DO i = 1,xfocus
													var3d_in(i,:,:) = (u_in(i,:,:)+u_in(i+1,:,:))/2.*cosalpha_in(:,:) - (v_in(i,:,:)+v_in(i+1,:,:))/2.*sinalpha_in(:,:) !rotate to earth grid
												END DO
											!print*, 'var3d_in(50,50,10)', var3d_in(50,50,10), var_cmip(ivar)
											ELSE IF ( (var_cmip(ivar) == "va850") .or. (var_cmip(ivar) == "va500") .or. (var_cmip(ivar) == "va200") ) THEN
												DO j = 1,yfocus
													var3d_in(:,j,:) = (v_in(:,j,:)+v_in(:,j+1,:))/2.*cosalpha_in(:,:) + (u_in(i,:,:)+u_in(i+1,:,:))/2.*sinalpha_in(:,:) !rotate to earth grid
												END DO
											!print*, 'var3d_in(50,50,10)', var3d_in(50,50,10), var_cmip(ivar)
											ELSE IF ( (var_cmip(ivar) == "zg500") .or. (var_cmip(ivar) == "zg200") ) THEN
												var3d_in(:,:,:) = ph_fl(:,:,:)
												!print*, 'var3d_in(50,50,10)', var3d_in(50,50,10), var_cmip(ivar)
											END IF
											
											var_pl = 1.e20
											
											DO i = 1,xfocus 
												DO j = 1,yfocus
													DO nl = 1,nz - 1
														IF (pout(np).lt.p_in(i,j,nl) .and. pout(np).gt.p_in(i,j,nl+1)) then
															
															slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ (p_in(i,j,nl)-p_in(i,j,nl+1))
															!print*, "pout(np)", np, pout(np)
															!print*, "var3d_in(i,j,nl+1)", var3d_in(i,j,nl+1)
															var_pl(i,j,np) = var3d_in(i,j,nl+1) + slope* (pout(np)-p_in(i,j,nl+1))
															
														END IF
													END DO
												END DO
											END DO
											!END DO
											
											temp_data(:,:) = var_pl(:,:,np)
											
										END IF
										
									END IF
									
									
									!       ***prw, clwvi, clivi***
									
									IF ( (var_cmip(ivar) == "prw") ) THEN     
										
										t_in(:,:,:) = (theta_in(:,:,:)+300.)*((pp_in(:,:,:)+pb_in(:,:,:))/100000.)**(R/cp)
										p_in = pp_in+pb_in
										
										prw(:,:) = 0.
										
										DO nl = 1,nz
											
											prw(:,:) = prw(:,:) + qv_in(:,:,nl) * p_in(:,:,nl)/(R*t_in(:,:,nl)) * ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+phb_in(:,:,nl)))/9.81
											
											temp_data(:,:) = prw(:,:)            
											
										END DO
										
									END IF
									
									
									IF ( (var_cmip(ivar) == "clwvi") ) THEN
										
										t_in(:,:,:) = (theta_in(:,:,:)+300.)*((pp_in(:,:,:)+pb_in(:,:,:))/100000.)**(R/cp)
										p_in = pp_in+pb_in          
										
										clwvi(:,:) = 0.
										
										DO nl = 1,nz
											
											clwvi(:,:) = clwvi(:,:) + (qc_in(:,:,nl) + qi_in(:,:,nl) + qr_in(:,:,nl) + qs_in(:,:,nl) ) * p_in(:,:,nl)/(R*t_in(:,:,nl)) * ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+phb_in(:,:,nl)))/9.81
											
										END DO
										
										temp_data(:,:) = clwvi(:,:) 
										
									END IF
									
									
									IF ( (var_cmip(ivar) == "clivi")) THEN
										
										t_in(:,:,:) = (theta_in(:,:,:)+300.)*((pp_in(:,:,:)+pb_in(:,:,:))/100000.)**(R/cp)
										p_in = pp_in+pb_in
										
										clivi(:,:) = 0.
										
										DO nl = 1,nz
											
											clivi(:,:) = clivi(:,:) + (qi_in(:,:,nl) + qs_in(:,:,nl)) * p_in(:,:,nl)/(R*t_in(:,:,nl)) * ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+phb_in(:,:,nl)))/9.81
											
										END DO
										
										temp_data(:,:) = clivi(:,:)
										
									END IF
									
									
									!       ***cape***
									
									IF ( (var_cmip(ivar) == "cape") ) THEN
										
										t_in(:,:,:) = (theta_in(:,:,:)+300.)*((pp_in(:,:,:)+pb_in(:,:,:))/100000.)**(287./1004.)
										
										p_in = pp_in+pb_in          
										
										t_p(:,:,1) = t_in(:,:,1)
										
										cape(:,:) = 0.
										cin(:,:) = 0.
										lcl(:,:) = -999.
										lfc(:,:) = -999.
										
										DO i = 1,xfocus
											DO j = 1,yfocus
												DO nl = 1,nz-1
													
													qvs(i,j,nl) = 0.622*a*exp(b*(t_p(i,j,nl)-c)/(t_p(i,j,nl)-d))/p_in(i,j,nl)
													
													IF (qvs(i,j,nl) .gt. qv_in(i,j,1)) THEN !dry adiabatic ascent
														
														t_p(i,j,nl+1) = (theta_in(i,j,1)+300.)*(p_in(i,j,nl+1)/100000.)**(R/cp)   
													
													ELSE IF (qvs(i,j,nl) .lt. qv_in(i,j,1)) THEN ! moist adiabatic ascent
														
														IF (lcl(i,j) .eq. -999) THEN    ! lifting condensation level
															lcl(i,j) = p_in(i,j,nl)
														END IF
														
														t_ii = t_p(i,j,nl)
														
														DO ii = 1,10  !solve iteratively
															
															qvs(i,j,nl+1) = 0.622*a*exp(b*(t_ii-c)/(t_ii-d))/p_in(i,j,nl+1)
															
															t_ii = t_ii - (t_ii*(100000./p_in(i,j,nl+1))**(R/cp)*exp(L*qvs(i,j,nl+1)/(cp*t_ii)) &
															- (t_p(i,j,nl)*(100000./p_in(i,j,nl))**(R/cp)*exp(L*qvs(i,j,nl)/(cp*t_p(i,j,nl))))) &
															/ ( (100000./p_in(i,j,nl+1))**(R/cp)*exp(n/(p_in(i,j,nl+1)*t_ii)*exp(b*(t_ii-c)/(t_ii-d))) * &
															(1 - (n/p_in(i,j,nl+1)*exp(b*(t_ii-c)/(t_ii-d))*(t_ii*(t_ii-b*c)+(b-2)*d*t_ii+d**2))/(t_ii*(d-t_ii)**2)) )
															
														END DO
														
														!print*, 'thetae(i,j,nl)',(t_p(i,j,nl)*(100000./p_in(i,j,nl))**(R/cp)*exp(L*qvs(i,j,nl)/(cp*t_p(i,j,nl))))
														!print*,'thetae(i.j.nl+1)',(t_ii*(100000./p_in(i,j,nl+1))**(R/cp)*exp(L*qvs(i,j,nl+1)/(cp*t_ii))) 
														
														
														t_p(i,j,nl+1) = t_ii
														
														!print*, nl, 'moist', t_p(i,j,nl+1), t_in(i,j,nl+1), (t_p(i,j,nl+1)-t_in(i,j,nl+1))
														
													END IF                 
													
													
													IF (t_p(i,j,nl) .gt. t_in(i,j,nl)) THEN
														
														IF (lfc(i,j) .eq. -999) THEN   ! level of free convection
															lfc(i,j) = p_in(i,j,nl)
														END IF
														
														cape(i,j) = cape(i,j) + (t_p(i,j,nl) - t_in(i,j,nl)) / t_in(i,j,nl) * ((phb_in(i,j,nl)+ph_in(i,j,nl))-(phb_in(i,j,nl-1)+ph_in(i,j,nl-1)))
														
														!print*, 'nl, cape(i,j)', nl, cape(i,j)
													
													ELSE IF ( (t_p(i,j,nl) .lt. t_in(i,j,nl)) .and. (cape(i,j) .eq. 0.) )  THEN   !convective inhibition 
														
														cin(i,j) = cin(i,j) + (t_in(i,j,nl) - t_p(i,j,nl)) / t_in(i,j,nl) * ((phb_in(i,j,nl)+ph_in(i,j,nl))-(phb_in(i,j,nl-1)+ph_in(i,j,nl-1))) 
														
														!print*, 'nl, cin(i,j)', nl, cin(i,j)
														
													END IF
													
													
												END DO
											END DO
										END DO
										
										temp_data(:,:) = cape(:,:)
										
									END IF
									
									!       ***clt***
									
									IF (var_cmip(ivar) == "clt") THEN
										
										ALLOCATE( cldfra_inv( xfocus, yfocus ), STAT=sts )
										
										cldfra_inv(:,:) = 1.
										
										!DO nl = 1,40 - 1
										DO i = 1,xfocus
											DO j = 1,yfocus
												IF (maxval(cldfra_in(i,j,:)) .lt. 0.99) THEN
													
													!print *, "maxval(cldfra_in(",i,",",j,",:))", maxval(cldfra_in(i,j,:))
													cldfra_inv(i,j) = 1.
													
													DO nl = 2,nz
														
														cldfra_inv(i,j) = cldfra_inv(i,j)*(1- max(cldfra_in(i,j,nl),cldfra_in(i,j,nl-1))/(1-cldfra_in(i,j,nl-1))) !unit [%] 
														!print *, "Checking final clt value:", cldfra_inv(i,j)
													
													END DO
												ELSE 
													!print *, "ELSE condition met (cldfra_inv(i,j) = 0)"
													cldfra_inv(i,j) = 0.  
														
												END IF
											END DO
										END DO
										!END DO
										
										temp_data(:,:) = (1 - cldfra_inv(:,:))*100.
										
										
										WHERE (temp_data .gt. 100.) temp_data = 100.
										WHERE (temp_data .lt. 0.) temp_data = 0.
										!print *, "Final temp_data for clt:", temp_data(:,:)
										
									END IF
									
									
									!       ***pr***
									IF (var_cmip(ivar) == "pr") THEN 
										
										temp_data(:,:) = ((rainnc_in(:,:,2) + rainc_in(:,:,2)) - (rainnc_in(:,:,1) + rainc_in(:,:,1)))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										!ATTENTION: implement adjustable time intervals that the differences are devided by
										
									END IF
									
									
									!       ***prc***
									IF (var_cmip(ivar) == "prc") THEN
										
										temp_data(:,:) = (rainc_in(:,:,2) - rainc_in(:,:,1))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***prsn***
									IF (var_cmip(ivar) == "prsn") THEN
										
										temp_data(:,:) = (snownc_in(:,:,2) - snownc_in(:,:,1))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***snm***
									IF (var_cmip(ivar) == "snm") THEN
										
										temp_data(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1))/(obs_interval*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspsbl***
									IF (var_cmip(ivar) == "evspsbl") THEN
										
										temp_data(:,:) = (sfcevp_in(:,:,2) - sfcevp_in(:,:,1))/(obs_interval*3600.) !unit [kg m-2 /3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***evspblpot**
									IF (var_cmip(ivar) == "evspsblpot") THEN
										
										temp_data(:,:) = (potevp_in(:,:,2) - (potevp_in(:,:,1)))/L   !unit [W m-2]/[J kg-1] -> [kg m-2 s-1]
										
										! THERE IS STH WRONG WITH THE UNITS: WRF's POTEVP is accumulated and declared to be in W m-2. 
										! It doesen't make sense to accumulate in W m-2, but even if assume it as J m-2 or derive kg m-2 
										! by using latent heat of vaporization you never get values that have a reasonable magnitude...
										
										
									END IF
									
									!       ***mrros***
									IF (var_cmip(ivar) == "mrros") THEN
										
										temp_data(:,:) = (sfroff_in(:,:,2) - sfroff_in(:,:,1))/(obs_interval*3600.)       !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									!       ***mrro***
									IF (var_cmip(ivar) == "mrro") THEN
										
										temp_data(:,:) = ((sfroff_in(:,:,2) - sfroff_in(:,:,1)) + (udroff_in(:,:,2) - udroff_in(:,:,1)))/(obs_interval*3600.) !unit [mm/3hr] to [kg m-2 s-1]
										
									END IF
									
									
									!       ***snc,sic***
									IF ( (var_cmip(ivar) == "snc") .or. (var_cmip(ivar) == "sic") ) THEN
										
										temp_data(:,:) = temp_data(:,:)*100. !unit [] to [%]
										
									END IF
									
									
									!       ***rsds, rlds, rsus, rlus***
									IF ( (var_cmip(ivar) == "rsds") .or. (var_cmip(ivar) == "rlds")      &
									.or. (var_cmip(ivar) == "rsus") .or. (var_cmip(ivar) == "rlus") &
									.or. (var_cmip(ivar) == "hfss") .or. (var_cmip(ivar) == "hfls")) THEN
										
										IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) - rad_in(:,:,1)) /(obs_interval*3600.)       ! take difference of accumulated values
										
										ELSE IF (TRIM(fnNMLvar(varnml)) == "runctrl.vars.nml_radiation_alternative") THEN
											
											temp_data(:,:) = (rad_in(:,:,2) + rad_in(:,:,1)) / 2.              ! take mean of instantaneous values
											
										END IF 
										
									END IF
									
									
									!       ***mrso***
									IF (var_cmip(ivar) == "mrso") THEN
										
										temp_data(:,:) = (((smois_in(:,:,1,1)*0.1 + smois_in(:,:,2,1)*0.3 + smois_in(:,:,3,1)*0.6 + smois_in(:,:,4,1)*1.0 ) + &
										(smois_in(:,:,1,2)*0.1 + smois_in(:,:,2,2)*0.3 + smois_in(:,:,3,2)*0.6 + smois_in(:,:,4,2)*1.0 ))/2.)*1000. 
										
									END IF
									
									!       ***sfcWind***
									IF (var_cmip(ivar) == "sfcWind") THEN
										
										temp_data(:,:) = (u10_in(:,:)**2 + v10_in(:,:)**2)**0.5 
										
									END IF
									
									
									!       ***uas***
									IF (var_cmip(ivar) == "uas") THEN 
										
										temp_data(:,:) = u10_in(:,:)*cosalpha_in(:,:) - v10_in(:,:)*sinalpha_in(:,:) ! rotate to earth grid
										
									END IF
									
									IF (var_cmip(ivar) == "vas")  THEN
										
										temp_data(:,:) = v10_in(:,:)*cosalpha_in(:,:) + u10_in(:,:)*sinalpha_in(:,:) ! rotate to earth grid
										
									END IF
									
									
									!Up to this point, temp_data has played the role of data_in... for 3-hr intervals (strictly!)
									
									temp_data_InTime(:,:,it) = temp_data(:,:)
									
								END DO !DO it = 1, ...
								!variable "it" will be InDimLenRec+1 at this point.
								
								
								
								! Perform averaging of data --> This is for one day if data are DAILY!
								temp_data(:,:) = SUM(temp_data_InTime,3)/SIZE(temp_data_InTime,3)
								
								! will not average if requested interval matches the observation interval... (i.e. division by 1)
								! This is always going to be a daily average because the original data are in 1-day files!
								print *, "it = ", it
								currentDayOfMonth = InDateTimeDay(it-1) !hopefully, daily files will be used for input...
								currentMonthOfYear = InDateTimeMonth(it-1) !...and so, the last record will be on the same single day..
								print *, "currentDayOfMonth = ", currentDayOfMonth
								
								IF (dtHours == 720) THEN !month average requested
									
									temp_data_Monthly(:,:,currentDayOfMonth) = temp_data !Add daily average
									
									IF (currentDayOfMonth == NoOfDaysPerMonth(currentMonthOfYear)) THEN !last day of month JUST passed..
										
										isEndOfMonth = 1
										!print *, "End of Month:", currentMonthOfYear			
										temp_data_Month_Average = SUM(temp_data_Monthly,3)/SIZE(temp_data_Monthly,3)
										
										print *, "How many days did the month have??? ", SIZE(temp_data_Monthly,3)
										!print *, "Month average: ", temp_data_Month_Average(1,2)
									
									ELSE
										
										isEndOfMonth = 0
										
									END IF
									
								
								ELSE IF (dtHours == 2880) THEN !season average requested
									
									IF (currentMonthOfYear == 12 .OR. &
									currentMonthOfYear == 3 .OR. &
									currentMonthOfYear == 6 .OR. &
									currentMonthOfYear == 9) THEN !it is the first month of the season
										
										monthIndex = 1
									
									ELSE IF (currentMonthOfYear == 1 .OR. &
									currentMonthOfYear == 4 .OR. &
									currentMonthOfYear == 7 .OR. &
									currentMonthOfYear == 10) THEN !it is the 2nd month of the season
										
										monthIndex = 2 
									
									ELSE !it is the last month of the season..
										
										monthIndex = 3 
										
									END IF
									
									!temp_data is a record of the daily data average. So one day has passed.
									temp_data_Seasonally(:,:,monthIndex,currentDayOfMonth) = temp_data
									dayCount = dayCount + 1
									print *, "the daily average is stored at (:,:,", monthIndex, currentDayOfMonth, ")"
									print *, "the daily average value was at (50,50): ", temp_data_Seasonally(50,50,monthIndex,currentDayOfMonth)
									print *, "temp_data(50,50) was: ", temp_data(50,50)
									
									IF (currentDayOfMonth == NoOfDaysPerMonth(currentMonthOfYear) .AND. &
									(currentMonthOfYear == 2 .OR. currentMonthOfYear == 5 .OR. currentMonthOfYear == 8 .OR. &
									currentMonthOfYear == 11) ) THEN !last months of season AND last day of month just passed
										
										isEndOfSeason = 1								
										
										!yearDays = CheckForLeapyear(InDateTimeYear(it-1), dataIncludesLeapYearDays)
										!IF (currentMonthOfYear == 2) THEN
										!	IF (yearDays == 365) THEN
										!		dayCount = 90
										!	ELSE !it's a leap year
										!		dayCount = 91
										!	END IF
										!	ELSE IF (currentMonthOfYear == 5 .OR. currentMonthOfYear == 8) THEN
										!		dayCount = 92
										!	ELSE !currentMonthOfYear == 11
										!		dayCount = 91
										!	END IF
										
										temp_data_Season_Average = (SUM(temp_data_Seasonally(:,:,1,:),3) + &
										SUM(temp_data_Seasonally(:,:,2,:),3) + &
										SUM(temp_data_Seasonally(:,:,3,:),3) ) / dayCount
										!(SIZE(temp_data_Seasonally(:,:,1,:),3) + 
										!SIZE(temp_data_Seasonally(:,:,2,:),3) + 
										!SIZE(temp_data_Seasonally(:,:,3,:),3) )
										!could not use size because not all months have 31 days and the array was allocated
										!with 31 days... the zeros in element 31 (e.g. for months such as November) won't pose
										!a problem, but the incorrect denominator will... even if it is going to be small...
										print *, "(1,1) of THE SEASON AVERAGE IS EQUAL TO ", temp_data_Season_Average(1,1)
										print *, "dayCount", dayCount
										print *, "tempDataSeasonally(1,1,1,:)", temp_data_Seasonally(1,1,1,:)
										print *, "tempDataSeasonally(1,1,2,:)", temp_data_Seasonally(1,1,2,:)
										print *, "tempDataSeasonally(1,1,3,:)", temp_data_Seasonally(1,1,3,:)
									
									ELSE
										
										isEndOfSeason = 0
										
									END IF
									print *, "is End of season??? ", isEndOfSeason, " and the days included are: ", dayCount
								END IF
								
								
								IF (record_count /= 1) THEN !determine the position of the record (counter)
									
									!PRINT *, "reading WRF sim. res. = ", TRIM(InVarDataRec(it)), it
									!PRINT *, SIZE(TimeRefArraySelYear,1)
									!PRINT *, SHAPE(TimeRefArraySelYear)
									!PRINT *, "current transferred input time: ", InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
									
									counter = 0
									DO i = 1, SIZE(TimeRefArraySelYear,1),1 ! time content of the WRF file
										
										counter = counter + 1
										!PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySelYear(i,1),TimeRefArraySelYear(i,2),TimeRefArraySelYear(i,3),TimeRefArraySelYear(i,4),TimeRefArraySelYear(i,5)
										
										IF (dtHours.eq.24) THEN !we are averaging per day.. so compare DAYS
											IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
											( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) .AND. &
											( TimeRefArraySelYear(i,4) == InDateTimeDay(1)   ) ) THEN
												
												EXIT
												
											END IF
										
										ELSE IF (dtHours.eq.720) THEN !we are averaging per month... so compare Months
											IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
											( TimeRefArraySelYear(i,3) == InDateTimeMonth(1) ) ) THEN
												
												EXIT
												
											END IF
											ELSE IF (dtHours.eq.2880) THEN !we are averaging over seasons.. so ..................
												
												!Setting counter manually to zero at this point, because when computing seasonal avg. then
												!counter_prev_year should come into play in end of last season of year,
												!and not in end of year (in contrast to daily/monthly avg.).
												!Otherwise counter will be added twice right after this loop.
												!I believe this creates no confusion in start of new year (line 5635).
												!Connected to line 7138...
												IF (currentMonthOfYear == 12) THEN
																							
													counter = 1
													EXIT
												
												ELSE IF ( ( TimeRefArraySelYear(i,2) == InDateTimeYear(1)  ) .AND. &
												( TimeRefArraySelYear(i,3) > InDateTimeMonth(1)  ) ) THEN
													
													EXIT
													
												END IF
												
											END IF
											
									END DO
									
									print*, "counter:", counter
									print*, "counter_prev_year:", counter_prev_year
									counter = counter + counter_prev_year
									PRINT *, "index where in the NC file the WRF data is sorted in = ", &
									counter
									
									!Adjusting counter here specifically, because when computing seasonal avg. then
									!previous counter should be stored and come into play in end of last season of year,
									!and not in end of year (in contrast to daily/monthly avg.).
									IF (currentDayOfMonth == NoOfDaysPerMonth(currentMonthOfYear) .AND. currentMonthOfYear == 11) THEN 
										
										counter_prev_year = counter
									
									END IF	
					
					
					
								END IF
								
								IF (record_count.eq.1) THEN
									
									write_iter_count = InDimLenRec
								
								ELSE !we are averaging, so write only once, the requested average!
									IF (isEndOfMonth .OR. isEndOfSeason .OR. dtHours.eq.24) THEN
										
										IF (isEndOfSeason) THEN	
										
											IF (dayCount < 90) THEN !worst case scenario is Winter -- 90 days usually
												write_iter_count = 0 ! Ignore the incomplete season.
											ELSE
												write_iter_count = 1 ! Full season available
											END IF
											
										ELSE
												write_iter_count = 1 ! In other cases, write the results!
										END IF
										
									ELSE		
										write_iter_count = 0 !not yet end of anything, so don't write.		
									END IF
									
								END IF
								
								!dtHours.eq.24 ONLY works because we assume that the WRF data are in DAILY files...
								!else it would be necessary to check whether it actually isEndOfDay or not... (somehow)
								
								DO writeIndex = 1,write_iter_count,1
									!-------------------------------------------------------------------------------
									! write data to NetCDF file
									
									IF (record_count.eq.1) THEN
										
										data_in = temp_data_InTime(:,:,writeIndex)
										counter = counter_array(writeIndex)
									
									ELSE IF (dtHours == 24) THEN !daily average
										
										data_in = temp_data
									
									ELSE IF (dtHours == 720) THEN !monthly average
										
										data_in = temp_data_Month_Average
									
									ELSE IF (dtHours == 2880) THEN !seasonal average (JESUS...)
										
										data_in = temp_data_Season_Average
										
									END IF
									!!!!!!! ARIS CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!
									
									print *,'write data to NetCDF file'
									print *,'fn_out',fn_out
									sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), NF90_WRITE, ncid )
									IF (sts/=0) EXIT
									print *, 'NF90_OPEN',  sts
									sts = NF90_INQ_VARID(ncid, TRIM(var_cmip(ivar)), x_varid)
									print *, 'NF90_INQ_VARID', ncid
									print *, 'var_cmip(ivar)', var_cmip(ivar)
									print *, 'x_varid', x_varid
									print *, 'counter', counter
									print *, 'xfocus', xfocus
									print *, 'yfocus', yfocus
									IF ( height(ivar) /= -999 ) THEN
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1, 1 /) )
									
									ELSE
										
										sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
										START=(/ 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1 /) )
										
									END IF
									
									print *,'NF90_PUT_VAR', sts
									print *, 'ncid', ncid
									print *, 'x_varid', x_varid
									print *, 'data_in(50:52,50:52)', data_in(50:52,50:52)
									sts = NF90_CLOSE(ncid)
									
									print *, pn_out//"/"//fn_out
									print *, TRIM(var_cmip(ivar)), xfocus, yfocus, counter, ncid, x_varid
									
									!-------------------------------------------------------------------------------
									
								END DO !writeIndex ...
								
								!-------------------------------------------------------------------------------
								! next WRF file contains different number of output intervals
								
								DEALLOCATE(InVarDataRec)
								
								DEALLOCATE(InDateTimeYear)
								DEALLOCATE(InDateTimeMonth)
								DEALLOCATE(InDateTimeDay)
								DEALLOCATE(InDateTimeHour)
								
								!-------------------------------------------------------------------------------
								
								PRINT *, "------------------------------------------------------------"
								CALL CPU_TIME(cpuTe)
								PRINT '("CPU timing for 1 WRF file (e.g. 1 month worth of data) = ",F6.3," sec")',cpuTe-cpuTs
								
							END DO !ifl - specific WRF input file, filelist loop

							!Hmmm.. now which is the order of the following two lines?? THIS OR THE OTHER WAY AROUND???
						END IF
						InDateTimeYearPrev = 0
						
					END DO !nvars - variable loop
					
		END DO ! nvarnml - namelist loop
		
		!END DO !ifrq - different temporal aggregations
		
		
		
		!===============================================================================
		
END PROGRAM ppWRFCMIP

!===============================================================================

SUBROUTINE generateFilelist
	
	USE flhandling
	
	IMPLICIT NONE
	
	CHARACTER (len = 200) :: ifl
	INTEGER :: i, IOstatus, nfl, AllocateStatus
	
	OPEN(2,FILE=tmpfileFL,STATUS='old')
	
	i = 0
	DO
		READ(2,FMT='(a)',IOSTAT=IOstatus) ifl
		IF (IOstatus/=0) EXIT
		i = i + 1
	END DO
	nfl = i
	PRINT *, "number of matching files contained in filelist = ", nfl
	
	IF (ft == 0) THEN
		
		ALLOCATE(fl_wrfout(nfl), STAT=AllocateStatus)
		
	END IF
	
	IF (ft == 1) THEN
		
		ALLOCATE(fl_wrfxtr(nfl), STAT=AllocateStatus)
		
	END IF
	
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	
	REWIND(2)
	DO i = 1,nfl
		IF (ft == 0) THEN
			
			READ(2,FMT='(a)') fl_wrfout(i)
			
		END IF
		
		IF (ft == 1) THEN
			
			READ(2,FMT='(a)') fl_wrfxtr(i)
			
		END IF
	END DO
	
	CLOSE(2)
	
END SUBROUTINE generateFilelist

!===============================================================================

SUBROUTINE CreateRefTimeArray(dt, dataIncludesLeapYearDays)
	
	USE RefTimeVecs
	USE NameListHandling
	
	IMPLICIT NONE
	
	INTEGER :: CheckForLeapyear
	
	CHARACTER (LEN = 3), INTENT(IN) :: dt
	
	INTEGER :: i, j, k, l, counter
	REAL :: dtDecDay
	INTEGER :: tstotYYYY, tstotMM, tstotDD, tstotHH
	INTEGER :: tetotYYYY, tetotMM, tetotDD, tetotHH
	!INTEGER :: temp_tstotYYYY ! Temporary variable to start from December of the
	!previous year than tstotYYYY. NOT HARDCODED TO 1949 ANYMORE)
	INTEGER, DIMENSION(12) :: ndpm
	INTEGER :: ndOverall = 31 ! these are the 31 days of Dec 1949 (Now works for any
	!December )
	!INTEGER :: ndOverall = 0
	INTEGER :: ntspd ! number of timesteps per time-interval
	
	REAL :: indexSum !sum the days to create the fractional index starting from 0.0000
	INTEGER :: nmOverall = 1 ! first December!! number of Months overall -- monthly Average
	INTEGER :: nsOverall !same as above... for season average..
	INTEGER :: seasonIndex !to count for every season... just a counter!
	INTEGER :: dataIncludesLeapYearDays !For 
	
	PRINT *, "CreateRefTimeArray"
	PRINT *, dt
	
	
	! "tstot" contains the absolute starting point, 1949-12-01_00:00:00 (Changes to
	! work for 12-01_00:00:00 of any other year)
	READ( tstot, '(I4,1X,I2,1X,I2,1X,I2)' ) tstotYYYY, tstotMM, tstotDD, tstotHH
	READ( tetot, '(I4,1X,I2,1X,I2,1X,I2)' ) tetotYYYY, tetotMM, tetotDD, tetotHH
	
	! one year before  "tstot"
	!temp_tstotYYYY =  tstotYYYY - 1
	
	! get the overall number of days within the considered timespan
	DO i=tstotYYYY+1,tetotYYYY,1
		
		!DO i=temp_tstotYYYY+1,tetotYYYY,1
		ndOverall = ndOverall + CheckForLeapyear(i, dataIncludesLeapYearDays)
		
	END DO
	!PRINT *, "number of days, overall = ", ndOverall
	
	
	
	
	
	IF (dt == 'mon') THEN
		
		! get the overall number of months within the considered timespan
		DO i=tstotYYYY+1,tetotYYYY,1
			
			!DO i=temp_tstotYYYY+1,tetotYYYY,1
			nmOverall = nmOverall + 12
			
		END DO
		
		!!!!!!!!!!! Allocate THAT many months (12 times the years involved..)
		ALLOCATE( TimeRefArray( nmOverall, 5 ) ) ! index, y, m, d, h         y,x
		!PRINT *, "size and shape of the TimRefArray = ", SIZE(TimeRefArray), &
		!  SHAPE(TimeRefArray)
		
		! fill up the decimal days
		indexSum = 0
		DO i=tstotYYYY,tetotYYYY,1
			
			IF ( CheckForLeapyear( i, dataIncludesLeapYearDays ) == 366 ) THEN
				
				ndpm = (/31,29,31,30,31,30,31,31,30,31,30,31/)
			
			ELSE
				
				ndpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
				
			END IF
			
			DO j = 1,12,1
				
				TimeRefArray( j + 12*(i-tstotYYYY) , 1 ) = indexSum
				indexSum = indexSum + ndpm(j)
				! TimeRefArray( i+1, 1 ) = i * dtDecDay
				
			END DO	
		END DO
		
		! handle the Dec 1949, too complicated to have this in the upcoming loop
		! overall start is at 1949-12-01_00:00:00
		TimeRefArray( 1, 2 ) = 1949
		TimeRefArray( 1, 3 ) = 12.
		TimeRefArray( 1, 4 ) = 1.
		TimeRefArray( 1, 5 ) = 0. !(/ (j, j=0, 24-24/ntspd , 24/ntspd) /) !00  03 06 09 12 15 18 21
		
		!	print*, "(/ (j, j=0, 24-24/ntspd , 24/ntspd) /)", (/ (j, j=0, 24-24/ntspd , 24/ntspd) /)
		! add the rest of the Y M D H information
		counter = 2
		DO i=tstotYYYY+1,tetotYYYY,1
			
			!DO i=temp_tstotYYYY+1,tetotYYYY,1			
			DO j=1,12,1
				
				! sort in on monthly basis
				TimeRefArray( counter, 2) = i
				TimeRefArray( counter, 3) = j
				TimeRefArray( counter, 4) = 1
				TimeRefArray( counter, 5) = 0. !(/(l, l=0, 24-24/ntspd , 24/ntspd) /)
				
				counter = counter + 1
				
			END DO
		END DO
	
	ELSE IF (dt == 'sem') THEN
		

		!Add 4 seasons for all intermediate spanned years...
		
		
		DO i=tstotYYYY+1,tetotYYYY,1
			!DO i=temp_tstotYYYY+1,tetotYYYY,1
			nsOverall = nsOverall + 4 !Adds 4 seasons for 1950, which includes December of 1949
		END DO
		
		!However, the last December cannot count as a whole season (December 2100)..
		! Because the end of the interval ends in December, do not count it as a season (hence the -1)
		ALLOCATE( TimeRefArray( nsOverall+1, 5 ) ) ! index, y, m, d, h         y,x
		!used to be nsOverall - 1 above.. no sense however, in taking SO many precautions..
		
		!PRINT *, "size and shape of the TimRefArray = ", SIZE(TimeRefArray), &
		!  SHAPE(TimeRefArray)
		
		! fill up the decimal days
		indexSum = 0 !These are the 31 days of the first December...
		
		TimeRefArray(1,1) = indexSum !Assign the first value correctly as 0.0 (December of 1949 - 1st season)
		! handle the Dec 1949
		! overall start is at 1949-12-01_00:00:00
		TimeRefArray(1,2) = 1949.
		
		TimeRefArray(1,3) = 12.
		
		TimeRefArray(1,4) = 1.
		
		TimeRefArray(1,5) = 0. !(/ (j, j=0, 24-24/ntspd , 24/ntspd) /) !00  03 06 09 12 15 18 21
		
		
		seasonIndex = 2
		
		DO i=tstotYYYY+1,tetotYYYY,1
			IF ( CheckForLeapyear( i, dataIncludesLeapYearDays ) == 366 ) THEN
				
				ndpm = (/31,29,31,30,31,30,31,31,30,31,30,31/)
			
			ELSE
				
				ndpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
				
			END IF
			
			!Assign end of Winter
			indexSum = indexSum + 31 + ndpm(1) + ndpm(2) !December + Jan + Feb
			TimeRefArray(seasonIndex, 1) = indexSum 
			TimeRefArray(seasonIndex, 2 ) = i
			TimeRefArray(seasonIndex, 3 ) = 3.
			TimeRefArray(seasonIndex, 4 ) = 1.
			TimeRefArray(seasonIndex, 5 ) = 0.
			
			seasonIndex = seasonIndex + 1
			
			DO j = 2,4,1 !Rest of seasons..
				indexSum = indexSum + ndpm((j-1)*3) + ndpm((j-1)*3 + 1) + ndpm((j-1)*3 + 2)
				TimeRefArray( seasonIndex , 1 ) = indexSum
				TimeRefArray(seasonIndex, 2 ) = i
				TimeRefArray(seasonIndex, 3 ) = j*3
				TimeRefArray(seasonIndex, 4 ) = 1.
				TimeRefArray(seasonIndex, 5 ) = 0.
				! TimeRefArray( i+1, 1 ) = i * dtDecDay
				seasonIndex = seasonIndex + 1
			END DO
			
		END DO
		
	
	ELSE !dt is not monthly or season average...
		
		SELECT CASE (dt)
			CASE ('3hr')
				
				dtDecDay = 0.125 ! (or 1.0/8.0)
				ntspd = 1.0 / dtDecDay
			
			CASE ('6hr')
				
				dtDecDay = 1.0 / 4.0 		!changed! --- this was mistakenly set as 1/6.. but there are 4 6-hr intervals
				ntspd = 1.0 / dtDecDay
			
			CASE ('1hr')
				
				dtDecDay = 1.0 / 24.0
				ntspd = 1.0 / dtDecDay
			
			CASE ('day')
				
				dtDecDay = 1.0
				ntspd = 1.0 / dtDecDay
				
				!CASE ('mon')
				!dtDecDay = 365.0/12.0
				!ntspd = 1.0 / dtDecDay
				
				!CASE ('sem')
				!dtDecDay = 365.0/91.0
				!ntspd = 1.0 / dtDecDay
			
			CASE DEFAULT
				
				PRINT *, "invalid time interval specified"
				STOP
				
			END SELECT
			
			
			!!!!!!!!!!! CHANGED ALL INTEGERS WITH NINT
			ALLOCATE( TimeRefArray( ndOverall*ntspd, 5 ) ) ! index, y, m, d, h         y,x
			!PRINT *, "size and shape of the TimRefArray = ", SIZE(TimeRefArray), &
			!SHAPE(TimeRefArray)
			
			! fill up the decimal days
			DO i=0,ndOverall*ntspd-1,1
				
				TimeRefArray( i+1, 1 ) = i / ntspd + mod(i,ntspd) * dtDecDay										
				! TimeRefArray( i+1, 1 ) = i * dtDecDay
				
			END DO
			
			! handle the Dec 1949, too complicated to have this in the upcoming loop
			! overall start is at 1949-12-01_00:00:00
			TimeRefArray( 1:31*ntspd, 2 ) = 1949
			TimeRefArray( 1:31*ntspd, 3 ) = 12.
			
			
			DO i=1,31,1 !possible error to occur here !!!!!!!!!!!!!!!!!!!
				
				TimeRefArray( i*ntspd-(ntspd-1):i*ntspd, 4 ) = i
				TimeRefArray( i*ntspd-(ntspd-1):i*ntspd, 5 ) = (/ (j, j=0, 24-24/ntspd , 24/ntspd) /) !00  03 06 09 12 15 18 21
				
			END DO
			
			print*, "(/ (j, j=0, 24-24/ntspd , 24/ntspd) /)", (/ (j, j=0, 24-24/ntspd , 24/ntspd) /)
			! add the rest of the Y M D H information
			counter = 1
			DO i=tstotYYYY+1,tetotYYYY,1
				
				!DO i=temp_tstotYYYY+1,tetotYYYY,1
				IF ( CheckForLeapyear( i, dataIncludesLeapYearDays ) == 366 ) THEN
					
					ndpm = (/31,29,31,30,31,30,31,31,30,31,30,31/)
				
				ELSE
					
					ndpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
					
				END IF
				
				DO j=1,12,1
					
					! sort in on daily basis
					DO k=1,ndpm(j)
						
						!TimeRefArray( 31*24/ntspd + counter*ntspd-(ntspd-1) : 31*24/ntspd + counter*ntspd , 2) = i
						!TimeRefArray( 31*24/ntspd + counter*ntspd-(ntspd-1) : 31*24/ntspd + counter*ntspd , 3) = j
						!TimeRefArray( 31*24/ntspd + counter*ntspd-(ntspd-1) : 31*24/ntspd + counter*ntspd , 4) = k
						!TimeRefArray( 31*24/ntspd + counter*ntspd-(ntspd-1) : 31*24/ntspd + counter*ntspd , 5) = (/(l, l=0, 24-24/ntspd , 24/ntspd) /)
						
						TimeRefArray( 31*ntspd + counter*ntspd-(ntspd-1) : 31*ntspd + counter*ntspd , 2) = i
						TimeRefArray( 31*ntspd + counter*ntspd-(ntspd-1) : 31*ntspd + counter*ntspd, 3) = j
						TimeRefArray( 31*ntspd + counter*ntspd-(ntspd-1) : 31*ntspd + counter*ntspd , 4) = k
						TimeRefArray( 31*ntspd + counter*ntspd-(ntspd-1) : 31*ntspd + counter*ntspd , 5) = (/(l, l=0, 24-24/ntspd , 24/ntspd) /)
						
						counter = counter + 1
						
					END DO
					
				END DO
			END DO
			
	END IF
	
	!PRINT '(5(F9.3))', TRANSPOSE(TimeRefArray(1:300,:))
	!PRINT '(5(F9.3))', TRANSPOSE(TimeRefArray(200000:200040,:))
	
	!-------------------------------------------------------------------------------
	
	!CONTAINS
	
	!  INTEGER FUNCTION CheckForLeapyear( year )
	
	!    IMPLICIT NONE
	
	!    INTEGER, INTENT(IN) :: year
	
	!    IF ( ( MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0 ) .OR. ( MOD(year, 400) == 0 ) ) THEN
	!      CheckForLeapyear = 366
	!    ELSE
	!      CheckForLeapyear = 365
	!    END IF
	
	!  END FUNCTION CheckForLeapyear
	
END SUBROUTINE CreateRefTimeArray


! General function
INTEGER FUNCTION CheckForLeapyear( year, leap )
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: year
	INTEGER :: leap
	
	IF (leap) THEN
	
		IF ( ( MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0 ) .OR. ( MOD(year, 400) == 0 ) ) THEN
		
			CheckForLeapyear = 366
		
		ELSE
			
			CheckForLeapyear = 365
			
		END IF
	
	ELSE
	
		CheckForLeapyear = 365
		
	END IF

END FUNCTION CheckForLeapyear																																							