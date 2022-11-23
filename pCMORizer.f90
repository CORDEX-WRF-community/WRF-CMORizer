!-------------------------------------------------------------------------------
! passing allocatable arrays between main program and external subroutine

MODULE FilelistHandling

  IMPLICIT NONE
  SAVE

  CHARACTER (len = 200):: tmpfileFL_std, tmpfileFL_xtrm, tmpfileFL_3d
  CHARACTER (len = 200), DIMENSION(:), ALLOCATABLE :: &
    fl_wrfout   , &
    fl_wrfxtr   , & 
    fl_wrfpres  , &
    fl_input
  INTEGER :: ft

END MODULE FilelistHandling

!-------------------------------------------------------------------------------
! index, Y, M, D, H, depending on the "frequency" of the dataset, i.e. time
! intervals

MODULE RefTimeVecs

  IMPLICIT NONE
  SAVE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: TimeRefArray
  !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TimeRefArray
  !INTEGER :: Tyear_start, Tyear_end

END MODULE RefTimeVecs

!-------------------------------------------------------------------------------
! namelist handling

MODULE NamelistHandling

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: nvars = 39 ! maximum number of vars per namelist

  CHARACTER (len = 300) :: Conventions, conventionsURL, contact, experiment, &
    driving_experiment, experiment_id, driving_model_id, driving_model_ensemble_member, &
    driving_experiment_name, institution, institute_id, model_id, &
    rcm_version_id, project_id, CORDEX_domain, product, references

  CHARACTER (len = 1000) :: title, institute_run_id, comment, &
    nesting_levels, comment_nesting, comment_1nest, comment_2nest

  CHARACTER (LEN = 100), DIMENSION(nvars) :: var_wrf, var_cmip, standard_name, &
    long_name, units, filetype, cm1hr, cm3hr, cm6hr, cmDay, cmMon, cmSea, positive
  INTEGER, DIMENSION(nvars):: height, cordexID
  LOGICAL, DIMENSION(nvars):: time1hr, time3hr, time6hr, timeDay, timeMon, timeSea, &
     interpolate

  CHARACTER (len = 300) :: DirInputSimResRoot, DirOutputPostProRoot, domain

  INTEGER ::  nx, ny, nz, xoffset, yoffset, xfocus, yfocus
  CHARACTER (len = 4) :: ts, te
  CHARACTER (len = 19) :: tstot, tetot

  CHARACTER (len = 300) :: PnFnGeo
  
  CHARACTER (len = 3) :: model

  LOGICAL :: aggregation_yearly, aggregation_monthly, aggregation_individually
  CHARACTER (len = 19) :: tsact, teact
  INTEGER :: nvar

  NAMELIST / globalvars / Conventions, contact, experiment_id, experiment, &
    driving_experiment, driving_model_id, driving_model_ensemble_member, &
    driving_experiment_name, institution, institute_id, model_id, &
    rcm_version_id, project_id, CORDEX_domain, product, references, &
    conventionsURL

  NAMELIST / globalvars_additional / comment, institute_run_id, title, &
    nesting_levels, comment_nesting, comment_1nest, comment_2nest

  NAMELIST / vars / var_wrf, var_cmip, standard_name, long_name, units, &
    height, time1hr, time3hr, time6hr, timeDay, timeMon, timeSea, filetype, &
    cm1hr, cm3hr, cm6hr, cmDay, cmMon, cmSea, interpolate, cordexID, positive

  NAMELIST / filesystem / DirInputSimResRoot, DirOutputPostProRoot, domain

  NAMELIST / model_config / ts, te, nx, ny, nz, xoffset, yoffset, xfocus, &
    yfocus, tstot, tetot

  NAMELIST / static_fields / PnFnGeo

  NAMELIST / tool_config / nvar, aggregation_yearly, aggregation_monthly, &
    aggregation_individually, tsact, teact, model

END MODULE NamelistHandling

!===============================================================================

PROGRAM WRFCMORizer

USE FilelistHandling
USE RefTimeVecs
USE NamelistHandling

USE netcdf

USE MPI
!USE OMP_LIB

IMPLICIT NONE

!===============================================================================

INTERFACE

  SUBROUTINE GenerateFilelist
  END SUBROUTINE GenerateFilelist

  SUBROUTINE CreateRefTimeArray( dt )
    IMPLICIT NONE
    CHARACTER (LEN = 3), INTENT(IN) :: dt
  END SUBROUTINE CreateRefTimeArray

  ! HTr: modified calculation of mean sea level pressure adopted from 
  ! wrf_interp.F90	
  SUBROUTINE calcslp(slp,pres,qv,tk1,ght,nz,ns,ew,T00)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nz,ns,ew
    REAL, DIMENSION(:,:), INTENT(INOUT) :: slp
    REAL, DIMENSION(:,:,:), INTENT(IN) :: pres,qv,tk1,ght
    REAL, INTENT(IN) :: T00
  END SUBROUTINE calcslp

  SUBROUTINE calcslptwo(slp, PP, P_s, PHI_s, T_L, nz, ns, ew)
    IMPLICIT NONE
    REAL, DIMENSION(:,:),   INTENT(OUT) :: slp
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: PP
    REAL, DIMENSION(:,:),   INTENT(IN)  :: P_s
    REAL, DIMENSION(:,:),   INTENT(IN)  :: PHI_s
    REAL, DIMENSION(:,:),   INTENT(IN)  :: T_L
    INTEGER,                INTENT(IN)  :: nz, ns, ew
  END SUBROUTINE calcslptwo

END INTERFACE

!===============================================================================
! filenames
! do not use the individual variable namelists here, they are set later on in a 
! loop; create symbolic links to have different base namelists, depending on the
! dataset and experiment

CHARACTER (len = *), PARAMETER :: fnNMLexp = "runctrl.current.nml"
CHARACTER (len = 100), DIMENSION(:), ALLOCATABLE :: fnNMLvar
CHARACTER (len = 200) :: pn_out, fn_out, iflWRFin

!-------------------------------------------------------------------------------

REAL, PARAMETER :: cp = 1004.0 ! [J kg-1 K-1]
REAL, PARAMETER :: R = 287.04 ! [J kg-1 K-1]
REAL, PARAMETER :: L = 2501000.0 ! [J kg-1]
REAL, PARAMETER :: a = 610.78 ! [Pa]
REAL, PARAMETER :: b = 17.27 ! 17.2693882
REAL, PARAMETER :: c = 273.15 !
REAL, PARAMETER :: d = 35.86 !
REAL, PARAMETER :: n = L*0.622*a/cp !
REAL, PARAMETER :: mv = 1.e20 ! missing value as specified
REAL, PARAMETER :: gr = 9.81
REAL, PARAMETER :: epsil = 0.6220

! auxilliary vars, just needed during development
! INTEGER, PARAMETER :: nt = 8

! constants for hurs calc
real(kind=8), parameter :: a1=-2.8365744e3,a2=-6.028076559e3,a3=19.54263612,a4=-2.737830188e-2
real(kind=8), parameter :: a5=1.6261698e-5,a6=7.0229056e-10,a7=-1.8680009e-13,a8=2.7150305
real(kind=8), parameter :: b1=-5.8666426e3,b2=22.32870244,b3=1.39387003e-2,b4=-3.4262402e-5
real(kind=8), parameter :: b5=2.7040955e-8,b6=6.7063522e-1
real(kind=8), parameter :: fi1=3.62183e-4,fi2=2.6061244e-5,fi3=3.8667770e-7,fi4=3.8268958e-9, &
  fi5=-10.7604,fi6=6.3987441e-2,fi7=-2.6351566e-4,fi8=1.6725084e-6
real(kind=8), parameter :: fw1=3.536240e-4,fw2=2.932836e-5,fw3=2.616898e-7,fw4=8.581361e-9, &
  fw5=-10.75880,fw6=6.326813e-2,fw7=-2.536893e-4,fw8=6.340529e-7

! new netCDF file
INTEGER :: ncid, ncidin, ncidin0
INTEGER :: lon_dimid, lat_dimid, rec_dimid, height_dimid, &
  nb2_dimid, x_dimid, y_dimid, plev_dimid, depth_dimid
INTEGER :: varid, x_varid, lon_varid, lat_varid, rlon_varid, rlat_varid, &
  rotated_pole_varid, height_varid, rec_varid, pp_varid, pb_varid, ph_varid, &
  phb_varid, qv_varid, qc_varid, qi_varid, qr_varid, qs_varid, &
  theta_varid, t2_varid, recbnds_varid, rainnc_varid, &
  rainc_varid, snownc_varid, u10_varid, v10_varid, u_varid, v_varid, w_varid, &
  sfcevp_varid, potevp_varid, sfroff_varid, udroff_varid, acsnom_varid, q2_varid, &
  sinalpha_varid, cosalpha_varid, plev_varid, plevbnds_varid, psfc_varid, &
  landmask_varid, xland_varid, swdown_varid, depth_varid, soillayerbnds_varid

! input data general query
INTEGER :: ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in !!!, formatp_in

! inputs, number of elements
INTEGER :: nvar_nml
! record variable in input data
INTEGER :: InVarIdRec, InDimLenRec !!!, InVarNdimsRec
CHARACTER (len = NF90_MAX_NAME) :: InDimNameRec !!!, InVarNameRec
!!!INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: InDimIdsRec
CHARACTER (len = 50) :: fl_filter
CHARACTER (len = 19), DIMENSION(:), ALLOCATABLE :: InVarDataRec

! data
REAL, DIMENSION(:,:), ALLOCATABLE :: &
  data_in       , &
  psl_in        , &
  t2_in         , &
  q2_in         , &
  cldfra_inv    , &
  u10_in        , &
  v10_in        , &
  cape          , &
  cin           , &
  lcl           , &
  lfc           , &
  prw           , &
  clwvi         , &
  clivi         , &
  sinalpha_in   , &
  cosalpha_in   , &
  !!var_pl      , &
  psfc_in       , &
  tmp_2d        , &
  tmp1_2d       , &
  tmp2_2d       , &
  tmp3_2d       , &
  tmp4_2d       , &
  rainc_max_in  , &
  rainnc_max_in , &
  landmask_in   , &
  xland_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
  data_in_3D  , &
  pp_in       , &
  pb_in       , &
  ph_in       , &
  phb_in      , &
  qv_in       , &
  qvs         , &
  qc_in       , &
  qr_in       , &
  qi_in       , &
  qs_in       , &
  theta_in    , &
  t_in        , &
  ph_fl       , &
  p_in        , &
  u_in        , &
  v_in        , &
  w_in        , &
  var3d_in    , &
  !var_pl      , &
  potevp_in   , &
  rainnc_in   , &
  rainc_in    , &
  rad_in      , &
  hflux_in    , &
  t_p         , &
  snownc_in   , &
  acsnom_in   , &
  GeoInLonLat , &
  sfcevp_in   , &
  sfroff_in   , &
  udroff_in   , &
  swdown_in   , &
  tmp_3d
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: &
  cldfra_in   , &
  smois_in
REAL, DIMENSION(:), ALLOCATABLE :: &
  GeoInRLat   , &
  GeoInRLon

! target pressure levels [Pa]
REAL, DIMENSION(6) :: pout = (/100000.,92500.,85000.,70000.,50000.,20000./)

! time vec stuff
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: &
  TimeRefArraySubset, &
  Time_bnds
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: &
  TimeRefArraySubsetMean

! bucket system
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
  i_rainnc_in , &
  i_rainc_in  , &
  i_rad_in
REAL :: bucket_mm, bucket_J

! (het): base state temperature is made flexible in newer versions of WRF. The
! actual value is stored in variable T00
! also the base state pressure (P00) is kept flexible now.
REAL, DIMENSION(1) :: T00, P00
INTEGER :: t00_varid, p00_varid

! variable for adopted vertical interpolation
REAL :: zg_pout

! soil layer thickness may vary from simulation to simulation
! either DZS is read from WRF outout, or it can be hardcoded, see below
REAL, DIMENSION(:), ALLOCATABLE :: DZS
REAL, DIMENSION(:), ALLOCATABLE :: DZSmp

! subsurface layer thickness prescribed, e.g., compatible to ParFlow, hardcoded (hc)
! any number of soil layers is possible
!REAL, DIMENSION(4) :: DZShc = (/0.10D0, 0.30D0, 0.60D0, 1.00D0/)
REAL, DIMENSION(4) :: DZShc = (/0.1, 0.3, 0.6, 1.0/)
REAL, DIMENSION(:,:), ALLOCATABLE :: soillayer_bnds

! meta information about the geographic projection used (coordinates of the rotated pole)
REAL :: GeoNPLat, GeoNPLon

REAL :: t_ii, dtHours

! time and date handling
CHARACTER (len = 3), DIMENSION(:), ALLOCATABLE :: frequency
!!!INTEGER :: WRFfileNyears, WRFfileNmonths
INTEGER, DIMENSION(:), ALLOCATABLE :: InDateTimeYear, InDateTimeMonth, &
  InDateTimeDay, InDateTimeHour, InDateTimeMinute, InDateTimeSecond      !, WRFfileIyears, WRFfileImonths
REAL, DIMENSION(:), ALLOCATABLE :: InDateTimeCombined
CHARACTER (LEN=4) :: InDateTimeYearStr
CHARACTER (LEN=2) :: &
  InDateTimeMonthStr    , &
  FirstHourStr          , &
  FirstMinuteStr        , &
  LastDayStr            , &
  LastHourStr           , &
  LastMinuteStr
INTEGER :: InDateTimeYearPrev = 0, InDateTimeMonthPrev = 0
CHARACTER (len = 12) :: FileNameStartDateTime, FileNameEndDateTime
INTEGER :: tsactYear, tsactMonth, tsactDay, tsactHour, tsactMinute, tsactSecond, & 
  teactYear, teactMonth, teactDay, teactHour, teactMinute, teactSecond
CHARACTER (LEN=4) :: tsactYearStr, teactYearStr
CHARACTER (LEN=2) :: tsactMonthStr, tsactDayStr, tsactHourStr, tsactMinuteStr, &
  tsactSecondStr, teactMonthStr, teactDayStr, teactHourStr, teactMinuteStr, &
  teactSecondStr
REAL(KIND=8) :: tsact_singlenumber, teact_singlenumber

REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp2D_singlenumber
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp2D

CHARACTER (LEN = 100), DIMENSION(nvars) :: cell_methods
LOGICAL, DIMENSION(nvars) :: procflag

LOGICAL :: fractSeaIce = .FALSE.

!-------------------------------------------------------------------------------
! statistics

REAL :: stat_mean, slope

!-------------------------------------------------------------------------------
! general
! some special switch for input handling

INTEGER :: i, j, k, sts, ivar, ifrq, ifl, it, counter, np, nl, ii, ivarnml, &
  prevpass = 0, AllocateStatus, DeAllocateStatus
LOGICAL :: FileExists, newpass, time_match, calc = .TRUE., inputtimesteptruncate = .FALSE.
REAL :: cpuTs, cpuTe
INTEGER, ALLOCATABLE :: ipos(:)
! choose how psl is calculated, 0 OK, 1 OK, 2 OK, no not change throughout
! see preamble for detailed information
INTEGER :: calc_slp_type = 2 
! (HTr): choose how geopotantial height and hus is vertically interpolated to given pressure levels
! see preamble for detailed information
INTEGER :: vint_type = 1

!INTEGER :: deflate_level, endianness
!logical :: contiguous, shuffle, fletcher32
!integer :: xtype, ndims, natts

CHARACTER (len = 3) :: FileNrStr

!-------------------------------------------------------------------------------
! system calls

CHARACTER (len = *), PARAMETER :: cmdUUID = "uuidgen -t > tmpfileUUID"
CHARACTER (len = 37) :: trackingID

CHARACTER (len = *), PARAMETER :: cmdDate = "date -u +%Y-%m-%d-T%H:%M:%SZ > tmpfileDate"
CHARACTER (len = 21) :: creationDate

!===============================================================================

! MPI

INTEGER, DIMENSION(:), ALLOCATABLE :: ivar_list
INTEGER :: ierr, rank, numtasks

! Initialise MPI evironment
  CALL MPI_INIT(ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_INIT"

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_COMM_RANK"

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_COMM_WORLD"

PRINT *, "============================================================"
PRINT *, "WRF RCM CMORizer"
PRINT *, "============================================================"

PRINT *, "============================================================"
PRINT *, "*** CENTRAL NAMELIST READING ***"
PRINT *, fnNMLexp

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

OPEN(2,FILE=fnNMLexp)
READ(UNIT=2,NML=tool_config)
CLOSE(2)

!-------------------------------------------------------------------------------
! allocate main data input array outside the loops based on nml entries
! dummy allocation of the ref time array -> has to maintain its values during
! several looping constructs and is de-allocated before the initial allocation

ALLOCATE( data_in( xfocus, yfocus ), STAT=sts )
IF (sts /= 0) STOP "*** Not enough memory on this device, stopping***"

! this is a dirty hack for the mrsol and tsl variables
! at this point it is not even clear whether this variable is needed at all
ALLOCATE( data_in_3D( xfocus, yfocus, 4 ), STAT=sts )
IF (sts /= 0) STOP "*** Not enough memory on this device to create data_in_3D, stopping***"

ALLOCATE( TimeRefArraySubset(2,2) )
ALLOCATE( TimeRefArraySubsetMean(2) )
ALLOCATE( Time_bnds(2,2) )

!-------------------------------------------------------------------------------
! as indices are relative to 1 and the model_config namelist sets the boundary
! relaxation zone with its total width, this adjustment is needed to get the 
! correct start index in x and y direction, i.e. usually xoffset=10 then the 
! 1st index to read is 11 see also https://www.unidata.ucar.edu/software/netcdf/
! docs-fortran/f90-variables.html#f90-reading-data-values-nf90_get_var
PRINT *, "xoffset, yoffset before: ", xoffset, yoffset
xoffset = xoffset + 1
yoffset = yoffset + 1
PRINT *, "xoffset, yoffset after: ", xoffset, yoffset

!-------------------------------------------------------------------------------
! get the invariant vars which have to be added all the time from seperate file
! lon, lat, rlon, rlat, mass grid, etc.
! normal order: x y z t

PRINT *, "============================================================"
PRINT *, "*** STATIC FIELDS ***"
PRINT *, TRIM(PnFnGeo)

ALLOCATE( GeoInLonLat(xfocus, yfocus, 2) )
PRINT *, SHAPE(GeoInLonLat)
ALLOCATE( landmask_in(xfocus, yfocus) )
PRINT *, SHAPE(landmask_in)
ALLOCATE( GeoInRLat(yfocus) )
ALLOCATE( GeoInRLon(xfocus) )

sts = NF90_OPEN(TRIM(PnFnGeo), IOR(NF90_NOWRITE, NF90_MPIIO), ncidin, &
      comm = MPI_COMM_WORLD, info = MPI_INFO_NULL )

sts = NF90_INQ_VARID(ncidin, "lon", varid)
sts = NF90_GET_VAR(ncidin, varid, GeoInLonLat(:, :, 1), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /))

sts = NF90_INQ_VARID(ncidin, "lat", varid)
sts = NF90_GET_VAR(ncidin, varid, GeoInLonLat(:, :, 2), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /))

sts = NF90_INQ_VARID(ncidin, "rlon", varid)
sts = NF90_GET_VAR(ncidin, varid, GeoInRLon(:), &
  START = (/ xoffset, 1, 1 /), COUNT = (/ xfocus, 1, 1 /))
!PRINT *, "GeoInRLon shape (EUR-15, 339): ", SHAPE(GeoInRLon)
!PRINT *, "GeoInRLon: ", GeoInRLon

sts = NF90_INQ_VARID(ncidin, "rlat", varid)
sts = NF90_GET_VAR(ncidin, varid, GeoInRLat(:), &
  START = (/ 1, yoffset, 1 /), COUNT = (/ 1, yfocus, 1 /))
!PRINT *, "GeoInRLat shape (EUR-15, 330): ", SHAPE(GeoInRLat)
!PRINT *, "GeoInRLat: ", GeoInRLat

!sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "POLE_LAT", GeoNPLat)
!sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "POLE_LON", GeoNPLon)
GeoNPLat=39.25
GeoNPLon=-162.0

! binary, land=1, water=0, compared to LANDUSEF class 16, water fraction, 
! >0.5 water fraction is water, including large lakes
sts = NF90_INQ_VARID(ncidin, "FR_LAND", varid)
sts = NF90_GET_VAR(ncidin, varid, landmask_in(:, :), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /))
! TODO: make this a binary mask
  
sts = NF90_CLOSE(ncidin)

!PRINT *, "rlon = "
!PRINT *, SHAPE(GeoInRLon)
!PRINT *, GeoInRLon

!PRINT *, "rlat = "
!PRINT *, SHAPE(GeoInRLat)
!PRINT *, GeoInRLat

!PRINT *, "GeoInLonLat shape = ", SHAPE(GeoInLonLat)
!PRINT *, "lon = ", GeoInLonLat(:,:,1)
!PRINT *, "lat = ", GeoInLonLat(:,:,2)

!PRINT *, 'landmask_in shape = ', SHAPE(landmask_in)

!-------------------------------------------------------------------------------
! setup of loop control vectors

ALLOCATE ( frequency(7) )
frequency(1) = "1hr"
frequency(2) = "3hr"
frequency(3) = "6hr"
frequency(4) = "day"
frequency(5) = "mon"
frequency(6) = "sem"
frequency(7) = "fx"

ALLOCATE ( fnNMLvar(17) )
fnNMLvar(1) = "runctrl.vars.nml_pr_tas_1hr_test" ! OK
fnNMLvar(2) = "runctrl.vars.nml_vars_on_plevels" ! OK
fnNMLvar(3) = "runctrl.vars.nml" ! OK
fnNMLvar(4) = "runctrl.vars.nml_pr_mrso" ! OK
fnNMLvar(5) = "runctrl.vars.nml_minmax" ! OK
fnNMLvar(6) = "runctrl.vars.nml_evp_roff" ! ok not yet tested: evspsbl, evspsblpot
fnNMLvar(7) = "runctrl.vars.nml_water_column" ! OK
fnNMLvar(8) = "runctrl.vars.nml_radiation" ! OK
fnNMLvar(9) = "runctrl.vars.nml_snow" ! ok not yet tested in winter: sic
fnNMLvar(10) = "runctrl.vars.nml_cape" ! OK
fnNMLvar(11) = "runctrl.vars.nml_weathertyping"  ! new from HTr, not implemented
fnNMLvar(12) = "runctrl.vars.nml_psl"            ! new from HTr, not implemented
fnNMLvar(13) = "runctrl.vars.nml_Alvaro"         ! new from HTr, not implemented
fnNMLvar(14) = "runctrl.vars.nml_alb"         ! new from HTr, not implemented
fnNMLvar(15) = "runctrl.vars.nml_KlausStefan"         ! new from HTr, not implemented
fnNMLvar(16) = "runctrl.vars.nml_Alvaro_200"         ! new from HTr, not implemented
fnNMLvar(17) = "runctrl.vars.nml_Sophie"         ! new from HTr, not implemented


!-------------------------------------------------------------------------------
! individual vars contain information on whether they are treated or not, i.e.
! the tool may loop over all frequencies and the namelist controls what is to 
! be done
! this is the outer loop as it also controls the averaging; the main processing
! is done during the 1st pass, i.e., the tier-2 processing at dt=1h or 3h
! ATTENTION: averaging functionality not it reintegrated
!            during testing too lazy to set all temporally aggregated vars to 
!            FALSE
!            -> just run the outer loop over 1 frequency only

!DO ifrq = 1, SIZE(frequency), 1 
DO ifrq = 1, 1, 1 ! 1hr
!DO ifrq = 4, 4, 1 ! check min/max day

  PRINT *, "============================================================"
  PRINT *, "freq = ", frequency(ifrq)

!-------------------------------------------------------------------------------
! - get a file list of all wrfout, wrfxtrm and wrfpress files -- if they exist
! - use regex to refine the ls output and the filelist
! - non-std F95, works for gfortran (fct & subroutine) + ifort
! - in order to not check the filesystem too often, this is done at this point
! - also the search results might be limited using the filter
! - this can be made more efficient, but then it is not so generic anymore
!   depending on the depth of the search path find takes some time
! - using find in combination with a data root directory should be OK with most
!   ways users store their files
  
  PRINT *, "============================================================"
  PRINT *, "*** FILELIST CREATION ***"
  
  tmpfileFL_std = "tmpfileFLstd" // TRIM(domain)
  tmpfileFL_3d = "tmpfileFL3d" // TRIM(domain)
  tmpfileFL_xtrm = "tmpfileFLxtrm" // TRIM(domain)

  ! creates a year range, can be expanded also to use months
  IF ( (ts == "0000") .AND. (te == "0000") ) THEN
    fl_filter = ""
  ELSE
    !fl_filter = "{" // ts // ".." // te // "}" ! does not work with SYSTEM call
    fl_filter = ts
  END IF
  PRINT *, "fl_filter from runctrl nml: ", fl_filter

  ! decisive for the filelist which in turn determines whether, e.g., the last
  ! timestep in the output file can be generated when averaging
  ! there is an interplay between the scripting and this pattern matching, which
  ! has to be adjusted individually in the end
  !CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/" // TRIM(domain) // "/" // TRIM(ts) // " -name wrfout*" // TRIM(domain) // "*_" // "*.nc | sort > " // tmpfileFL)
  !CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/" // TRIM(domain) // "/" // TRIM(ts) // " -name wrfout*" // TRIM(domain) // "*_" // TRIM(fl_filter) // "*.nc | sort > " // tmpfileFL)
  !CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/" // TRIM(domain) // " -name wrfout*" // TRIM(domain) // "*_" // TRIM(ts) // "*.nc -o -name wrfout*" // TRIM(domain) // "*_" // TRIM(te) // "*.nc | sort > " // tmpfileFL)

  IF ( rank == 0 ) THEN
    !CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'lfdd*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*.nc' | sort > " // tmpfileFL_std)
    CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'lffd*" // ".nc' | sort > " // tmpfileFL_std)
  END IF
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
  ft = 0 ! file type
  CALL GenerateFilelist

  IF ( rank == 0 ) THEN 
    CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfxtrm*" // TRIM(domain) // "*" // &
      TRIM(fl_filter) // "*.nc' | sort > " // tmpfileFL_xtrm)
  END IF
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
  ft = 1
  CALL GenerateFilelist

  IF ( rank == 0) THEN
    CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfpress*" // TRIM(domain) // "*" // &
      TRIM(fl_filter) // "*.nc' | sort > " // tmpfileFL_3d)
  END IF
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
  ft = 2
  CALL GenerateFilelist
 
  DO i=1,SIZE(fl_wrfout(:)),1
    PRINT '(100A)', fl_wrfout(i)
  END DO
  DO i=1,SIZE(fl_wrfxtr(:)),1
    PRINT '(100A)', fl_wrfxtr(i)
  END DO
  DO i=1,SIZE(fl_wrfpres(:)),1
    PRINT '(100A)', fl_wrfpres(i)
  END DO
  
!-------------------------------------------------------------------------------
! creation of the main reference array
! independent of the actual timespan under processing
! usually it starts earlier or at the same date/time and ends at the same 
! date/time or later
! input time information in WRF is 2009-06-20_00:00:00 for extremes

  PRINT *, "============================================================"
  PRINT *, "*** TIME REFERENCE ARRAY ***"
  
  CALL CreateRefTimeArray( frequency(ifrq) )
  
  PRINT *, "size & shape of the TimRefArray = ", SIZE(TimeRefArray), &
    SHAPE(TimeRefArray)
  PRINT *, "SIZE(TimeRefArray,1): ", SIZE(TimeRefArray,1) 
  !PRINT *, "SHAPE(TimeRefArray,1): ", SHAPE(TimeRefArray,1)

!-------------------------------------------------------------------------------
! loop over the different namelists, each containing a specific set of related
! variables, this offers more flexibility in using the tool than putting all
! vars into a single namelist; choose just specific namelists from list above if
! you want to postprocess just specific variables or create your own variable 
! combinations
  
!  DO ivarnml = 1, 9, 1 ! loop over all regular namelists
!  DO ivarnml = 2, 2, 1 ! recommended to all for first steps and testing: nml #1
!  DO ivarnml = 17, 17, 1 ! tas and precip, only. For Sophie
   DO ivarnml = 3, 3, 1 ! tas and precip, only
!  DO ivarnml = 13, 13, 1 ! Alvaro paper
!  DO ivarnml = 15, 15, 1 ! Klaus Stefan
!  DO ivarnml = 16, 16, 1 ! Alvaro 200 hPa
!  DO ivarnml = 14, 14, 1 ! albedo investigations
!  DO ivarnml = 1, 2, 1 ! ICTP paper data contrib
!  DO ivarnml = 5, 5, 1 ! test min/max
  
    PRINT *, "============================================================"
    PRINT *, "var. namelist nr. and name: ", ivarnml, TRIM(fnNMLvar(ivarnml))

    ! read the very specific namelist
    OPEN(2,FILE=TRIM(fnNMLvar(ivarnml)))
      READ(UNIT=2,NML=vars)
    CLOSE(2)

    ! must have read the namelist already
    ! get the call methods from the namelist, determines what has to be done with 
    ! a variable 
    SELECT CASE (frequency(ifrq))
    CASE ('1hr')
      cell_methods(:) = cm1hr(:)
      procflag(:) = time1hr(:)
      dtHours = 1.
    CASE ('3hr')
      cell_methods(:) = cm3hr(:)
      procflag(:) = time3hr(:)
      dtHours = 3.
    CASE ('6hr')
      cell_methods(:) = cm6hr(:)
      procflag(:) = time6hr(:)
      dtHours = 6.
    CASE ('day')
      cell_methods(:) = cmDay(:)
      procflag(:) = timeDay(:)
      dtHours = 24.
    CASE ('mon')
      STOP "monthly aggregation not yet implemented"
      cell_methods(:) = cmMon(:)
      procflag(:) = timeMon(:)
    CASE ('sem')
      STOP "seasonal aggregation not yet implemented"
      cell_methods(:) = cmSea(:)
      procflag(:) = timeSea(:)
    CASE DEFAULT
      PRINT *, "invalid time interval specified"
      STOP
    END SELECT

    ! nvar_nml might be determined through the namelist itself
    ! sequence here does not matter
    SELECT CASE (TRIM(fnNMLvar(ivarnml)))
    CASE ("runctrl.vars.nml") ! OK
      nvar_nml = nvar !36 38  ! 2 or 13
    CASE ("runctrl.vars.nml_evp_roff") ! ok -- see above
      nvar_nml = 4
    CASE ("runctrl.vars.nml_water_column") ! OK
      nvar_nml = 3
    CASE ("runctrl.vars.nml_vars_on_plevels") ! OK
      nvar_nml = 36
    CASE ("runctrl.vars.nml_pr_mrso") ! OK --> new: tsl mrlsl, from 4 to 6
      nvar_nml = 6
    CASE ("runctrl.vars.nml_snow") ! ok -- see above
      nvar_nml = 6 
    CASE ("runctrl.vars.nml_radiation") ! OK
      nvar_nml = 10
    CASE ("runctrl.vars.nml_cape") ! OK
      nvar_nml = 2
    CASE ("runctrl.vars.nml_pr_tas_1hr_test") ! OK
      nvar_nml = 12 
    CASE ("runctrl.vars.nml_weathertyping") ! new stuff from HTr, no nml yet available
      nvar_nml = 6
    CASE ("runctrl.vars.nml_psl") ! new stuff from HTr, no nml yet available
      nvar_nml = 1
    CASE ("runctrl.vars.nml_minmax") ! OK
      nvar_nml = 4
    CASE ("runctrl.vars.nml_Alvaro") 
      nvar_nml = 27
    CASE ("runctrl.vars.nml_alb") 
      nvar_nml = 8
    CASE ("runctrl.vars.nml_KlausStefan") 
      nvar_nml = 11
    CASE ("runctrl.vars.nml_Alvaro_200") 
      nvar_nml = 6
    CASE ("runctrl.vars.nml_Sophie") 
      nvar_nml = 2
    END SELECT
  
    PRINT *, "number of vars inside current namelist: nvar_nml = ", nvar_nml
  
!-------------------------------------------------------------------------------
! loop over all vars in the individual namelist
! OR for testing choose just specific variables from namelist 
! (look up var column position in individual namelist)
! better avoid this kind of filtering: create a new namelist or switch the vars
! properly on/off for the respective temporal aggregation level

  ALLOCATE( ivar_list(nvar_nml) )
  DO ivar = 1, nvar_nml, 1
    ivar_list(ivar) = ivar
  END DO


  CALL MPI_SCATTER( ivar_list, 1, MPI_INT, ivar, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_SCATTER"

!    DO ivar = 1, nvar_nml, 1
!    DO ivar = 8,9, 1
  
      PRINT *,"============================================================"
      PRINT *, "*** ", TRIM(var_cmip(ivar)), procflag(ivar), " ***"
 
      ! shall this variable be processed at all?
      ! either control by i) including/excluding brute-force from a namelist, by 
      ! customizing the namelist list or ii) by setting "timeXX" per temporal 
      ! resolution/aggregation to TRUE or FALSE (=recommended way)
      IF (procflag(ivar)) THEN

      !this is the place to put the ref time vec calculation to incorporate
      !the specifics of the respective variable, like point or mean data in
      !combinaiton with the ifrq  
 
!-------------------------------------------------------------------------------
! loop over the filelist
! content of filelist is defined by filename patterns in system call
! three types of filelists: fl_wrfout, fl_wrfxtr (all depends on the variable
! which filelist shall be used)
! the filelist loops does not care whether (a) the filelist spans multiple
! simulated years or just a single day or hour, neither does it care (b) how
! many files are contained and (c) whether a single file is oerlapping a 
! beginnning or end of the desired temporal unit (year, month, individual, etc.)

      IF ( filetype(ivar) == "s" ) THEN ! variable is in "s"tandard wrfout file
        IF (ALLOCATED(fl_input)) DEALLOCATE( fl_input )
        ALLOCATE( fl_input ( SIZE(fl_wrfout) ), STAT=sts )
        fl_input = fl_wrfout
      ELSE IF ( filetype(ivar) == "x" ) THEN ! variable is in e"x"tremes wrfxtrm file
        IF (ALLOCATED(fl_input)) DEALLOCATE( fl_input )
        ALLOCATE( fl_input ( SIZE(fl_wrfxtr) ), STAT=sts )
        fl_input = fl_wrfxtr
      ELSE IF ( filetype(ivar) == "p" ) THEN ! variable is in "p"ress wrfpress file
        IF (ALLOCATED(fl_input)) DEALLOCATE( fl_input )
        ALLOCATE( fl_input ( SIZE(fl_wrfpres) ), STAT=sts )
        fl_input = fl_wrfpres
      ELSE
        STOP "filetype in the variable namelist must be one of s, x, p"
      END IF

      !DO ifl = 1, 1, 1 ! testing: loop over specific entry in filelist (e.g. just January)
      DO ifl = 1, SIZE(fl_input), 1 ! operational: loop over complete filelist

        PRINT *,"============================================================"
        PRINT *, "filelist filetype = ", filetype(ivar) 
        PRINT *, "# files to process = ", SIZE(fl_input)

        ! measure CPU time for 1 variable and file, including all timesteps
        CALL CPU_TIME(cpuTs)
        PRINT *, "-----------------------------------------------------------"

        iflWRFin = fl_input(ifl)
        PRINT *, "this is the file to work on now:"
        PRINT '(100A)', TRIM(iflWRFin)
  
!-------------------------------------------------------------------------------
! which timespan is covered by the CURRENT WRF input file
! even in the same filelist file, each input file may cover a different timespan
! this determines how many times the tool has to loop over the inputs
! format of 'Times': 2009-06-20_08:00:00
! format of time

        SELECT CASE (model)
        CASE ('WRF')
          PRINT *, "WRF time handling"

          sts = NF90_OPEN(iflWRFin, IOR(NF90_NOWRITE, NF90_MPIIO), ncid_in, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
          sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
          sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
          sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
            NAME = InDimNameRec, LEN = InDimLenRec)
          ! with the wrfxtrm minimum and maximum fields, the 1st field per file is
          ! empty, but the last field has already the next day assigned; this 
          ! leads to a plus 1 day time shift in the minimum and maximum data
          ! by shortening the input time vector and by offsetting the time-dim
          ! by 1 in the output writing this is fixed
          ! see "variable to read/write with no additional processing" part
          PRINT *, "number of timesteps in the input data: ", InDimLenRec
! KGo: not in my data
          IF ( ( cell_methods(ivar) == "minimum" ) .OR. &
               ( cell_methods(ivar) == "maximum" ) ) THEN
            InDimLenRec = InDimLenRec - 1
            PRINT *, "fixing number of input timesteps for min/max: ", InDimLenRec
          END IF
! HTr: also "point" cell_methods need to reduced by 1 record
!         IF ( ( cell_methods(ivar) == "minimum" ) .OR. &
!              ( cell_methods(ivar) == "maximum" ) .OR. &
!              ( cell_methods(ivar) == "point" ) ) THEN
!           InDimLenRec = InDimLenRec - 1
!           PRINT *, "fixing number of input timesteps: ", InDimLenRec
!         END IF
          ALLOCATE(InVarDataRec(InDimLenRec))
          sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
          sts = NF90_CLOSE(ncid_in)
          PRINT *, "RCM input file time coverage:"
          PRINT *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)

        CASE ('COS')
          PRINT *, "COSMO time handling"
          ! read using CDI, and put into the same variables -> 
          ! then the time matching needs not to be changed
        CASE DEFAULT
          PRINT *, "invalid model specified"
          STOP
          ! more model veriants are possible
        END SELECT

        ALLOCATE(InDateTimeYear(InDimLenRec))
        ALLOCATE(InDateTimeMonth(InDimLenRec))
        ALLOCATE(InDateTimeDay(InDimLenRec))
        ALLOCATE(InDateTimeHour(InDimLenRec))
        ALLOCATE(InDateTimeMinute(InDimLenRec))
        ALLOCATE(InDateTimeSecond(InDimLenRec))
        ALLOCATE(InDateTimeCombined(InDimLenRec))
        
        PRINT *, SIZE(InDateTimeYear)
        PRINT *, SIZE(InVarDataRec)
        PRINT *, InVarDataRec(:)
  
        ! this is the temporal coverage of the RCM input data of the specific
        ! file; split the date/time information string into several vectors
        ! Format example of WRF Times variable: 2000-10-06_00:00:00
        DO i = 1, SIZE(InVarDataRec), 1
        
          READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)' ) &
            InDateTimeYear(i), &
            InDateTimeMonth(i), &
            InDateTimeDay(i), &
            InDateTimeHour(i), &
            InDateTimeMinute(i), &
            InDateTimeSecond(i)
            
          ! 2009062008 use for later-on comparisons
          ! sub-hourly is not possible
          InDateTimeCombined(i) = & 
            InDateTimeYear(i) * 1000000 + &
            InDateTimeMonth(i) * 10000 + &
            InDateTimeDay(i) * 100 + &
            InDateTimeHour(i)
            
        END DO

!-------------------------------------------------------------------------------
! loop over the individual timesteps in the individual RCM input file
! this may take some time but it is robust and also extremely large arrays
! might be used; this leads to many I/O operations, but only small amounts of
! data, but the tradeoff is robustness, flexibility and still far less I/O 
! operations as if using script-based solutions (cdo, nco, etc.)

        DO it = 1, InDimLenRec, 1
  
! HTr: skip last record for "sum" or "mean" to avoid re-initialisation of next file with nans, if 
!      the next file is computed in parallel
! KGo: this does not work with my input data 00-23 per wrfout file
!      if inputtimesteptruncate=T, then there is a gap in the data, i.e., one timestep
!      at the end / beginning of new day remains missing
          IF (inputtimesteptruncate) THEN 
            IF ( ( it .EQ. InDimLenRec ) .AND. ( ( cell_methods(ivar) == "mean" ) .OR. &
                 ( cell_methods(ivar) == "sum" ) ) ) THEN
              PRINT *, "skip last input timestep to avoid NaNs in the following aggregation file ", InDimLenRec
              EXIT
            END IF
          END IF

          PRINT *, "-----------------------------------------------------------"
          PRINT *, "inner loop"
          PRINT *, "working on date in WRF inp. file = ", TRIM(InVarDataRec(it))
          PRINT *, "working on timestep it = ", &
                   InDateTimeYear(it), "-", InDateTimeMonth(it), "-", &
                   InDateTimeDay(it), "_", &
                   InDateTimeHour(it), " ", InDateTimeCombined(it)

!-------------------------------------------------------------------------------
! performance issue: speed up the tool during subsequent passes, i.e. is the 
! storage file etc. existing already
! check whether file exists, highly likely as the RCM outputs usually have a
! higher frequency (e.g., monthly or daily) than the temporal granularity of the
! storage files (e.g., annual with the CORDEX protocol)
! also check whether this is the first pass on a new storage file or a 
! subsequent pass
! because data is read with a very fine-grained granularity, it's highly likely,
! that the file which is supposed to hold the data exists already, e.g., 8760
! timesteps per year at 1h temporal resolution vs. creation of 1 annual file,
! the tool checks before any operation whether the same temporal aggregation has
! been dealt with before (yearly, monthly, individual storage)
! if the InDateTimeMonthPrev or InDateTimeYearPrev are > 0 all variables are set
! already also
! in case of annual and monthly aggregation, the file content of the WRF inputs
! defines the storage files creation, in case of the individual definition, the
! opposite is true

          IF (aggregation_individually) THEN ! first pass, always needed, file may exist or not
            IF ( prevpass == 0 ) THEN
              newpass = .TRUE.
            END IF
          ELSE IF (aggregation_monthly) THEN
            IF ( InDateTimeMonthPrev /= InDateTimeMonth(it) ) THEN
              newpass = .TRUE.
            END IF
          ELSE IF (aggregation_yearly) THEN ! default, if the current year is different than the one previous > create a new file
            IF ( InDateTimeYearPrev /= InDateTimeYear(it) ) THEN
              newpass = .TRUE.
            END IF
          END IF
          
          IF (newpass) THEN
  
            InDateTimeYearPrev = InDateTimeYear(it)
            InDateTimeMonthPrev = InDateTimeMonth(it)
            prevpass = 1
            newpass = .FALSE.
  
            PRINT *, "start of processing or new year/month/timespan encountered"

!-------------------------------------------------------------------------------
! extract the time info from the ref time array which matches the respective
! file in which data is to be written and thereby matches also the input data
! ...as there is no "WHERE" the way I need it in F95, use loops
! this is needed whenever a new netCDF file is to be used and also if
! this file exists already
! still working on single variable 'ivar', file 'ifl', and timestep 'it'
! InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it), InDateTimeCombined(it)
! TimeRefArray, TimeRefArraySubset
!
! the time span covered by the storage file is determining the time vector subset
! to be extracted from the reference time vec which is again spanning the
! complete experiment; once this reference time vector subset (matching the 
! storage file) is determined the individual offset may be determined which is 
! used to sort the individual WRF field in

            PRINT *, "subsetting the TimeRefArray"

            !PRINT *, "TimeRefArray: index, y, m, d, h"
            !PRINT *, "size & shape of the main TimRefArray = ", &
            !  SIZE(TimeRefArray), SHAPE(TimeRefArray)
            !PRINT "(6F20.6)", TRANSPOSE(TimeRefArray(20000:20010,:))
  
            DEALLOCATE( TimeRefArraySubset )
            DEALLOCATE( TimeRefArraySubsetMean )
            DEALLOCATE( Time_bnds )
            
!            PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
!            PRINT *,'SHAPE(TimeRefArray, 1)', SHAPE(TimeRefArray, 1)
!            PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
!            PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
!            PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)

            ! only needed for aggregation_individually
            READ( tsact, '(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)' ) &
              tsactYear, &
              tsactMonth, &
              tsactDay, &
              tsactHour, &
              tsactMinute, &
              tsactSecond
            READ( teact, '(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)' ) &
              teactYear, &
              teactMonth, &
              teactDay, &
              teactHour, &
              teactMinute, &
              teactSecond

            ! number of date/time steps in the ref array which match the current 
            ! year or month, or individual time span
            ! the TimeRefArray must span the requested subsets, this is a 
            ! special case as requested e.g. in some FPS sideline and WL/CM 
            ! studies
            ALLOCATE( ipos(0) )
            counter = 0
            IF (aggregation_individually) THEN
              PRINT *, "aggregation_individually"
              PRINT *,tsactYear, tsactMonth, tsactDay, tsactHour, teactYear, teactMonth, teactDay, teactHour

              tsact_singlenumber = REAL(tsactYear,8)*1000000._8 + REAL(tsactMonth,8)*10000._8 + &
                REAL(tsactDay,8)*100._8 + REAL(tsactHour,8)
              teact_singlenumber = REAL(teactYear,8)*1000000._8 + REAL(teactMonth,8)*10000._8 + &
                REAL(teactDay,8)*100._8 + REAL(teactHour,8)
              PRINT *, tsact_singlenumber, teact_singlenumber
              DO i = 1, SIZE(TimeRefArray, 1), 1
                IF ( ( TimeRefArray(i,6) >= tsact_singlenumber ) .AND. & 
                     ( TimeRefArray(i,6) <= teact_singlenumber ) ) THEN ! CHECK/TODO <= or < -->>> wrong wit <= when doing daily min/max: one day too many in time vec
                  counter = counter + 1
                  ipos = [ipos, i]
                END IF
              END DO
              PRINT *, "size ipos", SIZE(ipos)
              PRINT *, "counter = ", counter
              IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                   ( cell_methods(ivar) == "sum" ) .OR. &
                   ( cell_methods(ivar) == "minimum" ) .OR. &
                   ( cell_methods(ivar) == "maximum" ) ) THEN
                counter = counter - 1
                ipos = ipos(1:counter)
                PRINT *, "size ipos", SIZE(ipos)
                PRINT *, "counter = ", counter
              END IF     

!              DO i = 1, SIZE(TimeRefArray, 1), 1
!                IF ( ( TimeRefArray(i,2) >= tsactYear ) .AND. & 
!                     ( TimeRefArray(i,2) <= teactYear ) .AND. &
!                     ( TimeRefArray(i,3) >= tsactMonth ) .AND. &
!                     ( TimeRefArray(i,3) <= teactMonth ) .AND. &
!                     ( TimeRefArray(i,4) >= tsactDay ) .AND. &
!                     ( TimeRefArray(i,4) < teactDay ) ) THEN !.AND. &
!!                   ( TimeRefArray(i,5) >= tsactHour ) .AND. &
!!                   ( TimeRefArray(i,5) <= teactHour ) ) THEN ! problem if the end of also 00UTC > then there is no match with the hours
!                  ! special case, e.g. 20090620_00:00:00 to 20090627_00:00:00 = defined timespan
!                  ! the last mean field is then 20090627_00:30:00, one less mean fields in the total time span covered
!                  counter = counter + 1
!                  ipos = [ipos, i]
!                END IF
!              END DO
!              ! see the problem with the ending, there is one timestep too
!              ! few in the file if this is not done
!              ! artifically add a further postion at the very end
!                  IF ( cell_methods(ivar) == "point" ) THEN ! and then also automatically sub-daily
!                    ipos = [ipos, ipos(counter)+1] 
!                    counter = counter + 1
!                  END IF
            ELSE IF (aggregation_monthly) THEN
              PRINT *, "aggregation_monthly"
              DO i = 1, SIZE(TimeRefArray, 1), 1
                IF (( TimeRefArray(i,2) == InDateTimeYear(it)) .AND. &
                   ( TimeRefArray(i,3) == InDateTimeMonth(it))) THEN
                  counter = counter + 1
                  ipos = [ipos, i]
                END IF
              END DO
              !!ipos = [ipos, ipos(counter)+1] 
              !!counter = counter + 1
            ELSE IF (aggregation_yearly) THEN ! default, CORDEX annual files
              PRINT *, "aggregation_annually"
              PRINT *, InDateTimeYear(it)
              DO i = 1, SIZE(TimeRefArray, 1), 1
                IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
                  counter = counter + 1
                  ipos = [ipos, i]
                END IF
              END DO
              !!ipos = [ipos, ipos(counter)+1] ! XXXXXXXXXX remove this
              !!counter = counter + 1
            END IF
            PRINT *, "timesteps in the time ref. subset = ", counter

            ALLOCATE( TimeRefArraySubset( counter, 5 ) ) ! index, y, m, d, h
            ALLOCATE( TimeRefArraySubsetMean( counter ) )
            ALLOCATE( Time_bnds( 2, counter ) )

            ! the TimeRefArraySubset is the time vec of the newly created (or
            ! already existing) file; it may be any size =< the original ref
            ! time vec; WRF data will be matched with this time vec
            ! this time vec is basically created for point data
            DO i = 1, counter, 1
              j = ipos(i)
              TimeRefArraySubset(i,1:5) = TimeRefArray(j,1:5)
!              PRINT *, TimeRefArraySubset(i,1:5) ! index, y, m, d, h
            END DO
            !PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySubset(:,:) )
            !print *, TimeRefArraySubset(0,:)

            DEALLOCATE(ipos)

!-------------------------------------------------------------------------------
! create path- and filenames according to the ruleset of the CORDEX data 
! protocol: 2 main options here:
! Non standard:
! 1. individual time range
! 2. monthly, determined by input file
! Standard (storage file is determined by the temporal aggregation):
! 3. yearly, for anything higher than daily
! 4. 5-yearly or less, for anything daily
! 5. 10-yearly or less, for anything monthly or seasonally
! TODO: not all is implemented as the averaging is not yet implemented

            PRINT *, "generate path and filename"

            ! only a single file is created and filled successively
            ! use the information from the namelist
            ! 2009-06-20_00:00:00, may be full days or end e.g. at 23UTC
            ! arbitrary timespan
            IF (aggregation_individually) THEN

              WRITE (tsactYearStr,'(I4.4)') tsactYear
              WRITE (tsactMonthStr,'(I2.2)') tsactMonth
              WRITE (tsactDayStr,'(I2.2)') tsactDay
              WRITE (tsactHourStr,'(I2.2)') tsactHour
              WRITE (teactYearStr,'(I4.4)') teactYear
              WRITE (teactMonthStr,'(I2.2)') teactMonth
              WRITE (teactDayStr,'(I2.2)') teactDay
              WRITE (teactHourStr,'(I2.2)') teactHour 

              IF ((frequency(ifrq) == '1hr') .OR. &
                  (frequency(ifrq) == '3hr') .OR. &
                  (frequency(ifrq) == '6hr')) THEN

                FileNameStartDateTime = tsactYearStr//tsactMonthStr//tsactDayStr//tsactHourStr
                FileNameEndDateTime = teactYearStr//teactMonthStr//teactDayStr//teactHourStr

              ELSE IF (frequency(ifrq) == 'day') THEN
  
                FileNameStartDateTime = tsactYearStr//tsactMonthStr//tsactDayStr
                FileNameEndDateTime = teactYearStr//teactMonthStr//teactDayStr

              ELSE IF ((frequency(ifrq) == 'mon') .OR. &
                       (frequency(ifrq) == 'sem')) THEN

                FileNameStartDateTime = tsactYearStr//tsactMonthStr
                FileNameEndDateTime = teactYearStr//teactMonthStr

              END IF

            ! determined by the date information in the wrf outputs, automatic
            ! monthly aggregation is also special
            ELSE IF ( (aggregation_monthly) .OR. (aggregation_yearly) ) THEN

              !READ( InDateTimeYear(it), '(4A)' ) InDateTimeYearStr
              !READ( InDateTimeMonth(it), '(2A)' ) InDateTimeMonthStr
              WRITE (InDateTimeYearStr,'(I4.4)') InDateTimeYear(it)
              WRITE (InDateTimeMonthStr,'(I2.2)') InDateTimeMonth(it)              
              WRITE (LastDayStr,'(I2.2)') INT(TimeRefArraySubset(0,5))
              IF ( (cell_methods(ivar) == "mean") .OR. (cell_methods(ivar) == "sum") ) THEN 
                WRITE (FirstHourStr,'(I2.2)') INT( FLOOR( ((dtHours/2.)*60.) / 60. ) )
                WRITE (FirstMinuteStr,'(I2.2)') INT( MOD( (dtHours/2.)*60., 60. ) )
                WRITE (LastHourStr,'(I2.2)') INT( FLOOR( ( (24.*60.) - (dtHours/2.)*60.)  / 60. ) )
                WRITE (LastMinuteStr,'(I2.2)') INT( MOD( (24.*60.) - (dtHours/2.)*60., 60. ) )
                PRINT *, "first and last h + min strings: ", FirstHourStr, FirstMinuteStr, LastHourStr, LastMinuteStr
              ELSE
                FirstHourStr = "00"
                WRITE (LastHourStr,'(I2.2)') 24-INT(dtHours)
              END IF
              PRINT *, "date/time information for the automatic output path- and filename generation: ", &
                InDateTimeYearStr, InDateTimeMonthStr, FirstHourStr, LastHourStr, LastDayStr

              ! last hour may be YYYY-<NextMonth>-01_00:00:00 for point
              ! last hour may be YYYY-<NextMonth>-01_00:00:00 for mean (last 
              ! time in bound value) -> this is the old spec
              ! e.g.: time_bnds 21:00-24:00 and for the midpoint 22:30 
              ! here the new standard is applied:
              ! always have the midpoint time value of the first and last time 
              ! element in the file in the filename
              IF (aggregation_monthly) THEN

                IF ((frequency(ifrq) == '1hr') .OR. &
                    (frequency(ifrq) == '3hr') .OR. &
                    (frequency(ifrq) == '6hr')) THEN

                  IF ( (cell_methods(ivar) == "mean") .OR. (cell_methods(ivar) == "sum") ) THEN
                    FileNameStartDateTime = InDateTimeYearStr//InDateTimeMonthStr//"01"//FirstHourStr//FirstMinuteStr
                    FileNameEndDateTime = InDateTimeYearStr//InDateTimeMonthStr//LastDayStr//LastHourStr//LastMinuteStr
                  ELSE
                    FileNameStartDateTime = InDateTimeYearStr//InDateTimeMonthStr//"01"//FirstHourStr
                    FileNameEndDateTime = InDateTimeYearStr//InDateTimeMonthStr//LastDayStr//LastHourStr
                  END IF  

                ELSE IF (frequency(ifrq) == 'day') THEN
  
                  FileNameStartDateTime = InDateTimeYearStr//InDateTimeMonthStr//"01"
                  FileNameEndDateTime = InDateTimeYearStr//InDateTimeMonthStr//LastDayStr

                END IF

              ! with the full default, the aggregation in files depends all on 
              ! the temporal resolution of the data to be stored, see the 
              ! beginning of this section
              ! default: either annual, 5 years or 10 years, depending on the 
              ! aggregation
              ELSE IF (aggregation_yearly) THEN 

                IF ((frequency(ifrq) == '1hr') .OR. &
                    (frequency(ifrq) == '3hr') .OR. &
                    (frequency(ifrq) == '6hr')) THEN

                  IF ( (cell_methods(ivar) == "mean") .OR. (cell_methods(ivar) == "sum") ) THEN
                    FileNameStartDateTime = InDateTimeYearStr//"0101"//FirstHourStr//FirstMinuteStr
                    FileNameEndDateTime = InDateTimeYearStr//"1231"//LastHourStr//LastMinuteStr
                  ELSE
                    FileNameStartDateTime = InDateTimeYearStr//"0101"//FirstHourStr
                    FileNameEndDateTime = InDateTimeYearStr//"1231"//LastHourStr
                  END IF

                ! TODO: 5 yearly
                ELSE IF (frequency(ifrq) == 'day') THEN
  
                  FileNameStartDateTime = InDateTimeYearStr//"0101"
                  FileNameEndDateTime = InDateTimeYearStr//"1231"

                ! TODO: 10 yearly
                ELSE IF ((frequency(ifrq) == 'mon') .OR. &
                         (frequency(ifrq) == 'sem')) THEN

                  STOP "not yet implemented functionality, monthly/seasonally"
                !  FileNameStartDateTime = tsactYearStr//tsactMonthStr
                !  FileNameEndDateTime = teactYearStr//teactMonthStr

                END IF

              END IF

            END IF

            ! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/
            ! EUR-44/CRPGL/ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
            pn_out = TRIM(project_id)                    // "/" // &
                     TRIM(product)                       // "/" // &
                     TRIM(CORDEX_domain)                 // "/" // &
                     TRIM(institute_id)                  // "/" // &
                     TRIM(driving_model_id)              // "/" // &
                     TRIM(driving_experiment_name)       // "/" // &
                     TRIM(driving_model_ensemble_member) // "/" // &
                     TRIM(model_id)                      // "/" // &
                     TRIM(rcm_version_id)                // "/" // &
                     TRIM(frequency(ifrq))               // "/" // &
                     TRIM(var_cmip(ivar))

            ! evspsbl_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_CRPGL-WRFARW331_v1_
            ! 3hr_1989010100-1989123121
            fn_out = TRIM(var_cmip(ivar))                // "_" // &
                     TRIM(CORDEX_domain)                 // "_" // &
                     TRIM(driving_model_id)              // "_" // &
                     TRIM(driving_experiment_name)       // "_" // &
                     TRIM(driving_model_ensemble_member) // "_" // &
                     TRIM(model_id)                      // "_" // &
                     TRIM(rcm_version_id)                // "_" // &
                     TRIM(frequency(ifrq))               // "_" // &
                     TRIM(FileNameStartDateTime) // "-" // TRIM(FileNameEndDateTime) // &
                     ".nc"
  
            PRINT *, "CORDEX compliant pathname, pn_out = ", TRIM(pn_out)
            PRINT *, "CORDEX compliant filename, fn_out = ", TRIM(fn_out)

!-------------------------------------------------------------------------------
! check for existance of the new output file and generate this file if needed
! could exist already from a previous run of tool and due to multiple months in 
! a file (i.e. one WRF output file may cover several months)
! i.e. this could be the first pass of a processing session, but a subsequent 
! pass over the storage file when continuously adding new data 
! could exist already from a previous run of tool and due to multiple
! months in a file (i.e. 1 WRF output file may cover several months)

            PRINT *, "*** CHECK FOR FILE EXISTANCE / CREATE FILE WITH BASIC STRUCTURE + METADATA***"

            INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )

            IF ( FileExists ) THEN

              PRINT *, "++++ path and file exist, continue filling"

            ELSE

              PRINT *, "++++ path and file do not yet exist, create path and netCDF file first"
              PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)

              CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )

!-------------------------------------------------------------------------------
! the SYSTEM call is non-std Fortran95, works for gfortran (fct & subroutine) 
! and ifort
! comment lines in the netCDF file global attribute definition
! turn standard checking in Makefile off
! trackingID = "xxxxxxxx-xxxx-Mxxx-Nxxx-xxxxxxxxxxxx"
! creationDate = "YYYY-MM-DD-THH:MM:SSZ"
! with MPI, avoid same UUID / creation_date for all files and overlaps
! generate an individual tmp file per MPI rank (=var)

!            IF ( rank == 0 ) THEN
!              CALL SYSTEM(cmdUUID)
!            END IF
!
!            CALL mpi_barrier(MPI_COMM_WORLD, ierr)
!            IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
!              
!              OPEN(1,FILE="tmpfileUUID",STATUS='old')
!              READ(1,*) trackingID
!              CLOSE(1)
!              PRINT *, "uuidgen externally generated trackingID = ", trackingID
!
!            IF ( rank == 0 ) THEN
!              CALL SYSTEM(cmdDate)
!            END IF 
!
!            CALL mpi_barrier(MPI_COMM_WORLD, ierr)
!            IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
!              
!              OPEN(1,FILE="tmpfileDate",STATUS='old')
!              READ(1,*) creationDate
!              CLOSE(1)
!              PRINT *, "date externally generated creation date = ", creationDate

              WRITE(FileNrStr,'(i3)') rank

              CALL SYSTEM("uuidgen -t > tmpfileUUID"//TRIM(domain)//TRIM(ADJUSTL(FileNrStr)))
              OPEN(1,FILE="tmpfileUUID"//TRIM(domain)//TRIM(ADJUSTL(FileNrStr)),STATUS='old')
              READ(1,*) trackingID
              CLOSE(1)
              PRINT *, "uuidgen externally generated trackingID = ", trackingID

              CALL SYSTEM("date -u +%Y-%m-%d-T%H:%M:%SZ > tmpfileDate"//TRIM(domain)//TRIM(ADJUSTL(FileNrStr)))
              OPEN(1,FILE="tmpfileDate"//TRIM(domain)//TRIM(ADJUSTL(FileNrStr)),STATUS='old')
              READ(1,*) creationDate
              CLOSE(1)
              PRINT *, "date externally generated creation date = ", creationDate

              CALL mpi_barrier(MPI_COMM_WORLD, ierr)
              IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
             
!-------------------------------------------------------------------------------
! create netCDF file, must be NetCDF4 'classic data model' and compression lvl=1
! NF90_CLASSIC_MODEL = NetCDF4_classic
! NF90_HDF5 = NetCDF4 based on HDF5
! NF90_CLOBBER = old netCDF

              !https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/NF90_005fCREATE.html
              PRINT *, "create netCDF file"
              sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // &
                TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, &
                NF90_CLASSIC_MODEL), ncid)

              !-----------------------------------------------------------------

              ! always included define dimensions
              sts = NF90_DEF_DIM(ncid, "x", xfocus, x_dimid)
              sts = NF90_DEF_DIM(ncid, "y", yfocus, y_dimid)
              sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
              sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid)
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) <= 10. ) ) THEN
                sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
              END IF
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) > 10. ) ) THEN
                sts = NF90_DEF_DIM(ncid, "plev", 1, plev_dimid)
              END IF
              ! special 4D vars, i.e. with real depth dimension
              ! hardcoding number of depth layers in LSM
              ! with the nomenclature in the namelist this special case cannot
              ! be resolved easily, standard does not foresee many 3D vars
              IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. &
                  (TRIM(var_cmip(ivar)) == 'tsl')) THEN
                sts = NF90_DEF_DIM(ncid, "soil_layer", SIZE(DZShc), depth_dimid)
                sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
              ENDIF
              sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
              IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                   ( cell_methods(ivar) == "sum" ) .OR. &
                   ( cell_methods(ivar) == "minimum" ) .OR. &
                   ( cell_methods(ivar) == "maximum" ) ) THEN
                sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
              ENDIF

              !-----------------------------------------------------------------
              ! no rule in CF1.7 standard, for capital letters or not > make a
              ! rule of my own to have long names start with capital letter
              ! nice in plots

              ! always included -- longitude field, unrotated
              sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ x_dimid, y_dimid /), lon_varid)
              sts = nf90_def_var_deflate(ncid, lon_varid, 1, 1, 1)
              sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
              sts = nf90_put_att(ncid, lon_varid, "long_name", "Longitude")
              sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
              sts = nf90_put_att(ncid, lon_varid, "_CoordinateAxisType", "Lon") ! special addon, not needed, but allowed

              ! always included -- latitude field, unrotated
              sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ x_dimid, y_dimid /), lat_varid)
              sts = nf90_def_var_deflate(ncid, lat_varid, 1, 1, 1)
              sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
              sts = nf90_put_att(ncid, lat_varid, "long_name", "Latitude")
              sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
              sts = nf90_put_att(ncid, lat_varid, "_CoordinateAxisType", "Lat") ! special addon, not needed, but allowed
  
              ! always included -- longitude vector, rotated
              sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid)
              sts = nf90_def_var_deflate(ncid, rlon_varid, 1, 1, 1)
              sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
              sts = nf90_put_att(ncid, rlon_varid, "long_name", "Longitude in rotated pole grid")
              sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
              sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
  
              ! always included -- latitude vector, rotated
              sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid)
              sts = nf90_def_var_deflate(ncid, rlat_varid, 1, 1, 1)
              sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
              sts = nf90_put_att(ncid, rlat_varid, "long_name", "Latitude in rotated pole grid")
              sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
              sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")

              ! additional and useful
              ! important for remapping with cdo (conservative remapping only 
              ! possible with this information)
              !        vertices = 4 ;
              !        float lat_vertices(rlat, rlon, vertices) ;
              !                lat_vertices:units = "degrees_north" ;
              !        float lon_vertices(rlat, rlon, vertices) ;
              !                lon_vertices:units = "degrees_east" ;

              ! always included, restriction to one domain only
              sts = nf90_def_var(ncid, "rotated_pole", NF90_CHAR, rotated_pole_varid)
              sts = nf90_put_att(ncid, rotated_pole_varid, "long_name", "Coordinates of the rotated North Pole")
              sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
              sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", GeoNPLat)
              sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", GeoNPLon)
  
              ! depends whether height is set in the nml, all between 1.5 and 10m
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) <= 10. ) ) THEN
                sts = nf90_def_var(ncid, "height", NF90_DOUBLE, (/ height_dimid /), height_varid)
                sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
                sts = nf90_put_att(ncid, height_varid, "long_name", "Height")
                sts = nf90_put_att(ncid, height_varid, "units", "m")
                sts = nf90_put_att(ncid, height_varid, "positive", "up")
                sts = nf90_put_att(ncid, height_varid, "axis", "Z")
              END IF
  
              ! just level definition, single number, like height
              ! there is no distinction between height and plev in the nml file
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) > 10. ) ) THEN
                sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, (/ plev_dimid /), plev_varid)
                sts = nf90_put_att(ncid, plev_varid, "standard_name", "air_pressure")
                sts = nf90_put_att(ncid, plev_varid, "long_name", "Pressure")
                sts = nf90_put_att(ncid, plev_varid, "units", "Pa")
                sts = nf90_put_att(ncid, plev_varid, "positive", "down")
                sts = nf90_put_att(ncid, plev_varid, "axis", "Z")
                IF ( cell_methods(ivar) == "vmean" ) THEN ! if this is layers over which there has been some everaging
                  sts = nf90_put_att(ncid, plev_varid, "bounds", "plev_bnds")
                END IF
              END IF

              ! special case of 3D vars
              ! what we enter here is the thickness, but not the actual depth
              ! neither staggered nor upper and lower bounds
              ! as defined in ancilliary document
              IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. &
                  (TRIM(var_cmip(ivar)) == 'tsl')) THEN
                sts = nf90_def_var(ncid, "soil_layer", NF90_DOUBLE, (/ depth_dimid /), depth_varid)
                sts = nf90_put_att(ncid, depth_varid, "standard_name", "depth")
                sts = nf90_put_att(ncid, depth_varid, "long_name", "Soil layer depth")
                sts = nf90_put_att(ncid, depth_varid, "units", "m")
                sts = nf90_put_att(ncid, depth_varid, "positive", "down")
                sts = nf90_put_att(ncid, depth_varid, "axis", "Z")
                sts = nf90_put_att(ncid, depth_varid, "bounds", "soil_layer_bnds")

                ! 1st version has correct shape but it empty due to sorting
                ! 2nd version ha sincorrect shape and is filled form sorting
                ! -> fix the sorting
                sts = nf90_def_var(ncid, "soil_layer_bnds", NF90_DOUBLE, (/ nb2_dimid, depth_dimid /), soillayerbnds_varid)
                !sts = nf90_def_var(ncid, "soil_layer_bnds", NF90_DOUBLE, (/ depth_dimid, nb2_dimid/), soillayerbnds_varid)
              ENDIF

              ! for vertically averaged variables need the plev bounds
              ! plev_bnds(2), always just single field per file: this is just two numbers 
              IF ( cell_methods(ivar) == "vmean" ) THEN
                sts = nf90_def_var(ncid, "plev_bnds", NF90_DOUBLE, (/ nb2_dimid /), plevbnds_varid)
              END IF
  
              ! always included
              sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid)
              sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
              sts = nf90_put_att(ncid, rec_varid, "long_name", "Time")
              sts = nf90_put_att(ncid, rec_varid, "units", "days since " // tstot(1:10) // "T" // tstot(12:19) // "Z" )
              sts = nf90_put_att(ncid, rec_varid, "calendar", "standard")
              sts = nf90_put_att(ncid, rec_varid, "axis", "T")
              IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                   ( cell_methods(ivar) == "sum" ) .OR. &
                   ( cell_methods(ivar) == "minimum" ) .OR. &
                   ( cell_methods(ivar) == "maximum" ) ) THEN
                sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")
              END IF
  
              ! for mean variables need the time bounds
              ! no further attributes, like plev_bnds > confirmed
              IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                   ( cell_methods(ivar) == "sum" ) .OR. &
                   ( cell_methods(ivar) == "minimum" ) .OR. &
                   ( cell_methods(ivar) == "maximum" ) ) THEN
                sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid)
              END IF

              !-----------------------------------------------------------------
              
              ! always included -- global attributes
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "conventionsURL", conventionsURL)
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
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "rcm_version_id", rcm_version_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "CORDEX_domain", CORDEX_domain)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id", trackingID)
             
              ! optional global variables
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "title", title)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institute_run_id", institute_run_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "nesting_levels", nesting_levels)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment_nesting", comment_nesting)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment_1nest", comment_1nest)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment_2nest", comment_2nest)

              !-----------------------------------------------------------------
              ! always included -- definition of the individual variable
              ! compression is always on variable level, use compression with 
              ! level = 1 (=> CORDEX protocol),  
              ! REAL PERFORMANCE ISSUES: shuffling and chunking
              ! https://earthscience.stackexchange.com/questions/12527/regarding-compression-shuffle-filter-of-netcdf4
              ! https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f90.html#Variables
              ! https://www.unidata.ucar.edu/software/netcdf/workshops/2011/nc4chunking/index.html
              ! https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_perf_chunking.html
              ! compression on/off: 19MB->12MB, 0.1s->0.5s
              ! chunking: needs some careful considerations

              ! not sure whether the height/plev dimension this is needed
              ! with CORDEX only one layer always per variable, see the
              ! coordinates attribute which makes the association
              ! this is also better for coding and data reading
              !IF ( height(ivar) /= -999 ) THEN
              ! HTr: only 2nd option
              IF ((var_cmip(ivar) == "mrsol") &
                .OR. (var_cmip(ivar) == "tsl")) THEN
                sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, depth_dimid, rec_dimid /), x_varid)  
                !sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, depth_dimid, rec_dimid /), x_varid, chunksizes = (/10, 10, 1, 8/), shuffle = .TRUE., fletcher32 = .FALSE., endianness = nf90_endian_little, deflate_level = 1)
              ELSE
                sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid)
                !sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid, chunksizes = (/10, 10, 8/), shuffle = .TRUE., fletcher32 = .FALSE., endianness = nf90_endian_little, deflate_level = 1)
              END IF

              ! TODO -> determine chunksizes vector beforehand otherwise
              ! this can only make it worse
              ! see https://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_why_it_matters
              !sts = nf90_def_var_chunking(ncid, x_varid, NF90_CHUNKED, XXXchunksizesXXX) 
              ! fill up with missing values, despite unlimited dim.
              ! not needed, at this point the unlimited time dim is not filled
              !sts = nf90_def_var_fill(ncid, x_varid, 0, mv) 
              ! shuffle ON, deflate ON, deflate_level lowest
              sts = nf90_def_var_deflate(ncid, x_varid, 1, 1, 1)
              ! default, change with care
              !sts = nf90_def_var_endian(ncid, x_varid, NF90_ENDIAN_NATIVE)
                            
              sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
              sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
              sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
              IF ( positive(ivar) /= '-999' ) THEN
                sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
              END IF
              sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) <= 10. ) ) THEN
                sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat height")
              ELSE IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) > 10. ) ) THEN
                sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat plev") 
              ELSE IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. &
                (TRIM(var_cmip(ivar)) == 'tsl')) THEN
                !sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat soil_layer") 
                sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat") 
              ELSE
                sts = nf90_put_att(ncid, x_varid, "coordinates", "lon lat")
              END IF
              sts = nf90_put_att(ncid, x_varid, "grid_mapping", "rotated_pole")
              sts = nf90_put_att(ncid, x_varid, "missing_value", mv)
              sts = nf90_put_att(ncid, x_varid, "_FillValue", mv)

              !-----------------------------------------------------------------
  
              sts = NF90_ENDDEF(ncid)

              !-----------------------------------------------------------------  
              ! SUPER IMPORTANT: write time information, the complete time vec
              ! for the respective file
              ! this is a bit of a trick, although time informaiton in used to 
              ! sort in the data, the time vector is only derived from that 
              ! subset time vector, one could base the complete tool entirely
              ! very strictly on time and time_bnds information > v2 idea

              IF ( cell_methods(ivar) == "point" ) THEN
  
                PRINT *, 'cell_methods:', cell_methods(ivar)
                PRINT *, TimeRefArraySubset(:,1)
                sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySubset(:,1) )

              END IF

              IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                   ( cell_methods(ivar) == "sum" ) .OR. &
                   ( cell_methods(ivar) == "minimum" ) .OR. &
                   ( cell_methods(ivar) == "maximum" ) ) THEN

                ! problem
                ! https://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00071.html
                ! use the old scheme for the time_bnds
                ! mid-point of time interval [decimal days]
                TimeRefArraySubsetMean (:) = TimeRefArraySubset(:,1) + ( 0.5_8 * (1._8 / (24._8/dtHours) ) )
                PRINT *, "cell_methods:", cell_methods(ivar)
                PRINT *, "dtHours",  dtHours, dtHours/2.
                PRINT *, "TimeRefArraySubset", TimeRefArraySubset(:,1)
                PRINT *, "TimeRefArraySubsetMean, mid-points of time intervals: ", TimeRefArraySubsetMean(:)
                sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySubsetMean(:) )
                PRINT *, 'sts NF90_PUT_VAR time', sts
    
                Time_bnds(1,:) = TimeRefArraySubset(:,1)
                Time_bnds(2,:) = TimeRefArraySubset(:,1) + ( 1._8 * (1._8 / (24._8/dtHours) )) 
                PRINT *, "time bnds lower", Time_bnds(1,:)
                PRINT *, "time bnds upper", Time_bnds(2,:)
                sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
                PRINT *, 'sts NF90_PUT_VAR time bnds', sts
  
              END IF
              !print *,'TimeRefArraySubset(:,1)', TimeRefArraySubset(:,1)

              ! TODO
              ! define plev_bnds, needed mainly for clh clm cli

              !-----------------------------------------------------------------
              ! write coordinates and height/level information
              
              ! add non-rotated coordinate fields
              sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
                START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
              sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
                START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )

              ! add rotated coordinates
              sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
              sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
    
              ! add height from NML
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) <= 10. ) ) THEN
                sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
              END IF

              ! add plev from NML, [hPa] -> [Pa]
              IF ( ( height(ivar) /= -999 ) .AND. ( height(ivar) > 10. ) ) THEN
                sts = NF90_PUT_VAR(ncid, plev_varid, height(ivar)*100. )
              END IF

              ! total with Noah LSM 2m depth thickness of the layers in Noahi LSM / 0.1D0, 0.3D0, 0.6D0, 1.0D0 /
              ! add the center points of the depth layers, the bnds gives then the info on the actual thickness
              IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. (TRIM(var_cmip(ivar)) == 'tsl')) THEN
!                ALLOCATE(DZSmp(SIZE(DZShc)))
!                ALLOCATE(soillayer_bnds(SIZE(DZShc),2))
!                PRINT *, 'soil layer calc'
!                DO i=1,SIZE(DZShc)
!                  soillayer_bnds(1,i) = SUM(DZShc(1:i))-DZShc(i)
!                  soillayer_bnds(2,i) = SUM(DZShc(1:i))
!                  PRINT '(2(1X,F4.2))', soillayer_bnds(:,i)
!                  DZSmp(i) = soillayer_bnds(1,i) + DZShc(i)/2
!                  PRINT '(1X,F4.2)', DZSmp(i)
!                END DO
!                PRINT *, 'shape soillayer_bnds', SHAPE(soillayer_bnds)
!                ! depth (mid-point) of the layers, orientation: Z down
!                sts = NF90_PUT_VAR(ncid, soillayerbnds_varid, soillayer_bnds(:,:))
!                ! upper and lower bound of the layers
!                sts = NF90_PUT_VAR(ncid, depth_varid, DZSmp(:))
!                ! alternatively hardcode here
                !sts = NF90_PUT_VAR(ncid, depth_varid, (/ 0.1D0, 0.3D0, 0.6D0, 1.0D0 /) ) ! wrong
                sts = NF90_PUT_VAR(ncid, depth_varid, (/ 0.05D0, 0.25D0, 0.7D0, 1.5D0 /) ) ! center points of the depth layers, hardcoded
                sts = NF90_PUT_VAR(ncid, soillayerbnds_varid, RESHAPE((/ 0.0D0, 0.1D0, 0.1D0, 0.4D0, &
                  0.4D0, 1.0D0, 1.0D0, 2.0D0 /), (/2,4/))) ! bnds hardcoded
              END IF 

              !-----------------------------------------------------------------
  
              sts = NF90_CLOSE(ncid)

              !-----------------------------------------------------------------
  
            END IF ! file exists y/n

!-------------------------------------------------------------------------------
  
          ELSE ! checking whether this is the first pass or previously run

            PRINT *, "this is a subsequent pass, just read, procecess and sort in stuff"

          END IF

!-------------------------------------------------------------------------------
! match timestep 'it' of WRFin with the subset of the ref time vec which belongs
! to the netCDF file of the year/month/arbitrary timespan currently open to
! receive data
! counter = offset in the netCDF file
! this is done once for the 1h/3h/6h input (point/mean) AND daily data (min/max)
! this match 'translates' the way WRF stores variables tempotally to ESGF
! e.g. in the wrfxtrm files the min/max are for day1 are stored in 00UTC field
! of day2
  
          PRINT *, "reading WRF sim. res. = ", TRIM(InVarDataRec(it)), it
          !PRINT *, SIZE(TimeRefArraySubset,1)
          !PRINT *, SHAPE(TimeRefArraySubset)
          !PRINT *, "current transferred input time: ", InDateTimeYear(it), InDateTimeMonth(it), InDateTimeDay(it), InDateTimeHour(it)
  
          counter = 0
          time_match = .FALSE.
          DO i = 1, SIZE(TimeRefArraySubset,1),1 ! time content of the WRF file
  
            counter = counter + 1
  
            !PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)',TimeRefArraySubset(i,1),TimeRefArraySubset(i,2),TimeRefArraySubset(i,3),TimeRefArraySubset(i,4),TimeRefArraySubset(i,5)
  
            IF ( ( TimeRefArraySubset(i,2) == InDateTimeYear(it)  ) .AND. &
                 ( TimeRefArraySubset(i,3) == InDateTimeMonth(it) ) .AND. &
                 ( TimeRefArraySubset(i,4) == InDateTimeDay(it)   ) .AND. &
                 ( TimeRefArraySubset(i,5) == InDateTimeHour(it)  ) ) THEN
              time_match = .TRUE.
              EXIT
            END IF
  
          END DO

          ! in case the aggregation_individually is set and the desired file
          ! temporal coverage is after the beginning or before the end
          ! of the timespans covered by the WRF input file(s) this functionality
          ! prevents that the last counter value from above is just written to
          ! the file, even without a match; if there is no match then just go to
          ! the next 'it'
          ! with the annual and monthly this is not possible as with every 'it'
          ! there is a check whether a storage file exists
          IF (time_match) THEN
  
          PRINT *, "index, single number, where in the NC file the WRF data is sorted in = ", &
            counter
  
!-------------------------------------------------------------------------------
! read orig WRF outputs
! there is always a corresponding time-slot in the NC file
! extracted time from above
! "it" controls it all: timestep in the individual WRF file
! there is only one variable at a time under processing

          PRINT *, "*** SOME VARS ALWAYS HAVE TO BE READ: (T00), P00 ***"
  
          sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)
  
          ! HTr: read actual value of base state temperature T00
          ! vector, 1 value /timestep, usually constant
          ! default changed over time with WRF versions, was 290K, is 300K
          ! but for DA: checked with registry and namelist -> 290 is used
          T00(1) = 290.0 !300.0
          !sts = NF90_INQ_VARID(ncidin, "T00", t00_varid)
          !IF ( sts /= NF90_NOERR ) THEN
          !  T00(1) = 300.0
          !ELSE
          !  sts = NF90_GET_VAR(ncidin, t00_varid, T00(:), &
          !    START = (/ it /), COUNT = (/ 1 /) )
          !END IF
  
          ! HTr: use actual value of base state pressure P00, if possible
          sts = NF90_INQ_VARID(ncidin, "P00", p00_varid)
          IF ( sts /= NF90_NOERR ) THEN
            P00(1) = 100000.
          ELSE
            sts = NF90_GET_VAR(ncidin, p00_varid, P00(:), &
              START = (/ it /), COUNT = (/ 1 /) )
          END IF
          
          PRINT *, T00(1), P00(1)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! data handling here is inefficient
! read too many variables
! read all vertical levels

          PRINT *, "*** READING OF VARIABLES ***"
          PRINT *, "variable to work on = ", TRIM(var_cmip(ivar))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! psl [Pa] i Sea Level Pressure
! select via height:
!  ua [m s-1] i Eastward Wind
!  va [m s-1] i Northward Wind
!  wa
!  ta [K] i Air Temperature
!  hus [1] i Specific Humidity
!  zg
! prw [kg m-2] i Water Vapor Path
! clwvi [kg m-2] i Condensed Water Path  
! clwvi [kg m-2] i Condensed Water Path  
! cape [J kg-1] i 2-D Maximum convective available potential energy
! cin [J kg-1] i 2-D Maximum convective inhibition
! CONSIDER GRID STAGGERING
! keeping the code compact and avoiding unnecessary allocations is sometimes 
! tricky and may lead to unallocated variables during processing
! if changes are made here, it must be tested with each variable seperate (all 
! others OFF in the variable namelist) whether it works, otherwise the 
! allocation may still persist from a previous pass of the code with a different
! variable

          IF ( (var_cmip(ivar) == "psl") &
                .OR. (height(ivar) == 1000) &
                .OR. (height(ivar) == 925) &
                .OR. (height(ivar) == 850) &
                .OR. (height(ivar) == 700) &
                .OR. (height(ivar) == 500) &
                .OR. (height(ivar) == 200) &
                .OR. (var_cmip(ivar) == "prw") &
                .OR. (var_cmip(ivar) == "clwvi") &
                .OR. (var_cmip(ivar) == "clivi") &
                .OR. (var_cmip(ivar) == "cape") &
                .OR. (var_cmip(ivar) == "cin") ) THEN

            PRINT *,'read 3D vars'

            !-------------------------------------------------------------------
            ! internal vars needed in the pressure level / 3D processing section

            PRINT *, "allocate p_in, t_in, ph_fl" 
            IF (.not. ALLOCATED(p_in)) ALLOCATE( p_in( xfocus, yfocus, nz ), STAT=sts ) ! always needed
            IF (.not. ALLOCATED(t_in)) ALLOCATE( t_in( xfocus, yfocus, nz ), STAT=sts )
            IF (.not. ALLOCATED(ph_fl)) ALLOCATE( ph_fl( xfocus, yfocus, nz ), STAT=sts )

            IF ( var_cmip(ivar) == "prw" ) THEN
              IF (.not. ALLOCATED(prw)) ALLOCATE( prw( xfocus, yfocus ), STAT=sts )
            END IF

            IF ( var_cmip(ivar) == "psl" ) THEN
              IF (.not. ALLOCATED(psl_in)) ALLOCATE( psl_in ( xfocus, yfocus ), STAT=sts )
            END IF

            IF ( (var_cmip(ivar) == "cape" ) .OR. (var_cmip(ivar) == "cin" ) ) THEN
              IF (.not. ALLOCATED(cape)) ALLOCATE( cape( xfocus, yfocus ), STAT=sts )
              IF (.not. ALLOCATED(cin)) ALLOCATE( cin( xfocus, yfocus ), STAT=sts )
              IF (.not. ALLOCATED(t_p)) ALLOCATE( t_p( xfocus, yfocus, nz ), STAT=sts )
              IF (.not. ALLOCATED(qvs)) ALLOCATE( qvs( xfocus, yfocus, nz ), STAT=sts )
              IF (.not. ALLOCATED(lcl)) ALLOCATE( lcl( xfocus, yfocus ), STAT=sts )
              IF (.not. ALLOCATED(lfc)) ALLOCATE( lfc( xfocus, yfocus ), STAT=sts )
            END IF

            IF ( var_cmip(ivar) == "clwvi" ) THEN
              IF (.not. ALLOCATED(clwvi)) ALLOCATE( clwvi( xfocus, yfocus ), STAT=sts )
            ENDIF

            IF ( var_cmip(ivar) == "clivi" ) THEN
              IF (.not. ALLOCATED(clivi)) ALLOCATE( clivi( xfocus, yfocus ), STAT=sts )
            END IF

            ! PROBLEM to make this efficient: selection of vars via height
            ! is not really efficient; if one defines a list of concrete
            ! variable names this list has to always be expanded if a new 
            ! variable is needed
            IF ( height(ivar) > 10 ) THEN
              PRINT *, "prep. int. 3D pres. level vars"
              ! IF (.not. ALLOCATED(var_pl)) ALLOCATE( var_pl( xfocus, yfocus, 6 ), STAT=sts ) ! working on one variable at a time
              !!IF (.not. ALLOCATED(var_pl)) ALLOCATE( var_pl( xfocus, yfocus ), STAT=sts ) ! get completely rid of this thing
              ! var_pl(:,:,:) = mv
              !!var_pl(:,:) = mv
              IF (.not. ALLOCATED(var3d_in)) ALLOCATE( var3d_in( xfocus, yfocus, nz ), STAT=sts )
              !!!IF (.not. ALLOCATED(pout)) ALLOCATE( pout( SIZE(pout) ), STAT=sts ) ! # pressure levels can be removed
            END IF

            !-------------------------------------------------------------------
            ! read data as needed
            ! these 5 vars are always calculated

            IF ( ( var_cmip(ivar) == "prw" )   .OR. &
                 ( var_cmip(ivar) == "psl" )   .OR. &
                 ( height(ivar) > 10  )        .OR. &
                 ( var_cmip(ivar) == "clwvi" ) .OR. &
                 ( var_cmip(ivar) == "clivi" ) .OR. &
                 ( var_cmip(ivar) == "cape" )  .OR. &
                 ( var_cmip(ivar) == "cin" ) ) THEN

              PRINT *, "read P"
              IF (.not. ALLOCATED(pp_in)) ALLOCATE( pp_in( xfocus, yfocus, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "P", pp_varid)
              sts = NF90_GET_VAR(ncidin, pp_varid, pp_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )

              PRINT *, "read PB"
              IF (.not. ALLOCATED(pb_in)) ALLOCATE( pb_in( xfocus, yfocus, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "PB", pb_varid)
              sts = NF90_GET_VAR(ncidin, pb_varid, pb_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )

              PRINT *, "read PH"
              IF (.not. ALLOCATED(ph_in)) ALLOCATE( ph_in( xfocus, yfocus, nz+1 ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "PH", ph_varid)
              sts = NF90_GET_VAR(ncidin, ph_varid, ph_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz+1, 1 /) )

              PRINT *, "read PHB"
              IF (.not. ALLOCATED(phb_in)) ALLOCATE( phb_in( xfocus, yfocus, nz+1 ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "PHB", phb_varid)
              sts = NF90_GET_VAR(ncidin, phb_varid, phb_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz+1, 1 /) )

              PRINT *, "read T" ! perturbation potential temperature (theta-t0)
              IF (.not. ALLOCATED(theta_in)) ALLOCATE( theta_in( xfocus, yfocus, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "T", theta_varid)
              sts = NF90_GET_VAR(ncidin, theta_varid, theta_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )

            END IF

            !-------------------------------------------------------------------

            IF ( ( var_cmip(ivar) == "prw" ) .OR. ( var_cmip(ivar) == "psl" ) .OR. &
                 ( SCAN(var_cmip(ivar),"hus") == 1 ) .OR. ( var_cmip(ivar) == "cape" ) .OR. &
                 ( var_cmip(ivar) == "cin" ) ) THEN
              PRINT *, "read QVAPOR"
              IF (.not. ALLOCATED(qv_in)) ALLOCATE( qv_in( xfocus, yfocus, nz  ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "QVAPOR", qv_varid)
              sts = NF90_GET_VAR(ncidin, qv_varid, qv_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
            END IF  

            IF ( var_cmip(ivar) == "clwvi" ) THEN
              PRINT *, "read QCLOUD"
              IF (.not. ALLOCATED(qc_in)) ALLOCATE( qc_in( xfocus, yfocus, nz  ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "QCLOUD", qc_varid)
              sts = NF90_GET_VAR(ncidin, qc_varid, qc_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
            END IF

            IF ( ( var_cmip(ivar) == "clwvi" ) .OR. ( var_cmip(ivar) == "clivi" ) ) THEN
              PRINT *, "read QICE"
              IF (.not. ALLOCATED(qi_in)) ALLOCATE( qi_in( xfocus, yfocus, nz  ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "QICE", qi_varid)
              sts = NF90_GET_VAR(ncidin, qi_varid, qi_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
            END IF

            IF ( var_cmip(ivar) == "clwvi" ) THEN
              PRINT *, "read QRAIN"
              IF (.not. ALLOCATED(qr_in)) ALLOCATE( qr_in( xfocus, yfocus, nz  ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "QRAIN", qr_varid)
              sts = NF90_GET_VAR(ncidin, qr_varid, qr_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
            END IF
  
            IF ( ( var_cmip(ivar) == "clwvi" ) .OR. ( var_cmip(ivar) == "clivi" ) ) THEN
              PRINT *, "read QSNOW"
              IF (.not. ALLOCATED(qs_in)) ALLOCATE( qs_in( xfocus, yfocus, nz  ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "QSNOW", qs_varid)
              sts = NF90_GET_VAR(ncidin, qs_varid, qs_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
            END IF

            IF ( (SCAN(var_cmip(ivar),"ua") == 1) .OR. (SCAN(var_cmip(ivar),"va") == 1) ) THEN
              PRINT *, "read U"
              IF (.not. ALLOCATED(u_in)) ALLOCATE( u_in( xfocus+1, yfocus, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "U", u_varid)
              sts = NF90_GET_VAR(ncidin, u_varid, u_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus+1, yfocus, nz, 1 /) )
            END IF
  
            IF ( (SCAN(var_cmip(ivar),"ua") == 1) .OR. (SCAN(var_cmip(ivar),"va") == 1) ) THEN
              PRINT *, "read V"
              IF (.not. ALLOCATED(v_in)) ALLOCATE( v_in( xfocus, yfocus+1, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "V", v_varid)
              sts = NF90_GET_VAR(ncidin, v_varid, v_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus+1, nz, 1 /) )
            END IF
            
            IF ( (SCAN(var_cmip(ivar),"wa") == 1) ) THEN
              PRINT *, "read W"
              IF (.not. ALLOCATED(w_in)) ALLOCATE( w_in( xfocus, yfocus, nz+1 ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "W", w_varid)
              sts = NF90_GET_VAR(ncidin, w_varid, w_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz+1, 1 /) )
            END IF

            IF ( (SCAN(var_cmip(ivar),"ua") == 1) .OR. (SCAN(var_cmip(ivar),"va") == 1) ) THEN
              PRINT *, "read SINALPHA"
              IF (.not. ALLOCATED(sinalpha_in)) ALLOCATE( sinalpha_in( xfocus, yfocus ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "SINALPHA", sinalpha_varid)
              sts = NF90_GET_VAR(ncidin, sinalpha_varid, sinalpha_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            END IF
  
            IF ( (SCAN(var_cmip(ivar),"ua") == 1) .OR. (SCAN(var_cmip(ivar),"va") == 1) ) THEN
              PRINT *, "read COSALPHA"
              IF (.not. ALLOCATED(cosalpha_in)) ALLOCATE( cosalpha_in( xfocus, yfocus ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "COSALPHA", cosalpha_varid)
              sts = NF90_GET_VAR(ncidin, cosalpha_varid, cosalpha_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            END IF
  
            IF ( var_cmip(ivar) == "psl" ) THEN
              PRINT *, "read PSFC"
              IF (.not. ALLOCATED(psfc_in)) ALLOCATE( psfc_in( xfocus, yfocus ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "PSFC", psfc_varid)
              sts = NF90_GET_VAR(ncidin, psfc_varid, psfc_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! clt [%] a Total Cloud Fraction
  
          ELSE IF (var_cmip(ivar) == "clt") THEN

            IF (.not. ALLOCATED(cldfra_in)) ALLOCATE( cldfra_in( xfocus, yfocus, nz, 2 ), STAT=sts )   
  
            IF (it /= InDimLenRec) THEN

              sts = NF90_INQ_VARID(ncidin, "CLDFRA", varid)
              sts = NF90_GET_VAR(ncidin, varid, cldfra_in(:,:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), &
                COUNT = (/ xfocus, yfocus, nz, 2 /) )
  
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN

              sts = NF90_INQ_VARID(ncidin, "CLDFRA", varid)
              sts = NF90_GET_VAR(ncidin, varid, cldfra_in(:,:,:,1), &
                START = (/ xoffset, yoffset, 1, it /), &
                COUNT = (/ xfocus, yfocus, nz, 1 /) )
  
              iflWRFin = fl_input(ifl+1) ! set to the previous wrfoutfile 
                                          ! if it is not the first
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)

              sts = NF90_INQ_VARID(ncidin0, "CLDFRA", varid)
              sts = NF90_GET_VAR(ncidin0, varid, cldfra_in(:,:,:,2), &
                START = (/ xoffset, yoffset, 1, 1 /), &
                COUNT = (/ xfocus, yfocus, nz, 1 /) )
  
              sts = NF90_CLOSE(ncidin0)
  
              iflWRFin = fl_input(ifl)  !set to the current wrfout file again

            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.

            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! pr [kg m-2 s-1] a Precipitation 
! EURO-CORDEX Jan/2018 meeting, conv+incl. snow graupel, hail, etc.
! TODO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! always need two adjacent output intervals to calculate the amount of precip
! between dates and consider whether the bucket_mm has been used:
! http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3.8/users_guide_chap5.htm#bucket
  
          ELSE IF (var_cmip(ivar) == "pr") THEN 
             
            PRINT *, "read iflWRFin " , iflWRFin
            PRINT *, "it = ", it
            PRINT *, "inDimLenRec, number of output times in input = ", InDimLenRec

            ! e.g. it=1..24, 00UTC-23UTC, #24 fields/file > only 23 average fields
            ! 00-01UTC -> 00:30UTC average 
            ! try to read if 23UTC the 00UTC from next file
            ! this is all independent of bucket or not
            ! if 00-day1 to 00-day2, then the last search finds 00UTC-day2 as a match
            ! for 24UTC from day1 > will be written to 00:30 from day2 > no problem
            IF (it < InDimLenRec) THEN 

              PRINT *, "read rainc and rainnc two dates from the same file"

              IF (.not. ALLOCATED(rainnc_in)) ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
              sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )

              IF (.not. ALLOCATED(rainc_in)) ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)  
              sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )

              sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
              IF ( bucket_mm > 0. ) THEN
                IF (.not. ALLOCATED(i_rainnc_in)) ALLOCATE( i_rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "I_RAINNC", varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rainnc_in(:,:,:), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 2 /) )   
                IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "I_RAINC", varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rainc_in(:,:,:), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 2 /) )
              END IF

            ! if the successor file exists, get the data from that file
            ! get the last available field from the current input file
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN

              PRINT *, "read rainc and rainnc from the current and the subsequent file"

              IF (.not. ALLOCATED(rainnc_in)) ALLOCATE( rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )        
              sts = NF90_INQ_VARID(ncidin, "RAINNC", rainnc_varid)
              sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus,yfocus, 1 /) )

              IF (.not. ALLOCATED(rainc_in)) ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
              sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1 /) )

              sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
              IF ( bucket_mm > 0. ) THEN
                IF (.not. ALLOCATED(i_rainnc_in)) ALLOCATE( i_rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "I_RAINNC", varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rainnc_in(:,:,1), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 1 /) )   
                IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "I_RAINC", varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rainc_in(:,:,1), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 1 /) )
              END IF

              ! just get the next input file
              ! read the first timestep of the subsequent wrfout file
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
              sts = NF90_INQ_VARID(ncidin0, "RAINNC", rainnc_varid)
              sts = NF90_GET_VAR(ncidin0, rainnc_varid, rainnc_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT = (/ xfocus, yfocus,1 /) )

              sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)  
              sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT = (/ xfocus, yfocus,1 /) )

              sts = NF90_GET_ATT(ncidin0, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
              IF ( bucket_mm > 0. ) THEN
                sts = NF90_INQ_VARID(ncidin0, "I_RAINNC", varid)
                sts = NF90_GET_VAR(ncidin0, varid, i_rainnc_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT = (/ xfocus, yfocus, 1 /) )   
                sts = NF90_INQ_VARID(ncidin0, "I_RAINC", varid)
                sts = NF90_GET_VAR(ncidin0, varid, i_rainc_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT = (/ xfocus, yfocus, 1 /) )
              END IF

              sts = NF90_CLOSE(ncidin0)

              ! set to the current wrfout file again
              iflWRFin = fl_input(ifl) 

            ! no further file, nothing can be done
            ! just pass missing values on
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.
  
            END IF
 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prc [kg m-2 s-1] a Convective Precipitation 
! read two timesteps to calculate 1h or 3h sum  

          ELSE IF (var_cmip(ivar) == "prc") THEN

            IF (.not. ALLOCATED(rainc_in)) ALLOCATE( rainc_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
              sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2/) ) 

              sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
              IF ( bucket_mm > 0. ) THEN
                IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "I_RAINC", varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rainc_in(:,:,:), &
                  START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
              END IF

            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN 
  
                sts = NF90_INQ_VARID(ncidin, "RAINC", rainc_varid)
                sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_in(:,:,1), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 1/))

                sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
                IF ( bucket_mm > 0. ) THEN
                  IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
                  sts = NF90_INQ_VARID(ncidin, "I_RAINC", varid)
                  sts = NF90_GET_VAR(ncidin, varid, i_rainc_in(:,:,1), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
                END IF

                iflWRFin = fl_input(ifl+1)
  
                sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
                sts = NF90_INQ_VARID(ncidin0, "RAINC", rainc_varid)
                sts = NF90_GET_VAR(ncidin0, rainc_varid, rainc_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT = (/ xfocus,yfocus,1 /) ) 

                sts = NF90_GET_ATT(ncidin0, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
                IF ( bucket_mm > 0. ) THEN
                  sts = NF90_INQ_VARID(ncidin0, "I_RAINC", varid)
                  sts = NF90_GET_VAR(ncidin0, varid, i_rainc_in(:,:,2), &
                    START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
                END IF
   
                sts = NF90_CLOSE(ncidin0)
  
                iflWRFin = fl_input(ifl)  !set to the current wrfout file again

            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prsn [kg m-2 s-1] ?>a Snowfall Flux
  
          ELSE IF (var_cmip(ivar) == "prsn") THEN
  
            IF (.not. ALLOCATED(snownc_in)) ALLOCATE( snownc_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
              sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
           
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SNOWNC", snownc_varid)
              sts = NF90_GET_VAR(ncidin, snownc_varid, snownc_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1/))
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
              sts = NF90_INQ_VARID(ncidin0, "SNOWNC", snownc_varid)
              sts = NF90_GET_VAR(ncidin0, snownc_varid, snownc_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT = (/xfocus, yfocus, 1 /) )
  
              sts = NF90_CLOSE(ncidin0)
  
              iflWRFin = fl_input(ifl)
  
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! snm [kg m-2 s-1] a Surface Snow Melt
  
          ELSE IF (var_cmip(ivar) == "snm") THEN
  
            IF (.not. ALLOCATED(acsnom_in)) ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
              sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
  
              print*, sts
  
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "ACSNOM", acsnom_varid)
              sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1 /))
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
              
              sts = NF90_INQ_VARID(ncidin0, "ACSNOM", acsnom_varid)
              sts = NF90_GET_VAR(ncidin0, acsnom_varid, acsnom_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT =(/xfocus, yfocus, 1 /) )
            
              sts = NF90_CLOSE(ncidin0)
            
              iflWRFin = fl_input(ifl)
              
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! evspsbl [kg m-2 s-1] a Evaporation
  
          ELSE IF (var_cmip(ivar) == "evspsbl") THEN
  
            IF (.not. ALLOCATED(sfcevp_in)) ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
              sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
            
            ELSE IF ( (it == InDimLenRec) .AND. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SFCEVP", sfcevp_varid)
              sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1/))
  
              iflWRFin = fl_input(ifl+1) 
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
              sts = NF90_INQ_VARID(ncidin0, "SFCEVP", sfcevp_varid)
              sts = NF90_GET_VAR(ncidin0, sfcevp_varid, sfcevp_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT =(/xfocus,yfocus, 1 /) )
  
              sts = NF90_CLOSE(ncidin0)
  
              iflWRFin = fl_input(ifl)  !set to the current wrfout file again
  
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF 
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! evspblpot [kg m-2 s-1] a Potential Evapotranspiration 
  
          ELSE IF (var_cmip(ivar) == "evspsblpot") THEN
  
            IF (.not. ALLOCATED(potevp_in)) ALLOCATE( potevp_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN         
  
              sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
              sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
  
            ELSE IF ( (it == InDimLenRec) .AND. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "POTEVP", potevp_varid)
              sts = NF90_GET_VAR(ncidin, potevp_varid, potevp_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus,1/))
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
                sts = NF90_INQ_VARID(ncidin0, "POTEVP", potevp_varid)
                sts = NF90_GET_VAR(ncidin0, potevp_varid, potevp_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT =(/xfocus,yfocus,1 /) )
              sts = NF90_CLOSE(ncidin0)
  
              iflWRFin = fl_input(ifl)  !set to the current wrfout file again
  
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF 
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrros [kg m-2 s-1] a Surface Runoff
  
          ELSE IF (var_cmip(ivar) == "mrros") THEN
  
            IF (.not. ALLOCATED(sfroff_in)) ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
  
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus,1/))
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
                sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
                sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT =(/xfocus,yfocus, 1 /) )
              sts = NF90_CLOSE(ncidin0)

              iflWRFin = fl_input(ifl)
  
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrro [kg m-2 s-1] a Total Runoff
  
          ELSE IF (var_cmip(ivar) == "mrro") THEN
  
            IF (.not. ALLOCATED(sfroff_in)) ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
            IF (.not. ALLOCATED(udroff_in)) ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
  
              sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
              sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
  
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "SFROFF", sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus,1/))
  
              sts = NF90_INQ_VARID(ncidin, "UDROFF", udroff_varid)
              sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT =(/xfocus,yfocus,1 /) )
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
              sts = NF90_INQ_VARID(ncidin0, "SFROFF", sfroff_varid)
              sts = NF90_INQ_VARID(ncidin0, "UDROFF", udroff_varid)
  
              sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT =(/xfocus, yfocus, 1 /) )
  
              sts = NF90_GET_VAR(ncidin0, udroff_varid, udroff_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT =(/xfocus, yfocus, 1 /) )
  
              sts = NF90_CLOSE(ncidin0)
  
              iflWRFin = fl_input(ifl)
  
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.  
  
            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! hfss [W m-2] a Surface Upward Latent Heat Flux
! hfls [W m-2] a Surface Upward Sensible Heat Flux
! this was with the radiation before until v05 including; but for hfss and hlfs 
! there is no bucket defined, i.e. BUCKET_J is not applicable and also the 
! calculation below was wrong: leads to gross mistake
  
          ELSE IF ((var_cmip(ivar) == "hfss") &
            .OR. (var_cmip(ivar) == "hfls")) THEN

            PRINT *, "get heat flux variables"

            IF (.not. ALLOCATED(hflux_in)) ALLOCATE( hflux_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
 
              PRINT *, "read hflux_in for two times from the same file"

              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, hflux_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
            
            ELSE IF ( (it == InDimLenRec) .AND. (ifl /= SIZE(fl_input)) ) THEN

              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, hflux_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1/))

              iflWRFin = fl_input(ifl+1)

              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)

              sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin0, varid, hflux_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT = (/ xfocus, yfocus, 1/))

              sts = NF90_CLOSE(ncidin0)

              iflWRFin = fl_input(ifl)

            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.
  
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! rsds [W m-2] a Surface Downwelling Shortwave Radiation
! rlds [W m-2] a Surface Downwelling Longwave Radiation
! rsus [W m-2] a Surface Upwelling Shortwave Radiation 
! rlus [W m-2] a Surface Upwelling Longwave Radiation
! rlut [W m-2] a TOA Outgoing Longwave Radiation
! rsdt [W m-2] a TOA Incident Shortwave Radiation
! rsut [W m-2] a TOA Outgoing Shortwave Radiation
! only method in use finally: accumulated values with bucket
! restriction: assume bucket is used
  
          ELSE IF ((var_cmip(ivar) == "rsds") &
            .OR. (var_cmip(ivar) == "rlds") &
            .OR. (var_cmip(ivar) == "rsus") &
            .OR. (var_cmip(ivar) == "rlus") &
            .OR. (var_cmip(ivar) == "rlut") &
            .OR. (var_cmip(ivar) == "rsdt") &
            .OR. (var_cmip(ivar) == "rsut")) THEN

            PRINT *, "get radiation variables incl. bucket if set"

            IF (.not. ALLOCATED(rad_in)) ALLOCATE( rad_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
 
              PRINT *, "read rad_in for two times from the same file"

              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
 
              sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_J", bucket_J)
              IF ( bucket_J > 0. ) THEN
                PRINT *, "BUCKET_J = ", bucket_J
                IF (.not. ALLOCATED(i_rad_in)) ALLOCATE( i_rad_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, TRIM('I_' // var_wrf(ivar)), varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rad_in(:,:,:), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 2 /) )
              END IF

            ELSE IF ( (it == InDimLenRec) .AND. (ifl /= SIZE(fl_input)) ) THEN
  
              PRINT *, "read rad_in for one timestep from two files each"

              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, rad_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1/))

              sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_J", bucket_J)
              IF ( bucket_J > 0. ) THEN
                PRINT *, "BUCKET_J = ", bucket_J
                IF (.not. ALLOCATED(i_rad_in)) ALLOCATE( i_rad_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, TRIM('I_' // var_wrf(ivar)), varid)
                sts = NF90_GET_VAR(ncidin, varid, i_rad_in(:,:,1), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 1 /) )
              END IF

              iflWRFin = fl_input(ifl+1)

              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)

              sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin0, varid, rad_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT = (/ xfocus, yfocus, 1/))

              sts = NF90_GET_ATT(ncidin0, NF90_GLOBAL, "BUCKET_J", bucket_J)
              IF ( bucket_J > 0. ) THEN
                sts = NF90_INQ_VARID(ncidin0, TRIM('I_' // var_wrf(ivar)), varid)
                sts = NF90_GET_VAR(ncidin0, varid, i_rad_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT = (/ xfocus, yfocus, 1 /) )
              END IF

              sts = NF90_CLOSE(ncidin0)

              iflWRFin = fl_input(ifl)
  
            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.
  
            END IF
  
            PRINT *, "values in middle of domain = ", var_cmip(ivar), rad_in(xfocus/2,yfocus/2,1), rad_in(xfocus/2,yfocus/2,2)
            PRINT *, "difference [J m-2] = ", (rad_in(xfocus/2,yfocus/2,2) - rad_in(xfocus/2,yfocus/2,1))
            PRINT *, "mean [W m-2] = ", (rad_in(xfocus/2,yfocus/2,2) - rad_in(xfocus/2,yfocus/2,1))/ (dtHours*3600.)
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sund [s] s Duration of Sunshine
! WMO definition >= 120 W m-2

          ELSE IF (var_cmip(ivar) == "sund") THEN

            IF (.not. ALLOCATED(swdown_in)) ALLOCATE( swdown_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
 
              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, swdown_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )

            ELSE IF ( (it == InDimLenRec) .AND. (ifl /= SIZE(fl_input)) ) THEN

              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, swdown_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1/))

              iflWRFin = fl_input(ifl+1)

              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)

              sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin0, varid, swdown_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT = (/ xfocus, yfocus, 1/))

              sts = NF90_CLOSE(ncidin0)

              iflWRFin = fl_input(ifl)

            ELSE

              PRINT *, "no data available for average calculation any more"

              calc = .FALSE.
  
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrso [kg m-2] i Total Soil Moisture Content
! was beforehand treated as averaged variable, but is defined instantaneous
! new, also using smois_in, would be better soil_in
! mrsol [kg m-2] i Water Content of Soil Layer
! tsl [K] i Temperature of Soil
  
          ELSE IF ((var_cmip(ivar) == "mrso") &
            .OR. (var_cmip(ivar) == "mrsol") &
            .OR. (var_cmip(ivar) == "tsl")) THEN
 
!            IF (.not. ALLOCATED(DZS)) ALLOCATE( DZS( 4 ), STAT=sts )
!            sts = NF90_INQ_VARID(ncidin, "DZS", varid)
!            IF ( sts /= NF90_NOERR ) THEN
!              DZS = (/ 0.1, 0.3, 0.6, 1.0 /) ! total with Noah LSM 2m depth
!              PRINT*,'attention: using hardwired layer thickness for Noah LSM'
!            ELSE
!              sts = NF90_GET_VAR(ncidin, varid, DZS, &
!              START = (/ 1, it /), COUNT = (/ 4, 1 /) )
!              PRINT*,'reading layer thickness from file'
!            END IF
!           use hardcoded
            IF (.not. ALLOCATED(DZS)) ALLOCATE( DZS( SIZE(DZShc) ), STAT=sts )
            DZS(:) = DZShc(:)
            PRINT*,'DZS ', DZS(:)
 
!            IF (.not. ALLOCATED(smois_in)) ALLOCATE( smois_in( xfocus, yfocus, 4, 2 ), STAT=sts )
            IF (.not. ALLOCATED(smois_in)) ALLOCATE( smois_in( xfocus, yfocus, 4, 1 ), STAT=sts )
  
!            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), &
                COUNT = (/ xfocus, yfocus, SIZE(DZS), 1 /) )
!                COUNT = (/ xfocus, yfocus, 4, 2 /) )
  
!            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
!              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
!              sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:,1), &
!                START = (/ xoffset, yoffset, it /), &
!                COUNT = (/ xfocus, yfocus,1/))
  
!              iflWRFin = fl_input(ifl) ! set to the previous wrfoutfile 
!                                        ! if it is not the first 
  
!              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
!              sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
!              sts = NF90_GET_VAR(ncidin0, varid, smois_in(:,:,:,2), &
!                START = (/ xoffset, yoffset, 1 /), &
!                COUNT =(/xfocus,yfocus,1 /) )
  
!              sts = NF90_CLOSE(ncidin0)
  
!              iflWRFin = fl_input(ifl)  !set to the current wrfout file again
  
!            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! hurs [%] i Near-Surface Relative Humidity

          ELSE IF (var_cmip(ivar) == "hurs") THEN

            IF (.not. ALLOCATED(t2_in)) ALLOCATE( t2_in ( xfocus, yfocus ), STAT=sts )
            sts = NF90_INQ_VARID(ncidin, "T2", t2_varid)
            sts = NF90_GET_VAR(ncidin, t2_varid, t2_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            IF (.not. ALLOCATED(q2_in)) ALLOCATE( q2_in ( xfocus, yfocus ), STAT=sts )
            sts = NF90_INQ_VARID(ncidin, "Q2", q2_varid)
            sts = NF90_GET_VAR(ncidin, q2_varid, q2_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            IF (.not. ALLOCATED(psfc_in)) ALLOCATE( psfc_in ( xfocus, yfocus ), STAT=sts )
            sts = NF90_INQ_VARID(ncidin, "PSFC", psfc_varid)
            sts = NF90_GET_VAR(ncidin, psfc_varid, psfc_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sfcWind [m s-1] i Near-Surface Wind Speed

          ELSE IF (var_cmip(ivar) == "sfcWind") THEN
  
            IF (.not. ALLOCATED(u10_in)) ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts ) 
            IF (.not. ALLOCATED(v10_in)) ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts )
  
            sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
            sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)
  
            sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! uas [m s-1] i Eastward Near-Surface Wind
! vas [m s-1] i Northward Near-Surface Wind

          ELSE IF ((var_cmip(ivar) == "uas") .OR. (var_cmip(ivar) == "vas")) THEN
  
            IF (.not. ALLOCATED(u10_in)) ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts )
            IF (.not. ALLOCATED(v10_in)) ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts )
            IF (.not. ALLOCATED(sinalpha_in)) ALLOCATE( sinalpha_in( xfocus, yfocus ), STAT=sts )
            IF (.not. ALLOCATED(cosalpha_in)) ALLOCATE( cosalpha_in( xfocus, yfocus ), STAT=sts )
  
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
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prhmax [kg m-2 s-1] m Daily Maximum Hourly Precipitation Rate (using wrfxtrm)
! reuse variables form rainc and rainnc
! special offset, as needed with minimum and maximum per file 
! as 'it' counts per input file this is OK

          ELSE IF (var_cmip(ivar) == "prhmax") THEN
 
            IF (.not. ALLOCATED(rainc_max_in)) ALLOCATE( rainc_max_in ( xfocus, yfocus ), STAT=sts ) 
            IF (.not. ALLOCATED(rainnc_max_in)) ALLOCATE( rainnc_max_in ( xfocus, yfocus ), STAT=sts )
  
            sts = NF90_INQ_VARID(ncidin, "RAINCVMAX", rainc_varid)
            sts = NF90_INQ_VARID(ncidin, "RAINNCVMAX", rainnc_varid)
  
            sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_max_in(:,:), &
              START = (/ xoffset, yoffset, it+1 /), COUNT = (/ xfocus, yfocus, 1 /) )
            sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_max_in(:,:), &
              START = (/ xoffset, yoffset, it+1 /), COUNT = (/ xfocus, yfocus, 1 /) )
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! two possibilities with sea ice: either in SEAICE variable as fractional sea 
! ice, or contained as binary sea ice in the landmask field during runtime
! extent: XLAND = LU_INDEX

!          ELSE IF ( (var_cmip(ivar) == "sic") .AND. (fractSeaIce .EQV. .FALSE.) ) THEN

!            IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in ( xfocus, yfocus ), STAT=sts ) 
!            IF (.not. ALLOCATED(xland_in)) ALLOCATE( xland_in ( xfocus, yfocus ), STAT=sts )
  
!            sts = NF90_INQ_VARID(ncidin, "LANDMASK", landmask_varid)
!            sts = NF90_INQ_VARID(ncidin, "XLAND", xland_varid)
  
!            sts = NF90_GET_VAR(ncidin, landmask_varid, landmask_in(:,:), &
!              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
!            sts = NF90_GET_VAR(ncidin, xland_varid, xland_in(:,:), &
!              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! this is all other variables which are solely based on a namelist for whom 
! no processing is needed whatsoever, e.g. tas, ps, huss, ...
! they are read into data_in and also passed on in this 2D array for writing

          ELSE

            PRINT *, "variable to read/write with no additional processing = ", var_wrf(ivar)
  
            sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
 
!  sts = nf90_inquire_variable(ncidin, varid, contiguous = contiguous,  &
!       deflate_level = deflate_level, shuffle = shuffle, &
!fletcher32 = fletcher32, endianness = endianness)
!
!PRINT *, xtype, ndims, natts, contiguous, deflate_level, shuffle, fletcher32, endianness
!  
  
            IF ( ( cell_methods(ivar) == "minimum" ) .OR. &
                 ( cell_methods(ivar) == "maximum" ) ) THEN
              sts = NF90_GET_VAR(ncidin, varid, data_in(:,:), &
                    START = (/ xoffset, yoffset, it+1 /), COUNT = (/ xfocus, yfocus, 1 /) )
            ELSE 
              sts = NF90_GET_VAR(ncidin, varid, data_in(:,:), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            END IF

!            PRINT *, 'some sample output 3x3 in the middle of the domain', &
!              data_in(xfocus/2:(xfocus/2+2),yfocus/2:(yfocus/2+2))

          END IF
  
          sts = NF90_CLOSE(ncidin)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! some analysis of the data
!  
!         PRINT *, "*** STATISTICS BEFORE PROCESSING OF VARIABLES ***"
!         PRINT *, "(useless if data is stored in other vars than 'data_in')"
!
!         print *, "shape of array" , SHAPE(data_in)
!         print *, "size of array" , SIZE(data_in)
!         stat_mean = SUM(data_in(:,:))/SIZE(data_in(:,:))
!         PRINT *, "mean of array", stat_mean
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! this is where the real processing takes place 
! if nothing is to be calculated or scaled, etc., then the variables are just
! passed on to the write section in data_in
  
          PRINT *, "*** PROCESSING OF VARIABLES ***"

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! vars on pressure levels > for simplicity and as we like to use these data
! instead of the original model outputs: extract and convert more than specified
! at full temporal resolution, selective aggregation later
! use different interpolation methods for different variables
! psl [Pa] i Sea Level Pressure
! ua [m s-1] i Eastward Wind
! va [m s-1] i Northward Wind
! ta [K] i Air Temperature
! hus [1] i Specific Humidity
! also needed:
! wa
! zg
! levels needed CORDEX:            850       500  200
! levels needed FPS   : 1000  925   "   700   "    "

          IF ( (var_cmip(ivar) == "psl") .OR. &
               (height(ivar) == 1000) .OR. &
               (height(ivar) == 925) .OR. &
               (height(ivar) == 850) .OR. &
               (height(ivar) == 700) .OR. &
               (height(ivar) == 500) .OR. &
               (height(ivar) == 200) .OR. &
               (var_cmip(ivar) == "prw") .OR. &
               (var_cmip(ivar) == "clwvi") .OR. &
               (var_cmip(ivar) == "clivi") .OR. &
               (var_cmip(ivar) == "cape") .OR. &
               (var_cmip(ivar) == "cin") ) THEN

            PRINT *, var_cmip(ivar), height(ivar)

            ! needs this initialisation, might be possible that for very low
            ! high levels calculation cannot be done, then it is set to official
            ! missing value
            data_in(:,:) = mv

            ! total pressure [Pa]
            ! perturbation pressure + base state pressure
            ! see: http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3.8/users_guide_chap5.htm
            p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)

            ! PH and PHB are on nz+1 levels
            ! vertically destaggering
            ! total geopotential height [m]
            ! see: http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3.8/users_guide_chap5.htm
            DO nl = 1,nz 
              ph_fl(:,:,nl) = ((ph_in(:,:,nl)+phb_in(:,:,nl))+ &
                              (ph_in(:,:,nl+1)+phb_in(:,:,nl+1)))/2./gr
            END DO

            ! transfer theta-t0 to total potential temperature [K]
            ! and then convert potential temperature theta to absolute temperature
            ! this is the correct version, using 290K as base temp
            ! wrfpress contains data that are still calculated with 300K
            t_in(:,:,:) = ( theta_in(:,:,:) + T00(1) ) * &
                          ( p_in(:,:,:) / P00(1) )**(R/cp)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            IF (height(ivar) == 1000) THEN
              np = 1
            ELSE IF (height(ivar) == 925) THEN
              np = 2
            ELSE IF (height(ivar) == 850) THEN
              np = 3
            ELSE IF (height(ivar) == 700) THEN
              np = 4
            ELSE IF (height(ivar) == 500) THEN
              np = 5
            ELSE IF (height(ivar) == 200) THEN
              np = 6
            END IF
 
            ! vars needed: 3x int.: t_in, p_in, ph_fl 
            !              -> need: ph_in, phb_in, theta_in, pp_in, pb_in
            ! var_pl, var3d_in
            ! pout, slope, zg_pout
            ! comparison with wrfpress, same time: distribution OK, but ours is
            ! much colder, old vs new scheme about identical, ta500
            IF ( (var_cmip(ivar) == "ta1000") .OR. &
                 (var_cmip(ivar) == "ta925") .OR.  &
                 (var_cmip(ivar) == "ta850") .OR.  &
                 (var_cmip(ivar) == "ta700") .OR.  &
                 (var_cmip(ivar) == "ta500") .OR.  &
                 (var_cmip(ivar) == "ta200") ) THEN
              PRINT *, "calc ta..."
              var3d_in(:,:,:) = t_in(:,:,:)
              ! linear in p
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN
                      ! HTr: calculate zg at pout, first (linear in log(p)
                      !slope = (ph_fl(i,j,nl)-ph_fl(i,j,nl+1))/ (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                      !zg_pout = ph_fl(i,j,nl+1) + slope* (LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      ! HTr: interpolate linearly in zg
                      !slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ &
                      !        (ph_fl(i,j,nl)-ph_fl(i,j,nl+1))
                      !var_pl(i,j) = var3d_in(i,j,nl+1) + &
                      !              slope*(zg_pout-ph_fl(i,j,nl+1))
                      slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                              (p_in(i,j,nl)-p_in(i,j,nl+1))
                      data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO

            ! vars needed: int.: p_in (> in: pp_in, pb_in)
            !              out: var_pl, var3d_in
            !              tmp: pout, slope, zg_pout
            !              in: qv_in
            ! comparison with wrfpress, same time: about identical, also old and 
            ! new scheme, simple vs. sophisticated, hus500
            ELSE IF ( (var_cmip(ivar) == "hus1000") .OR. &
                 (var_cmip(ivar) == "hus925") .OR.       &
                 (var_cmip(ivar) == "hus850") .OR.       &
                 (var_cmip(ivar) == "hus700") .OR.       &
                 (var_cmip(ivar) == "hus500") .OR.       &
                 (var_cmip(ivar) == "hus200") ) THEN
              PRINT *, "calc hus..."
              var3d_in(:,:,:) = qv_in(:,:,:)
              ! either linear in p or log(p)
              !$OMP PARALLEL DO
              DO i = 1,xfocus
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN
                      !slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ &
                      !        (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                      !var_pl(i,j) = var3d_in(i,j,nl+1) + &
                      !              slope*(LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      SELECT CASE (vint_type)
                      CASE(0)
                      ! (HTr) - linear in p
                        slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                                (p_in(i,j,nl)-p_in(i,j,nl+1))
                        data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                      CASE(1)
                      ! (HTr) - linear in log(p)
                        slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                                (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                        data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      CASE DEFAULT
                        PRINT *, 'CAUTION: unknown setting for vint_type'
                        STOP 'vint_type not properly set'
                      END SELECT
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO

            ! check as is: strange data holes, all negative speeds are affected
            ! compare with wrfpress orig run, looks good, distribution OK
            ! compare with heimo old, identical
            ! insert SKn routine, simplified interpolation: looks OK
            ! could also compare with p_interp: processeding and code: not done
            ! comparison with wrfpress, same time: apart from the holes where
            ! negative values are with the more sophisticated scheme, all three
            ! datasets about identical: pattern and range, wrfpress, old scheme,
            ! new scheme
            ELSE IF ( (var_cmip(ivar) == "ua1000") .OR. &
                 (var_cmip(ivar) == "ua925") .OR.       &
                 (var_cmip(ivar) == "ua850") .OR.       &
                 (var_cmip(ivar) == "ua700") .OR.       &
                 (var_cmip(ivar) == "ua500") .OR.       &
                 (var_cmip(ivar) == "ua200") ) THEN
              ! rotate to earth grid and destagger
              PRINT *, "calc ua..."
              DO i = 1,xfocus
                 var3d_in(i,:,:) = ((u_in(i,:,:)+u_in(i+1,:,:))/2.)*cosalpha_in(:,:) - &
                                   ((v_in(i,:,:)+v_in(i+1,:,:))/2.)*sinalpha_in(:,:) 
!                var3d_in(i,:,:) = ((u_in(i,:,:)+u_in(i+1,:,:))/2.)
              END DO
              ! linear in p
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN
                      ! HTr: calculate zg at pout, first (linear in log(p)
                      !slope = (ph_fl(i,j,nl)-ph_fl(i,j,nl+1))/ &
                      !        (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                      !zg_pout = ph_fl(i,j,nl+1) + slope* &
                      !          (LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      ! HTr: interpolate logarithmic in zg
                      !slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ &
                      !        (ph_fl(i,j,nl)-ph_fl(i,j,nl+1))
                      !var_pl(i,j) = var3d_in(i,j,nl) * &
                      !                 EXP( (zg_pout - ph_fl(i,j,nl)) * &
                      !                      (LOG(var3d_in(i,j,nl+1)) - &
                      !                       LOG(var3d_in(i,j,nl))) / &
                      !                 (ph_fl(i,j,nl+1)-ph_fl(i,j,nl)))
                      slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                              (p_in(i,j,nl)-p_in(i,j,nl+1))
                      data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO

            ! comparison with wrfpress, same time: simple scheme about the same
            ! results
            ELSE IF ( (var_cmip(ivar) == "va1000") .OR. &
                 (var_cmip(ivar) == "va925") .OR.       &
                 (var_cmip(ivar) == "va850") .OR.       &
                 (var_cmip(ivar) == "va700") .OR.       &
                 (var_cmip(ivar) == "va500") .OR.       &
                 (var_cmip(ivar) == "va200") ) THEN
              ! rotate to earth grid and destagger
               DO j = 1,yfocus
                 var3d_in(:,j,:) = (v_in(:,j,:)+v_in(:,j+1,:))/2.*cosalpha_in(:,:) + &
                                   (u_in(:,j,:)+u_in(:,j+1,:))/2.*sinalpha_in(:,:) 
!                var3d_in(i,:,:) = ((v_in(i,:,:)+v_in(i+1,:,:))/2.)
               END DO
              ! linear in p
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN                
                      ! HTr: calculate zg at pout, first (linear in log(p)
                      !slope = (ph_fl(i,j,nl)-ph_fl(i,j,nl+1))/ (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                      !zg_pout = ph_fl(i,j,nl+1) + slope* (LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      ! HTr: interpolate logarithmic in zg
                      !slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ (ph_fl(i,j,nl)-ph_fl(i,j,nl+1))
                      !var_pl(i,j) = var3d_in(i,j,nl) * EXP( (zg_pout - ph_fl(i,j,nl)) * (LOG(var3d_in(i,j,nl+1)) - LOG(var3d_in(i,j,nl))) / (ph_fl(i,j,nl+1)-ph_fl(i,j,nl)))
                      slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                              (p_in(i,j,nl)-p_in(i,j,nl+1))
                      data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO

            ! comparison with orig model data between level 16 and 17
            ! OK
            ELSE IF ( (var_cmip(ivar) == "wa1000") .OR. &
                 (var_cmip(ivar) == "wa925") .OR.       &
                 (var_cmip(ivar) == "wa850") .OR.       &
                 (var_cmip(ivar) == "wa700") .OR.       &
                 (var_cmip(ivar) == "wa500") .OR.       &
                 (var_cmip(ivar) == "wa200") ) THEN
              ! vertically destagger:
              DO nl = 1,nz
                var3d_in(:,:,nl) = ( w_in(:,:,nl) + w_in(:,:,nl+1) ) / 2.
              END DO
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).le.p_in(i,j,nl) .and. pout(np).gt.p_in(i,j,nl+1)) then                
                      slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                              (p_in(i,j,nl)-p_in(i,j,nl+1))
                      data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                    END IF
                  END DO
                END DO
              END DO              
              !$OMP END PARALLEL DO

            ! vars needed: int.: ph_fl (> in: ph_in, phb_in); p_in (> in: pp_in, pb_in)
            !              out: var_pl, var3d_in
            !              tmp: pout, slope
            ! comparison with wrfpress, same time: wrfpress, old and new scheme
            ! about identical
            ELSE IF ( (var_cmip(ivar) == "zg1000") .OR. &
                 (var_cmip(ivar) == "zg925") .OR.       &
                 (var_cmip(ivar) == "zg850") .OR.       &
                 (var_cmip(ivar) == "zg700") .OR.       &
                 (var_cmip(ivar) == "zg500") .OR.       &    
                 (var_cmip(ivar) == "zg200") ) THEN
              PRINT *, "calc zg..."
              var3d_in(:,:,:) = ph_fl(:,:,:)
              ! either linear in p or log(p)
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN               
                      !slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ &
                      !        (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                      !var_pl(i,j) = var3d_in(i,j,nl+1) + &
                      !              slope*(LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      SELECT CASE (vint_type)
                      CASE(0)
                      ! (HTr) - linear in p
                        slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                                (p_in(i,j,nl)-p_in(i,j,nl+1))
                        data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                      CASE(1)
                      ! (HTr) - linear in log(p)
                        slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                                (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
                        data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(LOG(pout(np))-LOG(p_in(i,j,nl+1)))
                      CASE DEFAULT
                        PRINT *, 'CAUTION: unknown setting for vint_type'
                        STOP 'vint_type not properly set'
                      END SELECT
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO
            END IF

            ! alternative: use one routine for all
            ! using np here is rubbish and wrong > originally trickeled up as well...
!            DO i = 1,xfocus 
!              DO j = 1,yfocus
!                DO nl = 1,nz - 1
!                  IF (pout(np).le.p_in(i,j,nl) .and. pout(np).gt.p_in(i,j,nl+1)) then                
!                    ! SKn: original version
!                    !slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ (p_in(i,j,nl)-p_in(i,j,nl+1))
!                    !var_pl(i,j,np) = var3d_in(i,j,nl+1) + slope*(pout(np)-p_in(i,j,nl+1))
!                    ! HTr: change interpolation to be linear in log(p)
!                    slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/ (LOG(p_in(i,j,nl))-LOG(p_in(i,j,nl+1)))
!                    var_pl(i,j,np) = var3d_in(i,j,nl+1) + slope*(LOG(pout(np))-LOG(p_in(i,j,nl+1)))
!                  END IF
!                END DO
!              END DO
!            END DO

            !!data_in(:,:) = var_pl(:,:) ! superfluous, remove var_pl TODO was some old idea of someone but does not fit into the overall concept

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! psl [Pa] i Sea Level Pressure

            ! several options to calculate slp, either from formula
            ! internal vars needed: ph_fl, t_in, p_in
            IF ( var_cmip(ivar) == "psl" ) THEN

              SELECT CASE (calc_slp_type)
              CASE(0)
                ! HTr this calculation of mean sea level pressure gives quite high
                ! anomalies in mountainous areas
                PRINT *, "psl calc option 0: standard"
                psl_in(:,:) = (pp_in(:,:,1)+pb_in(:,:,1))*((t_in(:,:,1)*(1.+0.61*qv_in(:,:,1))+ &
                  0.0065*ph_fl(:,:,1))/(t_in(:,:,1)*(1+0.61*qv_in(:,:,1))))**(gr/(R*0.0065))
              CASE(1)
                ! modified version from wrf_interp.F90
                PRINT *, "psl calc option 1: p_interp"
                CALL calcslp(psl_in,p_in,qv_in,theta_in,ph_fl,nz,yfocus,xfocus,T00(1))
              CASE(2)
                ! officially agreed upon variant using ECMWF methodology as
                ! decided during autumn 2018 Trieste CORDEX FPS CEM meeting
                ! fullpos_cy38 implementation (ERA/IFS/APACHE/ALADIN)
                PRINT *, "psl calc option 2: fullpos_cy38 method"
                ! the lowest layer of ph+phb (nz+1 layers!) is exactly the surface)
                !            (slp,    PP,   P_s,     PHI_s,                      T_L,         nz, ns,     ew    )
                !slp        !(output) sea level pressure
                !PP         !(input) 3D pressure
                !P_s        !(input) pressure at surface
                !PHI_s      !(input) 2D geopotential of the surface
                !T_L        !(input) temperature at lowest level
                !nz, ns, ew !(input) dimensions: vertical, north-south, east-west
                CALL calcslptwo(psl_in, p_in, psfc_in, ph_in(:,:,1)+phb_in(:,:,1), t_in(:,:,1), nz, yfocus, xfocus)
                !                     int c    read   read          read         int c
              CASE DEFAULT
                PRINT *, 'CAUTION: unknown setting for calcslp, proceed with default method'
                STOP 'calc_slp_type not properly set'
              END SELECT

              data_in(:,:) = psl_in(:,:)

            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! prw [kg m-2] i Water Vapor Path

            IF ( (var_cmip(ivar) == "prw") ) THEN
  
              !t_in(:,:,:) = ( theta_in(:,:,:)+T00(1) ) * ( (pp_in(:,:,:)+pb_in(:,:,:))/P00(1) )**(R/cp)
              !p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)
  
              prw(:,:) = 0.
              DO nl = 1, nz
                prw(:,:) = prw(:,:) + &
                  (qv_in(:,:,nl) * p_in(:,:,nl)/(R*t_in(:,:,nl)) * &
                  ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+ &
                  phb_in(:,:,nl)))/gr)
              END DO
              data_in(:,:) = prw(:,:)
  
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! clwvi [kg m-2] i Condensed Water Path  

            IF ( (var_cmip(ivar) == "clwvi") ) THEN
    
              !t_in(:,:,:) = (theta_in(:,:,:)+T00(1))*((pp_in(:,:,:)+pb_in(:,:,:))/P00(1))**(R/cp)
              !p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)

              clwvi(:,:) = 0.
              DO nl = 1,nz - 1
                clwvi(:,:) = clwvi(:,:) + (qc_in(:,:,nl) + qi_in(:,:,nl) + &
                             qr_in(:,:,nl) + qs_in(:,:,nl) ) * p_in(:,:,nl)/ &
                             (R*t_in(:,:,nl)) * ((ph_in(:,:,nl+1)+ &
                             phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+ &
                             phb_in(:,:,nl)))/gr
              END DO
              data_in(:,:) = clwvi(:,:)

            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! clivi  [kg m-2] i Ice Water Path

            IF ( (var_cmip(ivar) == "clivi")) THEN
    
              !t_in(:,:,:) = (theta_in(:,:,:)+T00(1))*((pp_in(:,:,:)+pb_in(:,:,:))/P00(1))**(R/cp)
              !p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)
    
              clivi(:,:) = 0.
              DO nl = 1,nz - 1
                clivi(:,:) = clivi(:,:) + (qi_in(:,:,nl) + qs_in(:,:,nl)) * &
                             p_in(:,:,nl)/(R*t_in(:,:,nl)) * &
                             ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - &
                             (ph_in(:,:,nl)+phb_in(:,:,nl)))/gr
              END DO
              data_in(:,:) = clivi(:,:)
    
            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [ ], FPS CEM [X]
! cape [J kg-1] i 2-D Maximum Convective Available Potential Energy
! cin [J kg-1] i 2-D Maximum Convective Inhibition

            IF ( (var_cmip(ivar) == "cape") .OR. (var_cmip(ivar) == "cin") ) THEN
    
              !t_in(:,:,:) = (theta_in(:,:,:)+T00(1))*((pp_in(:,:,:)+pb_in(:,:,:))/P00(1))**(R/cp)
              !p_in = pp_in+pb_in          
    
              t_p(:,:,1) = t_in(:,:,1)
    
              cape(:,:) = 0.
              cin(:,:) = 0.
              lcl(:,:) = -999.
              lfc(:,:) = -999.
    
              !$OMP PARALLEL DO
              DO i = 1,xfocus
                DO j = 1,yfocus
                  DO nl = 1,nz-1
    
                    qvs(i,j,nl) = 0.622*a*exp(b*(t_p(i,j,nl)-c)/(t_p(i,j,nl)-d))/p_in(i,j,nl)
    
                    IF (qvs(i,j,nl) .gt. qv_in(i,j,1)) THEN ! dry adiabatic ascent
                   
                      t_p(i,j,nl+1) = (theta_in(i,j,1)+T00(1))*(p_in(i,j,nl+1)/P00(1))**(R/cp)
    
                    ELSE IF (qvs(i,j,nl) .lt. qv_in(i,j,1)) THEN ! moist adiabatic ascent
    
                      IF (lcl(i,j) .eq. -999) THEN ! lifting condensation level
                        lcl(i,j) = p_in(i,j,nl)
                      END IF
    
                      t_ii = t_p(i,j,nl)
    
                      DO ii = 1,10 ! solve iteratively
    
                        qvs(i,j,nl+1) = 0.622*a*exp(b*(t_ii-c)/(t_ii-d))/p_in(i,j,nl+1)
    
                        t_ii = t_ii - (t_ii*(P00(1)/p_in(i,j,nl+1))**(R/cp)*exp(L*qvs(i,j,nl+1)/(cp*t_ii)) - &
                               (t_p(i,j,nl)*(P00(1)/p_in(i,j,nl))**(R/cp)*exp(L*qvs(i,j,nl)/(cp*t_p(i,j,nl))))) / &
                               ( (P00(1)/p_in(i,j,nl+1))**(R/cp)*exp(n/(p_in(i,j,nl+1)*t_ii)*exp(b*(t_ii-c)/(t_ii-d))) * &
                               (1 - (n/p_in(i,j,nl+1)*exp(b*(t_ii-c)/(t_ii-d))*(t_ii*(t_ii-b*c)+(b-2)* &
                                 d*t_ii+d**2))/(t_ii*(d-t_ii)**2)) )
  
                      END DO
    
                      !print*, 'thetae(i,j,nl)',(t_p(i,j,nl)*(100000./p_in(i,j,nl))**(R/cp)*exp(L*qvs(i,j,nl)/(cp*t_p(i,j,nl))))
                      !print*,'thetae(i.j.nl+1)',(t_ii*(100000./p_in(i,j,nl+1))**(R/cp)*exp(L*qvs(i,j,nl+1)/(cp*t_ii))) 
  
                      t_p(i,j,nl+1) = t_ii
    
                      !print*, nl, 'moist', t_p(i,j,nl+1), t_in(i,j,nl+1), (t_p(i,j,nl+1)-t_in(i,j,nl+1))
    
                    END IF                 
  
                    IF (t_p(i,j,nl) .gt. t_in(i,j,nl)) THEN
                   
                      IF (lfc(i,j) .eq. -999) THEN ! level of free convection
                        lfc(i,j) = p_in(i,j,nl)
                      END IF
    
                      cape(i,j) = cape(i,j) + (t_p(i,j,nl) - t_in(i,j,nl)) / t_in(i,j,nl) * &
                        ((phb_in(i,j,nl)+ph_in(i,j,nl))-(phb_in(i,j,nl-1)+ph_in(i,j,nl-1)))
    
                      !print*, 'nl, cape(i,j)', nl, cape(i,j)
    
                    ELSE IF ( (t_p(i,j,nl) .lt. t_in(i,j,nl)) .and. (cape(i,j) .eq. 0.) ) THEN ! convective inhibition 
                 
                      cin(i,j) = cin(i,j) + (t_in(i,j,nl) - t_p(i,j,nl)) / t_in(i,j,nl) * &
                        ((phb_in(i,j,nl)+ph_in(i,j,nl))-(phb_in(i,j,nl-1)+ph_in(i,j,nl-1))) 
    
                      !print*, 'nl, cin(i,j)', nl, cin(i,j)
    
                    END IF
    
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO
  
              IF ( var_cmip(ivar) == "cape") THEN
                data_in(:,:) = cape(:,:)
              END IF
              IF ( var_cmip(ivar) == "cin") THEN
                data_in(:,:) = cin(:,:)
              END IF
    
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

          END IF ! pressure levels processing

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! clt [%] a Total Cloud Fraction
! not instantaneous, but average, have to recode here
  
          IF (var_cmip(ivar) == "clt") THEN

            IF (calc) THEN
  
              IF (.not. ALLOCATED(cldfra_inv)) ALLOCATE( cldfra_inv( xfocus, yfocus ), STAT=sts )

              IF (ALLOCATED(tmp_3d)) DEALLOCATE( tmp_3d )
              ALLOCATE( tmp_3d( xfocus, yfocus, 2 ), STAT=sts )

              DO k =1,2,1

                cldfra_inv(:,:) = 1.
  
                !$OMP PARALLEL DO
                DO i = 1,xfocus
                  DO j = 1,yfocus
                    IF (maxval(cldfra_in(i,j,:,k)) .lt. 0.99) THEN
                      cldfra_inv(i,j) = 1.
                      DO nl = 2,nz
                        cldfra_inv(i,j) = cldfra_inv(i,j) * &
                          (1- max(cldfra_in(i,j,nl,k),cldfra_in(i,j,nl-1,k)) / &
                          (1-cldfra_in(i,j,nl-1,k))) 
                      END DO
                    ELSE 
                      cldfra_inv(i,j) = 0.  
                    END IF
                  END DO ! j
                END DO ! i
                !$OMP END PARALLEL DO
  
                tmp_3d(:,:,k) = (1 - cldfra_inv(:,:))*100.
  
              END DO ! k

              WHERE (tmp_3d > 100.) tmp_3d = 100.
              WHERE (tmp_3d < 0.) tmp_3d = 0.

              data_in(:,:) = ( tmp_3d(:,:,1) + tmp_3d(:,:,2) ) / 2.

              DEALLOCATE(tmp_3d)

            ! this is the last field, i.e. no subsequent field can be found
            ! any more in the staged data
            ELSE

                data_in(:,:) = mv

            END IF

            calc = .TRUE.
  
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CORDEX [X], FPS CEM [X]
! pr [kg m-2 s-1] a Precipitation  
! [mm dtHours-1] -> [kg m-2 s-1]
! accumulated quantity (using a bucket)

          IF (var_cmip(ivar) == "pr") THEN

            IF (calc) THEN
  
              IF ( bucket_mm > 0. ) THEN

                !data_in(:,:) = ((rainnc_in(:,:,2) + rainc_in(:,:,2)) - (rainnc_in(:,:,1) + rainc_in(:,:,1)) + &
                !                (i_rainnc_in(:,:,2) + i_rainc_in(:,:,2) - i_rainnc_in(:,:,1) - i_rainc_in(:,:,1))*bucket_mm)/(dtHours*3600.) 

                data_in(:,:) = (((rainnc_in(:,:,2)+i_rainnc_in(:,:,2)*bucket_mm) + &
                                 (rainc_in(:,:,2)+i_rainc_in(:,:,2)*bucket_mm  )) - &
                                ((rainnc_in(:,:,1)+i_rainnc_in(:,:,1)*bucket_mm) + &
                                 (rainc_in(:,:,1)+i_rainc_in(:,:,1)*bucket_mm  ))) / &
                               (dtHours*3600.) 

              ELSE

                data_in(:,:) = ( ( rainnc_in(:,:,2) + rainc_in(:,:,2) ) - &
                                 ( rainnc_in(:,:,1) + rainc_in(:,:,1) ) ) / &
                                 ( dtHours * 3600. )

              END IF

            ! this is the last field, i.e. no subsequent field can be found
            ! any more in the staged data
            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prc [kg m-2 s-1] a Convective Precipitation 
! unit [mm/3hr] to [kg m-2 s-1]
  
          IF (var_cmip(ivar) == "prc") THEN

            IF (calc) THEN

              IF ( bucket_mm > 0. ) THEN
              
                data_in(:,:) = ( ( rainc_in(:,:,2)+i_rainc_in(:,:,2)*bucket_mm ) - &
                                 ( rainc_in(:,:,1)+i_rainc_in(:,:,1)*bucket_mm ) ) / &
                                 (dtHours*3600.) 

              ELSE

                data_in(:,:) = ( rainc_in(:,:,2) - rainc_in(:,:,1) ) / &
                               (dtHours*3600.)
 
              END IF

            ! this is the last field, i.e. no subsequent field can be found
            ! any more in the staged data
            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prsn [kg m-2 s-1] ? Snowfall Flux
! unit [mm/3hr] to [kg m-2 s-1]
  
          IF (var_cmip(ivar) == "prsn") THEN
  
            IF (calc) THEN
            
              data_in(:,:) = (snownc_in(:,:,2) - snownc_in(:,:,1)) / &
                             (dtHours*3600.)
  
            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! snm [kg m-2 s-1] a Surface Snow Melt
! unit [kg m-2 /3hr] to [kg m-2 s-1]
  
          IF (var_cmip(ivar) == "snm") THEN
  
            IF (calc) THEN
              data_in(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1)) / (dtHours*3600.)
              WHERE (landmask_in == 0) data_in = mv
            ELSE
              data_in(:,:) = mv
            END IF
            calc = .TRUE.

          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! evspsbl [kg m-2 s-1] a Evaporation
! unit [kg m-2 /3hr] to [kg m-2 s-1]
! TO CHECK STILL, in principle OK
  
          IF (var_cmip(ivar) == "evspsbl") THEN
  
            IF (calc) THEN
            
              data_in(:,:) = ( sfcevp_in(:,:,2) - sfcevp_in(:,:,1) ) / &
                             (dtHours*3600.)
  
            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! evspblpot [kg m-2 s-1] a Potential Evapotranspiration 
! unit [W m-2]/[J kg-1] -> [kg m-2 s-1]
! TO CHECK STILL, in principle OK
! TODO / CHECK
! THERE IS STH WRONG WITH THE UNITS: WRF's POTEVP is accumulated and declared to be in W m-2. 
! It doesen't make sense to accumulate in W m-2, but even if assume it as J m-2 or derive kg m-2 
! by using latent heat of vaporization you never get values that have a reasonable magnitude...
! CHECK with registry on the units
  
          IF (var_cmip(ivar) == "evspsblpot") THEN
  
            IF (calc) THEN
            
              data_in(:,:) = (potevp_in(:,:,2) - (potevp_in(:,:,1))) / L
  
            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrros [kg m-2 s-1] a Surface Runoff
! unit [mm/3hr] to [kg m-2 s-1]
  
          IF (var_cmip(ivar) == "mrros") THEN
  
            IF (calc) THEN
              data_in(:,:) = (sfroff_in(:,:,2) - sfroff_in(:,:,1)) / &
                             (dtHours*3600.)
              WHERE (landmask_in == 0) data_in = mv
            ELSE
              data_in(:,:) = mv
            END IF
            calc = .TRUE.

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrro [kg m-2 s-1] a Total Runoff
! unit [mm/3hr] to [kg m-2 s-1]
! SFROFF, UDROFF
! UDROFF is accumulated!
  
          IF (var_cmip(ivar) == "mrro") THEN
  
            IF (calc) THEN
              data_in(:,:) = ( (sfroff_in(:,:,2) - sfroff_in(:,:,1)) + &
                               (udroff_in(:,:,2) - udroff_in(:,:,1)) ) / &
                               (dtHours*3600.) 
              WHERE (landmask_in == 0) data_in = mv
            ELSE
              data_in(:,:) = mv
            END IF
            calc = .TRUE.

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! hfss [W m-2] a Surface Upward Latent Heat Flux
! hfls [W m-2] a Surface Upward Sensible Heat Flux
  
          IF ( (var_cmip(ivar) == "hfss") &
            .OR. (var_cmip(ivar) == "hfls") ) THEN
    
            IF (calc) THEN

              data_in(:,:) = (hflux_in(:,:,1) + hflux_in(:,:,2)) / 2.

            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! rsds [W m-2] a Surface Downwelling Shortwave Radiation
! rlds [W m-2] a Surface Downwelling Longwave Radiation
! rsus [W m-2] a Surface Upwelling Shortwave Radiation 
! rlus [W m-2] a Surface Upwelling Longwave Radiation
! rlut [W m-2] a TOA Outgoing Longwave Radiation
! rsdt [W m-2] a TOA Incident Shortwave Radiation
! rsut [W m-2] a TOA Outgoing Shortwave Radiation
  
          IF ( (var_cmip(ivar) == "rsds") &
            .OR. (var_cmip(ivar) == "rlds") &
            .OR. (var_cmip(ivar) == "rsus") &
            .OR. (var_cmip(ivar) == "rlus") &
            .OR. (var_cmip(ivar) == "rlut") &
            .OR. (var_cmip(ivar) == "rsdt") &
            .OR. (var_cmip(ivar) == "rsut") ) THEN
    
            IF (calc) THEN

              IF ( bucket_J > 0. ) THEN

                data_in(:,:) = ( ( rad_in(:,:,2) + ( i_rad_in(:,:,2) * bucket_J ) ) - &
                               ( rad_in(:,:,1) + ( i_rad_in(:,:,1) * bucket_J ) ) ) / &
                               (dtHours*3600.)
               
              ELSE

                data_in(:,:) = (rad_in(:,:,2) - rad_in(:,:,1)) / (dtHours*3600.)

              END IF

            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sund [s] s Duration of Sunshine
! WMO definition >= 120 W m-2

          IF ( var_cmip(ivar) == "sund" ) THEN
    
            IF (calc) THEN

              DO i = 1,xfocus
                DO j = 1,yfocus

                  ! if both neighbouring data points are above the threshold
                  IF ( (swdown_in(i,j,1) >= 120.) .AND. (swdown_in(i,j,2) >= 120.) ) THEN 
                    data_in(i,j) = dtHours*3600.
                  ELSE IF ( (swdown_in(i,j,1) < 120.) .AND. (swdown_in(i,j,2) < 120.) ) THEN 
                    data_in(i,j) = 0.
                  ELSE
                    slope = ( swdown_in(i,j,2) - swdown_in(i,j,1) ) / &
                            ( dtHours * 3600. ) 
                    IF ( slope > 0 ) THEN 
                      data_in(i,j) = ( dtHours * 3600. ) - ( (120. - swdown_in(i,j,1)) / slope )
                    ELSE IF ( slope < 0 ) THEN
                      data_in(i,j) = (120. - swdown_in(i,j,1)) / slope
                    END IF
                  END IF
               
                END DO
              END DO

            ELSE

              data_in(:,:) = mv

            END IF

            calc = .TRUE.

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrso [kg m-2] i Total Soil Moisture Content 
! see comment in the variable aquisition section
! m3 m-3 -> kg m-2
! double checked in July 2019 
! SMOIS_i [m3 m-3] * DZS_i [m] * 1000 [kg m-3] = [kg m-2]
! assume standard area 1m2, assume 1l water = 1kg, 1mm m-2 water depth = 1l
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrso") THEN
    
            !data_in(:,:) = ((smois_in(:,:,1,1)*DZS(1) + smois_in(:,:,2,1)*DZS(2) + smois_in(:,:,3,1)*DZS(3) + smois_in(:,:,4,1)*DZS(4) ) + &
            !                (smois_in(:,:,1,2)*DZS(1) + smois_in(:,:,2,2)*DZS(2) + smois_in(:,:,3,2)*DZS(3) + smois_in(:,:,4,2)*DZS(4) ))/2.*1000. 
  
            data_in(:,:) =  (smois_in(:,:,1,1)*DZS(1) + &
                             smois_in(:,:,2,1)*DZS(2) + &
                             smois_in(:,:,3,1)*DZS(3) + &
                             smois_in(:,:,4,1)*DZS(4)) * 1000.
            WHERE (landmask_in == 0) data_in = mv

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrsol
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrsol") THEN

            DO i = 1, 4
              WHERE (landmask_in(:,:) == 1)
                data_in_3D(:,:,i) = smois_in(:,:,i,1) * DZS(i) * 1000.
              ELSEWHERE
                data_in_3D(:,:,i) = mv
              END WHERE
            END DO

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! tsl
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "tsl") THEN

            DO i = 1, 4
              WHERE (landmask_in(:,:) == 1)
                data_in_3D(:,:,i) = smois_in(:,:,i,1)
              ELSEWHERE
                data_in_3D(:,:,i) = mv
              END WHERE
            END DO

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! snc [%] i Snow Area Fraction
! unit [] to [%]
! masking over ocean, with mv initialized

          IF ( (var_cmip(ivar) == "snc") ) THEN

            !PRINT *, 'shape data_in with snc', SHAPE(data_in)
            !PRINT *, 'landmask_in test snc middle', data_in(180, :)*100. 
            !PRINT *, 'landmask_in', SHAPE(landmask_in)
            !PRINT *, 'landmask_in middle', landmask_in(180, :)
            !data_in = landmask_in

            data_in(:,:) = data_in(:,:)*100. 
            WHERE (landmask_in == 0) data_in = mv

            !WHERE (landmask_in == 1)
            !  data_in = data_in*100. 
            !ELSEWHERE
            !  data_in = mv
            !END WHERE

          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! snd [m] i Snow Depth
! snw [kg m-2] i Surface Snow Amount
! masking over ocean, with mv initialized
! variables are just passed through and need some masking

          IF ( (var_cmip(ivar) == "snd") &
            .OR. (var_cmip(ivar) == "snw") ) THEN
 
            WHERE (landmask_in == 0) data_in = mv
            !WHERE (landmask_in(:,:) == 1)
            !  data_in(:,:) = data_in(:,:)
            !ELSEWHERE
            !  data_in(:,:) = mv
            !END WHERE
  
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sic [%] ?>i Sea Ice Area Fraction
! no temporal aggregation defined, treat as tier-2 instantaneous data, most
! reasonable; also offers the possibility of some proper sea-ice treatment

          IF ( (var_cmip(ivar) == "sic") ) THEN
 
            ! fractional sea ice 
!           IF (fractSeaIce) THEN

              data_in(:,:) = data_in(:,:) * 100. 

            ! binary sea ice mask, SEAICE is empty variable
!           ELSE
!
!             WHERE ( xland_in == 2. )
!
!               xland_in = 0.
!
!             END WHERE
!
!              data_in(:,:) = ( landmask_in(:,:) - xland_in(:,:) ) * 100.
!
!            END IF
  
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! -1] i Near-Surface Wind Speed
  
          IF (var_cmip(ivar) == "sfcWind") THEN
  
            data_in(:,:) = (u10_in(:,:)**2 + v10_in(:,:)**2)**0.5 
  
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! uas [m s-1] i Eastward Near-Surface Wind
  
          IF (var_cmip(ivar) == "uas") THEN 
  
            data_in(:,:) = u10_in(:,:)*cosalpha_in(:,:) - v10_in(:,:)*sinalpha_in(:,:) ! rotate to earth grid
  
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! vas [m s-1] i Northward Near-Surface Wind
  
          IF (var_cmip(ivar) == "vas")  THEN
  
            data_in(:,:) = v10_in(:,:)*cosalpha_in(:,:) + u10_in(:,:)*sinalpha_in(:,:) ! rotate to earth grid
  
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prhmax [kg m-2 s-1] m Daily Maximum Hourly Precipitation Rate (using wrfxtrm)

          IF (var_cmip(ivar) == "prhmax") THEN

            data_in(:,:) = rainnc_max_in(:,:) + rainc_max_in(:,:)

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! hurs [%] i Near-Surface Relative Humidity
! for compatibility, copy the approach as implemented by IDL (Soares, Cardoso, Careto): RCM_convection_sfc_hum.f90 line 305ff 
! Following WMO, supercooled water is assumed for temperatures below 0ºC and esat is always calculated in reference to water.
! Enhancement factor is used calculate the effective esat in the presence of other gases
! REFERENCES
! Wexler, A., Vapor Pressure Formulation for Water in Range 0 to 100°C. A Revision, Journal of Research of the National Bureau of Standards  A. Physics and Chemistry, September  December 1976, Vol. 80A, Nos.5 and 6, 775-785
! Wexler, A., Vapor Pressure Formulation for Ice, Journal of Research of the National Bureau of Standards  A. Physics and Chemistry, January  February 1977, Vol. 81A, No. 1, 5-19
! Goff, J. A., Standardization of Thermodynamic Properties of Moist Air, Heating, Piping, and Air Conditioning, 1949, Vol. 21, 118.

          IF (var_cmip(ivar) == "hurs") THEN

            ! t2=t2_in 
            ! mr=q2_in 
            ! psf=psfc_in
            !        TK                    mr_sat                 e_sfc                  esat                   fact
            ALLOCATE(tmp_2d(xfocus,yfocus),tmp1_2d(xfocus,yfocus),tmp2_2d(xfocus,yfocus), &
              tmp3_2d(xfocus,yfocus),tmp4_2d(xfocus,yfocus))

            tmp_2d = t2_in-c
            tmp2_2d = (q2_in / (epsil+q2_in)) * psfc_in
            tmp3_2d = exp( a1*t2_in**(-2) + a2/t2_in + a3 + a4*t2_in + a5*t2_in**2 + a6*t2_in**3 + a7*t2_in**4 + a8*log(t2_in) )

            WHERE (tmp_2d < 0.)
               tmp4_2d = exp((fi1+fi2*tmp_2d+fi3*tmp_2d**2+fi4*tmp_2d**3)*(1.-(tmp3_2d/psfc_in))+ &  
                          exp(fi5+fi6*tmp_2d+fi7*tmp_2d**2+fi8*tmp_2d**3)*((tmp3_2d/psfc_in)-1.d0))
            ELSEWHERE
               tmp4_2d = exp((fw1+fw2*tmp_2d+fw3*tmp_2d**2+fw4*tmp_2d**3)*(1.-(tmp3_2d/psfc_in))+ &  
                          exp(fw5+fw6*tmp_2d+fw7*tmp_2d**2+fw8*tmp_2d**3)*((tmp3_2d/psfc_in)-1.d0))
            END WHERE

            tmp1_2d = epsil*(tmp3_2d*tmp4_2d) / (psfc_in-(tmp3_2d*tmp4_2d))
            data_in = max(min(q2_in/tmp1_2d,1.),0.0)*100.

            DEALLOCATE(tmp_2d,tmp1_2d,tmp2_2d,tmp3_2d,tmp4_2d)

            ! too simplified
            !data_in(:,:) = q2_in(:,:) / ( (379.90516 / psfc_in(:,:)) * exp(17.2693882*(t2_in(:,:)-273.15) / (t2_in(:,:) - 35.86)) )
            !data_in(:,:) = 0.01 * (psfc_in(:,:)*q2_in(:,:)/(q2_in(:,:)*(1.-0.622) + 0.622)) / (611.2*exp(17.67*(t2_in(:,:)-273.15)/(t2_in(:,:)-29.65)))

          END IF

!-------------------------------------------------------------------------------
! write data to netCDF file

          PRINT *, "*** WRITE DATA TO netCDF ***"
          PRINT *, TRIM(pn_out) // "/" // TRIM(fn_out)
          PRINT *, TRIM(var_cmip(ivar)), xfocus, yfocus, counter, ncid, x_varid
  
          sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) &
            // "/" // TRIM(fn_out), NF90_WRITE, ncid )
          PRINT *, "NF90_OPEN",  sts

          ! file must exist, just from the logic of the code, nevertheless: test
          IF (sts/=0) THEN
            PRINT *, "NF90_OPEN",  sts
            PRINT *, TRIM(pn_out) // "/" // TRIM(fn_out) // "  FAILED - EXIT"
            EXIT
          END IF
          
          sts = NF90_INQ_VARID(ncid, TRIM(var_cmip(ivar)), x_varid)
          PRINT *, "NF90_OPEN",  sts
  
            PRINT *, 'NF90_INQ_VARID', ncid
            PRINT *, 'var_cmip(ivar)', var_cmip(ivar)
            PRINT *, 'x_varid', x_varid
            PRINT *, 'counter/offset', counter
            PRINT *, 'xfocus', xfocus
            PRINT *, 'yfocus', yfocus
    
            ! most CMOR vars are not 3D, pressure level data is one field and
            ! level per file; but there are exceptions, where data is stored 3D
            ! including the repective z-coordinates, all according to CF
            ! convention
            IF ((var_cmip(ivar) == "mrsol") &
              .OR. (var_cmip(ivar) == "tsl")) THEN
              sts = NF90_PUT_VAR( ncid, x_varid, data_in_3D(:,:,:), &
                START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
            !IF ( height(ivar) /= -999 ) THEN
            !  sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
            !    START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1, 1 /) )
            ELSE
              sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
                START=(/ 1, 1, counter /), COUNT = (/ xfocus, yfocus, 1 /) )
            END IF
    
            PRINT *, 'NF90_PUT_VAR', sts
            PRINT *, 'ncid', ncid
            PRINT *, 'x_varid', x_varid
            PRINT *, 'some sample output 3x3 in the middle of the domain', &
              data_in(xfocus/2:(xfocus/2+2),yfocus/2:(yfocus/2+2))

          sts = NF90_CLOSE(ncid)
          PRINT *, "NF90_CLOSE",  sts

!-------------------------------------------------------------------------------

          ELSE

            ! tested this: just define file with a temporal coverage outside the
            ! range of the input file: create file, loop over the inputs,
            ! nothing happens, no false sorting in 
            PRINT *, "no time match", time_match

          END IF

!-------------------------------------------------------------------------------
  
        END DO ! it i-time WRF indiv file loop
  
!-------------------------------------------------------------------------------
! next WRF file may contain different number of output intervals

        DEALLOCATE(InVarDataRec)

        DEALLOCATE(InDateTimeYear)
        DEALLOCATE(InDateTimeMonth)
        DEALLOCATE(InDateTimeDay)
        DEALLOCATE(InDateTimeHour)
        DEALLOCATE(InDateTimeMinute)
        DEALLOCATE(InDateTimeSecond)

        DEALLOCATE(InDateTimeCombined)

!-------------------------------------------------------------------------------

        PRINT *, "-----------------------------------------------------------"
        CALL CPU_TIME(cpuTe)
        PRINT '("CPU timing for the processing of one WRF file and one output variable = ",F10.1," s")',cpuTe-cpuTs

      END DO ! ifl - specific WRF input file, filelist loop

      ! must be reset, when moving to new var, always check file existance
      InDateTimeYearPrev = 0
      InDateTimeMonthPrev = 0
      prevpass = 0 ! ?????????????????????? one level inside????

      ENDIF ! procflag T/F TODO fix indention

    END DO ! ivar - variable loop

    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
    CALL mpi_finalize(ierr)
    IF ( ierr /= MPI_SUCCESS ) stop "Problem with MPI_FINALIZE" 

    DEALLOCATE( ivar_list ) 

  END DO ! ivarnml - namelist loop 

!END DO ! ifrq - different temporal aggregations

!===============================================================================

END PROGRAM WRFCMORizer

!===============================================================================
! HTr routines for psl calculation

SUBROUTINE calcslp(slp, pres, qv, tk1, ght, nz, ns, ew, T00)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nz,ns,ew
REAL, DIMENSION(:,:), INTENT(INOUT) :: slp
REAL, DIMENSION(:,:,:), INTENT(IN) :: pres,qv,tk1,ght
REAL, INTENT(IN) :: T00
INTEGER :: i,j,k,klo,khi
INTEGER, DIMENSION(ew,ns) :: level
REAL, DIMENSION(ew, ns) :: t_sea_level, t_surf
REAL, DIMENSION(ew, ns, nz) :: tk
REAL :: rgas,grav,gamma
REAL :: plo , phi , tlo, thi , zlo , zhi
REAL :: p_at_pconst , t_at_pconst , z_at_pconst
REAL :: z_half_lowest
REAL :: tc,pconst
LOGICAL :: l1, l2, l3, found, ridiculous_mm5_test

rgas = 287.04
grav = 9.81
gamma = 0.0065
tc = 273.16+17.5

! in wrf_interp pconst is set to 10000 Pa, but this gives too low pressure 
! values at sea level. So, I've reduced pconst to 5000 Pa, which reduces this 
! underestimation.
!pconst = 10000.
pconst = 5000.
ridiculous_mm5_test = .TRUE. 

! Find least zeta level that is PCONST Pa above the surface. We later use this
! level to extrapolate a surface pressure and temperature, which is supposed
! to reduce the effect of the diurnal heating cycle in the pressure field.
tk = tk1 + T00
!$OMP PARALLEL DO
DO j = 1, ns
  DO i = 1, ew

    level(i,j) = -1
    k = 1
    found = .false.

    DO WHILE( (.NOT. found) .AND. (k .LE. nz) )
      IF ( pres(i, j, k) .LT. pres(i, j, 1)-PCONST ) THEN
        level(i, j) = k
        found = .true.
      END IF
      k = k+1
    END DO 

    IF ( level(i,j) .EQ. -1 ) THEN
      PRINT '(A,I4,A)','Troubles finding level ',NINT(PCONST)/100,' above ground.'
      PRINT '(A,I4,A,I4,A)','Problems first occur at (',i,',',j,')'
      PRINT '(A,F6.1,A)','Surface pressure = ',pres(i,j,1)/100,' hPa.'
      STOP 'Error_in_finding_100_hPa_up'
    END IF

  END DO ! ew
END DO ! ns
!$OMP END PARALLEL DO
  
! Get temperature PCONST Pa above surface. Use this to extrapolate 
! the temperature at the surface and down to sea level.
!$OMP PARALLEL DO
DO j = 1, ns
  DO i = 1, ew

    klo = MAX ( level(i,j) - 1 , 1      )
    khi = MIN ( klo + 1        , nz - 1 )
     
    IF ( klo .EQ. khi ) THEN
      PRINT '(A)','Trapping levels are weird.'
      PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi,': and they should not be equal.'
      STOP 'Error_trapping_levels'
    END IF

    plo = pres(i,j,klo)
    phi = pres(i,j,khi)
    tlo = tk(i,j,klo) * (1. + 0.608 * qv(i,j,klo) )
    thi = tk(i,j,khi) * (1. + 0.608 * qv(i,j,khi) )
    zlo = ght(i,j,klo)         
    zhi = ght(i,j,khi)

    p_at_pconst = pres(i,j,1) - pconst
    t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
    z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

    t_surf(i,j) = t_at_pconst*(pres(i,j,1)/p_at_pconst)**(gamma*rgas/grav)
    t_sea_level(i,j) = t_at_pconst+gamma*z_at_pconst

  END DO ! ew
END DO !ns
!$OMP END PARALLEL DO

! If we follow a traditional computation, there is a correction to the sea level 
! temperature if both the surface and sea level temnperatures are *too* hot.
IF ( ridiculous_mm5_test ) THEN
  !$OMP PARALLEL DO
  DO j = 1, ns
    DO i = 1, ew
      l1 = t_sea_level(i,j) .LT. TC 
      l2 = t_surf     (i,j) .LE. TC
      l3 = .NOT. l1
      IF ( l2 .AND. l3 ) THEN
        t_sea_level(i,j) = TC
      ELSE
        t_sea_level(i,j) = TC - 0.005*(t_surf(i,j)-TC)**2
      END IF
    END DO
  END DO
  !$OMP END PARALLEL DO
END IF

! The grand finale: ta da!
!$OMP PARALLEL DO
DO j = 1, ns
  DO i = 1, ew
    z_half_lowest=ght(i,j,1)
    slp(i,j) = pres(i,j,1) *EXP((2.*grav*z_half_lowest)/(rgas*(t_sea_level(i,j)+t_surf(i,j))))
    !slp(i,j) = slp(i,j) * 0.01
  END DO ! ew
END DO ! ns
!$OMP END PARALLEL DO

END SUBROUTINE calcslp

!===============================================================================

!(lah): fullpos_cy38 routine
!(lah): this code resembles the algorithm described in chapter 4 of "FULL-POS in
!       the cycle 38 of ARPEGE/IFS" by Yessad K. (Meteo-France/CNRM/GMAP/ALGO),
!       November 10, 2011. The equation can be found especially in the section
!       (4.3.1) "Mean sea level pressure PP_msl (routine PPPMER)"

SUBROUTINE calcslptwo(slp, PP, P_s, PHI_s, T_L, nz, ns, ew)

IMPLICIT NONE

real, DIMENSION(:,:),   INTENT(OUT) :: slp        !(output) sea level pressure
real, DIMENSION(:,:,:), INTENT(IN)  :: PP         !(input) 3D pressure
real, DIMENSION(:,:),   INTENT(IN)  :: P_s        !(input) pressure at surface
real, DIMENSION(:,:),   INTENT(IN)  :: PHI_s      !(input) 2D geopotential of the surface
real, DIMENSION(:,:),   INTENT(IN)  :: T_L        !(input) temperature at lowest level
INTEGER,                INTENT(IN)  :: nz, ns, ew !(input) dimensions: vertical, north-south, east-west

real, DIMENSION(ew,ns)              :: P_L        !(calculated) pressure at lowest level
real                                :: T_surf     !(calculated) surface temperature
real                                :: gamma_mod  !(calculated) modified lapse rate that is actually used for calculations
real                                :: x          !(calculated) expansion coefficient of (eq. 9)
real                                :: T_0        !(calculated) auxiliary variable

real, PARAMETER                     :: Rd    = 287.04 ![J kg-1 K-1] (const.) dry air constant
real, PARAMETER                     :: g     = 9.81   ![m s-2] (const.) acceleration due to gravity
real, PARAMETER                     :: gamma = 0.0065 ![K m-1] (const.) lapse rate at const. 0.0065 K/m, also denoted as (dT/dz)_st
integer                             :: j,k !loop parameters

!print *,'this is calcslptwo !!!'

P_L = PP(:,:,1)  !extract lowest layer
  
!$OMP PARALLEL DO
DO j = 1 , ew
  DO k = 1 , ns
        
    !always assume none of the IF conditions trigger, then (according to step 5)
    gamma_mod = gamma
    
    !(0) if abs(PHI_s)<0.001 ("sea" grid cells) then set slp to surface pressure
    IF (ABS(PHI_s(j,k)).lt.0.001) THEN
      
      slp(j,k) = P_s(j,k) !in this case we are done with this grid cell
       
    ELSE !else apply the following algorithm

      !always assume none of the following IF conditions trigger
      !then we use the constant camma (according to step 5)
      gamma_mod = gamma
      
      !(1) compute T_surf according to eq (1)
      !T_surf = T_L + gamma*(Rd/g)*(P_s/P_L-1.0)*T_L
      T_surf = T_L(j,k) + gamma * (Rd/g) * ( P_s(j,k)/P_L(j,k) - 1.0 ) * T_L(j,k)
      
      !(2) compute T_0=T_surf+gamma*PHI_s/g
      T_0 = T_surf + gamma * PHI_s(j,k) / g
      
      !(3) to avoid extrapolation of too low pressures over high and warm surfaces:
      ! if (T_0 > 290.5):
      !    if (T_surf <= 290.5): gamma_mod=(290.5-T_surf)*g/PHI_s (eq. 7)
      !    else: gamma_mod=0.0, T_surf=0.5*(290.5+T_surf)
      IF (T_0 .gt. 290.5) THEN
        IF (T_surf .le. 290.5) THEN
          gamma_mod = ( 290.5 - T_surf ) * g / PHI_s(j,k)
        ELSE
          gamma_mod = 0.0
          T_surf = 0.5 * ( 290.5 + T_surf )
        END IF
      END IF
         
      !(4) to avoid extrapolation of too high pressures over cold surfaces:
      ! if T_surf < 255: gamma_mod=gamma, T_surf=0.5*(255+T_surf)
      IF (T_surf .lt. 255.0) THEN
        gamma_mod = gamma
        T_surf = 0.5 * ( 255.0 + T_surf )
      END IF
      
      !(5) in other cases set gamma_mod=gamma
      !this was already done in the beginning of the loop!
      
      !(6) compute mean sea level pressure (eq. 8) using the above determined parameters
      !x=gamma_mod*PHI_s/(g*T_surf)
      !slp=P_s*exp(PHI_s/(R_d*T_surf)*(1-x/2.+x*x/3.))
      x = gamma_mod * PHI_s(j,k) / ( g * T_surf )
      slp(j,k) = P_s(j,k) * EXP( PHI_s(j,k) / (Rd * T_surf ) * (1.0 - x/2. + x*x/3.) )

      !done.  
    
    END IF
 
  END DO
END DO
!$OMP END PARALLEL DO
  
END SUBROUTINE calcslptwo

!===============================================================================

SUBROUTINE GenerateFilelist

USE FilelistHandling

IMPLICIT NONE

CHARACTER (len = 200) :: ifl
INTEGER :: i, IOstatus, nfl, AllocateStatus

IF (ft == 0) THEN
  OPEN(2,FILE=tmpfileFL_std,STATUS='old')
END IF
IF (ft == 1) THEN
  OPEN(2,FILE=tmpfileFL_xtrm,STATUS='old')
END IF
IF (ft == 2) THEN
  OPEN(2,FILE=tmpfileFL_3d,STATUS='old')
END IF

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
IF (ft == 2) THEN
  ALLOCATE(fl_wrfpres(nfl), STAT=AllocateStatus)
END IF
!IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

REWIND(2)
DO i = 1,nfl
  IF (ft == 0) THEN
    READ(2,FMT='(a)') fl_wrfout(i)
  END IF
  IF (ft == 1) THEN
    READ(2,FMT='(a)') fl_wrfxtr(i)
  END IF
  IF (ft == 2) THEN
    READ(2,FMT='(a)') fl_wrfpres(i)
  END IF
END DO

CLOSE(2)

END SUBROUTINE GenerateFilelist

!===============================================================================
! TODO: expand for monthly and seasonal data

SUBROUTINE CreateRefTimeArray(dt)

USE RefTimeVecs
USE NamelistHandling

IMPLICIT NONE

CHARACTER (LEN = 3), INTENT(IN) :: dt ! 1hr 3hr 6hr day mon sem

INTEGER :: i, j, k, l, counter
REAL(KIND=8) :: dtDecDay
INTEGER :: tstotYYYY, tstotMM, tstotDD, tstotHH
INTEGER :: tetotYYYY, tetotMM, tetotDD, tetotHH
REAL(KIND=8) :: tstot_singlenumber, tetot_singlenumber
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp2D_singlenumber
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp2D
INTEGER, DIMENSION(12) :: ndpm ! number of days per month
INTEGER :: ndOverall = 31 ! number of days initialized with the 31 days of Dec 1949, this is DIRTY
INTEGER :: ntspd  ! (n)umber (t)ime(s)teps (p)er (d)ay
REAL :: cpuTs, cpuTe

CALL CPU_TIME(cpuTs)

PRINT *, "CreateRefTimeArray subroutine"
PRINT *, "Overall start end end: ", tstot, " to ", tetot
PRINT *, "temporal resolution working on = ", dt
PRINT *, "temporal aggregation working on = "

SELECT CASE (dt)
CASE ('1hr')
  dtDecDay = 1.0 / 24.0_8
CASE ('3hr')
  dtDecDay = 0.125
CASE ('6hr')
  dtDecDay = 0.25
CASE ('day') ! always means or sums, centered at 12UTC per day
  dtDecDay = 1.0
CASE DEFAULT
  PRINT *, "invalid time interval specified"
  STOP
END SELECT
ntspd = INT(1.0 / dtDecDay)
PRINT *, "number of time steps per day (ntspd) = ", ntspd

! "tstot" contains the absolute starting point, protocol = 1949-12-01_00:00:00
READ( tstot, '(I4,1X,I2,1X,I2,1X,I2)' ) tstotYYYY, tstotMM, tstotDD, tstotHH
READ( tetot, '(I4,1X,I2,1X,I2,1X,I2)' ) tetotYYYY, tetotMM, tetotDD, tetotHH
PRINT *, tstotYYYY, tstotMM, tstotDD, tstotHH
PRINT *, tetotYYYY, tetotMM, tetotDD, tetotHH

! e.g., 1949-2101, too many steps here, cut later
DO i=tstotYYYY,tetotYYYY,1
  ndOverall = ndOverall + CheckForLeapyear( i )
END DO

! index, YYYY, MM, DD, hh
ALLOCATE( tmp2D( ndOverall*ntspd, 5 ) )
PRINT *, "number of date/time steps in the raw time ref array", SIZE(tmp2D) / 5

! create timesteps [decimal days], starting with 0.0 
DO i=0_8,ndOverall*ntspd-1,1
  tmp2D( i+1, 1 ) = i * dtDecDay
END DO

! fill up Y M D h date/time information
counter = 0
DO i=tstotYYYY,tetotYYYY,1
  ! if there is no-leapyear dataset, switch this here off
  ! here, also a 360day calendar can be implemented, just set ndpm(:) = 30
  IF ( CheckForLeapyear( i ) == 366 ) THEN 
    ndpm = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  ELSE
    ndpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  END IF
  DO j=1,12,1
    ! sort in on daily basis
    DO k=1,ndpm(j)
      tmp2D( counter*ntspd+1 : counter*ntspd+ntspd , 2) = i
      tmp2D( counter*ntspd+1 : counter*ntspd+ntspd , 3) = j
      tmp2D( counter*ntspd+1 : counter*ntspd+ntspd , 4) = k
      tmp2D( counter*ntspd+1 : counter*ntspd+ntspd , 5) = (/(l, l=0, 24-24/ntspd , 24/ntspd) /) ! 0,18,6 -> 0, 6, 12, 18
      counter = counter + 1
    END DO
  END DO
END DO

!PRINT "(5F10.6)", TRANSPOSE(tmp2D(1:10,:))

! the ref time vec as is now is most likely too big, starting, e.g., 19490-01-01_00, ending 2101-12-31_23 with 1hr data
! here do some subsetting
tstot_singlenumber = REAL(tstotYYYY,8)*1000000._8 + REAL(tstotMM,8)*10000._8 + REAL(tstotDD,8)*100._8 + REAL(tstotHH,8)
tetot_singlenumber = REAL(tetotYYYY,8)*1000000._8 + REAL(tetotMM,8)*10000._8 + REAL(tetotDD,8)*100._8 + REAL(tetotHH,8)
PRINT *, tstot_singlenumber, tetot_singlenumber

ALLOCATE( tmp2D_singlenumber( ndOverall*ntspd ) )
DO i=1,ndOverall*ntspd,1
  tmp2D_singlenumber(i) = tmp2D(i,2)*1000000._8 + tmp2D(i,3)*10000._8 + tmp2D(i,4)*100._8 + tmp2D(i,5)
END DO

counter=1_4
l=0_4
DO i=1_4,ndOverall*ntspd,1
  IF ( tmp2D_singlenumber(i) == tstot_singlenumber ) THEN
    j=counter
  END IF
  IF ( tmp2D_singlenumber(i) == tetot_singlenumber ) THEN
    k=counter
  END IF
  IF ( ( tmp2D_singlenumber(i) >= tstot_singlenumber ) .AND. &
       ( tmp2D_singlenumber(i) <= tetot_singlenumber ) ) THEN
    l = l + 1_4
  END IF
  counter = counter + 1_4
END DO

PRINT *, "subset bounds lower upper, # elements ", j, k, l
! (HTr) - deallocate TimeRefArray in the case it has been allocated in a previous loop
IF ( ALLOCATED( TimeRefArray ) ) DEALLOCATE( TimeRefArray )
ALLOCATE( TimeRefArray( l , 6 ) )
TimeRefArray(:,1:5) = tmp2D(j:k,1:5)
TimeRefArray(:,6) = tmp2D_singlenumber(j:k)
TimeRefArray(:,1) = TimeRefArray(:,1) - TimeRefArray(1,1) ! offset the time information [decimal days] used later on in the netCDF file
PRINT *, "start = ", TimeRefArray(1,:)
PRINT *, "end = ", TimeRefArray(l,:)

CALL CPU_TIME(cpuTe)
PRINT '("CPU timing for reference date/time vector generation = ",F10.1," s")',cpuTe-cpuTs

DEALLOCATE(tmp2D, tmp2D_singlenumber)

!-------------------------------------------------------------------------------

CONTAINS

  INTEGER FUNCTION CheckForLeapyear( year )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: year

    IF ( ( MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0 ) .OR. ( MOD(year, 400) == 0 ) ) THEN
      CheckForLeapyear = 366
    ELSE
      CheckForLeapyear = 365
    END IF

  END FUNCTION CheckForLeapyear

!-------------------------------------------------------------------------------

END SUBROUTINE CreateRefTimeArray
