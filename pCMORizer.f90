!-------------------------------------------------------------------------------
! passing allocatable arrays between main program and external subroutine

MODULE FilelistHandling

  IMPLICIT NONE
  SAVE

  CHARACTER (len = 200):: tmpfileFL_std, tmpfileFL_xtrm, tmpfileFL_3d
  CHARACTER (len = 200), DIMENSION(:), ALLOCATABLE :: &
    fl_wrfout,  &
    fl_wrfxtr,  & 
    fl_wrfpres, &
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

END MODULE RefTimeVecs

!-------------------------------------------------------------------------------
! namelist handling

MODULE NamelistHandling

  IMPLICIT NONE
  SAVE

! reading from runctrl.current.nml 
  INTEGER, PARAMETER :: nvars = 39 ! 39 maximum number of vars per namelist, keep const at max number

  CHARACTER (len = 300) :: activity_id, contact, Conventions, domain_id, att_domain, &
  	driving_experiment_id,driving_experiment, &
  	driving_institution_id, driving_source_id, driving_variant_label, grid, institution, &
  	institution_id, license, mip_era, product, project_id, source, &
  	source_id, source_type, version, version_realization, references, tracking_id, &
	variable_id

  CHARACTER (len = 1000) :: comment

  CHARACTER (len = 300) :: DirInputSimResRoot, DirOutputPostProRoot, domain, calendar
  
! Josipa: Add npl variables to read number of pressure levels in wrfpress file 
  INTEGER ::  nx, ny, nz, npl, xoffset, yoffset, xfocus, yfocus, nfiles
  CHARACTER (len = 4) :: ts, te
  CHARACTER (len = 19) :: tstot, tetot
  CHARACTER (len = 300) :: PnFnGeo
! Josipa: Add projection variable to distingush rotated from lambert  
  CHARACTER (len = 3) :: projection

  LOGICAL :: aggregation_yearly, aggregation_monthly, aggregation_individually
  CHARACTER (len = 19) :: tsact, teact
  INTEGER :: nvar
                     

! reading from runctrl.vars.nml*
  CHARACTER (LEN = 100), DIMENSION(nvars) :: var_wrf, var_cmip, standard_name, &
    long_name, units, filetype, cmfx, cm1hr, cm3hr, cm6hr, cmDay, cmMon, cmSea, positive
  INTEGER, DIMENSION(nvars):: height, plevel, cordexID
  LOGICAL, DIMENSION(nvars):: time1hr, time3hr, time6hr, timeDay, timeMon, timeSea, &
     interpolate, timefx

  NAMELIST / globalvars / activity_id, contact, Conventions, domain_id, att_domain, &
  driving_experiment_id,driving_experiment, &
  driving_institution_id, driving_source_id, driving_variant_label, grid, institution, &
  institution_id, license, mip_era, product, project_id, source, &
  source_id, source_type, version, version_realization, references, &
  comment

  NAMELIST / filesystem / DirInputSimResRoot, DirOutputPostProRoot, domain, calendar

  NAMELIST / model_config / ts, te, nx, ny, nz, npl, xoffset, yoffset, xfocus, &
    yfocus, tstot, tetot

  NAMELIST / static_fields / PnFnGeo, projection

  NAMELIST / tool_config / nvar, aggregation_yearly, aggregation_monthly, &
    aggregation_individually, tsact, teact
    
  NAMELIST / vars / var_wrf, var_cmip, standard_name, long_name, units, &
    plevel, height, time1hr, time3hr, time6hr, timeDay, timeMon, timeSea, timefx, &
    filetype, cmfx, cm1hr, cm3hr, cm6hr, cmDay, cmMon, cmSea, interpolate, cordexID, positive

END MODULE NamelistHandling

!===============================================================================

PROGRAM WRFCMORizer

USE FilelistHandling
USE RefTimeVecs
USE NamelistHandling

USE netcdf

#ifndef SERIAL
USE MPI
#endif
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

  SUBROUTINE calcslptwo(slp, PP, P_s, PHI_s, T_L, ns, ew)
    IMPLICIT NONE
    REAL, DIMENSION(:,:),   INTENT(OUT) :: slp
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: PP
    REAL, DIMENSION(:,:),   INTENT(IN)  :: P_s
    REAL, DIMENSION(:,:),   INTENT(IN)  :: PHI_s
    REAL, DIMENSION(:,:),   INTENT(IN)  :: T_L
    INTEGER,                INTENT(IN)  :: ns, ew
  END SUBROUTINE calcslptwo

! Josipa: Add subrotine to calculate wint at different levels 
  SUBROUTINE var_zwind(nz, u, v, z, u10, v10, sa, ca, newz, unewz, vnewz)
    IMPLICIT NONE
    INTEGER, INTENT(in)                 :: nz
    REAL, DIMENSION(nz), INTENT(in)     :: u,v,z
    REAL, INTENT(in)                    :: u10, v10, sa, ca, newz
    REAL, INTENT(out)                   :: unewz, vnewz
  END SUBROUTINE var_zwind 
  
  
! Josipa: Add subrotine for linear interpolation between levels 
  SUBROUTINE linear_int(nz, var3d, z, var2d, newz, varout)
    IMPLICIT NONE
    INTEGER, INTENT(in)                 :: nz
    REAL, DIMENSION(nz), INTENT(in)     :: var3d, z
    REAL, INTENT(in)                    :: var2d, newz
    REAL, INTENT(out)                   :: varout
  END SUBROUTINE linear_int 
  
  
END INTERFACE

!===============================================================================
! filenames
! do not use the individual variable namelists here, they are set later on in a 
! loop; create symbolic links to have different base namelists, depending on the
! dataset and experiment

CHARACTER (LEN=*), PARAMETER :: fnNMLexp = "runctrl.current.nml"
CHARACTER (LEN=*), PARAMETER :: fnVersion = "VERSION.txt"
CHARACTER (LEN=100), DIMENSION(:), ALLOCATABLE :: fnNMLvar
CHARACTER (LEN=200) :: pn_out, fn_out, iflWRFin
CHARACTER (LEN=100) :: cmorizer_version

! Choose how psl is calculated, 0 OK, 1 OK, 2 OK, no not change throughout
INTEGER :: calc_slp_type = 2 

! (HTr): Choose how geopotantial height and hus is vertically interpolated to given pressure levels
INTEGER :: vint_type = 1

!-------------------------------------------------------------------------------
! constant base temperature and pressure
REAL, PARAMETER :: T00 = 300.0
REAL, PARAMETER :: P00 = 100000.0

!-------------------------------------------------------------------------------
! constants
REAL, PARAMETER :: cp = 1004.0      ! [J kg-1 K-1]
REAL, PARAMETER :: R = 287.04       ! [J kg-1 K-1]
REAL, PARAMETER :: L = 2501000.0    ! [J kg-1]
REAL, PARAMETER :: a = 610.78       ! [Pa]
REAL, PARAMETER :: b = 17.27        ! 17.2693882
REAL, PARAMETER :: c = 273.15  
REAL, PARAMETER :: d = 35.86  
REAL, PARAMETER :: n = L*0.622*a/cp 
REAL, PARAMETER :: mv = 1.e20       ! missing value as specified
REAL, PARAMETER :: gr = 9.81
REAL, PARAMETER :: epsil = 0.6220
REAL, PARAMETER :: erad=6370000.

! constants for hurs calc
REAL(KIND=8), PARAMETER :: a1=-2.8365744e3,a2=-6.028076559e3,a3=19.54263612,a4=-2.737830188e-2
REAL(KIND=8), PARAMETER :: a5=1.6261698e-5,a6=7.0229056e-10,a7=-1.8680009e-13,a8=2.7150305
REAL(KIND=8), PARAMETER :: b1=-5.8666426e3,b2=22.32870244,b3=1.39387003e-2,b4=-3.4262402e-5
REAL(KIND=8), PARAMETER :: b5=2.7040955e-8,b6=6.7063522e-1
REAL(KIND=8), PARAMETER :: fi1=3.62183e-4,fi2=2.6061244e-5,fi3=3.8667770e-7,fi4=3.8268958e-9,fi5=-10.7604,fi6=6.3987441e-2,fi7=-2.6351566e-4,fi8=1.6725084e-6
REAL(KIND=8), PARAMETER :: fw1=3.536240e-4,fw2=2.932836e-5,fw3=2.616898e-7,fw4=8.581361e-9,fw5=-10.75880,fw6=6.326813e-2,fw7=-2.536893e-4,fw8=6.340529e-7

! new netCDF file
INTEGER :: ncid, ncidin, ncidin0
INTEGER :: lon_dimid, lat_dimid, rec_dimid, height_dimid, &
  nb2_dimid, x_dimid, y_dimid, plev_dimid, depth_dimid

INTEGER :: varid, x_varid, lon_varid, lat_varid, rlon_varid, rlat_varid, hgt_varid, &
  rotated_pole_varid, lambert_varid, height_varid, rec_varid, pp_varid, pb_varid, ph_varid, &
  phb_varid, pl_varid, pl_varid_u, pl_varid_v, qv_varid, qc_varid, qi_varid, qr_varid, &
  qg_varid, qh_varid, qs_varid, theta_varid, t2_varid, recbnds_varid, rainnc_varid, &
  rainc_varid, u10_varid, v10_varid, u_varid, v_varid, w_varid, &
  sfcevp_varid, potevp_varid, sfroff_varid, udroff_varid, acsnom_varid, q2_varid, &
  sinalpha_varid, cosalpha_varid, plev_varid, plevbnds_varid, psfc_varid, &
  depth_varid, soillayerbnds_varid, cd_varid, xlon_varid, ylat_varid, t00_varid, p00_varid, &
  mask_varid, sh2o_varid

! input data general query
INTEGER :: ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in
INTEGER :: InVarIdRec, InDimLenRec 
INTEGER :: nvar_nml 

! n of lowest vertical levels to calculate wind speed at heights > 10m, to avoid loading complete set of levels
INTEGER, PARAMETER :: nz_lowest = 10

CHARACTER (LEN=NF90_MAX_NAME) :: InDimNameRec 
CHARACTER (LEN=50) :: fl_filter
CHARACTER (LEN=19), DIMENSION(:), ALLOCATABLE :: InVarDataRec

!-------------------------------------------------------------------------------
! data
REAL :: &
  ! geographical metadata info
  GeoNPLon      , &
  GeoNPLat      , &
  GeoLat1       , &
  GeoLat2       , &
  GeoCenLon     , &
  GeoCenLat     , &
  dx_distance   , &  
  dy_distance   , &
  ! bucket variables
  bucket_mm     , &
  bucket_J      , &
  t_ii          , & 
  t_tmp         , & 
  !variable to fill the missing data in variabels extracted from wrfpress
  p_miss        , &          
  ! vertical interpolation
  slope         , &
  low_lev       , &
  high_lev

! Time vec stuff
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: &
  TimeRefArraySubsetMean
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: &
  TimeRefArraySubset      , &   
  Time_bnds               , &
  tmp1_2d                 , &
  tmp3_2d                 , &
  tmp4_2d                 , &
  data_in  

! 1D variables 
REAL, DIMENSION(:), ALLOCATABLE :: &
  GeoInRLat     , & 
  GeoInRLon     , &
  GeoInYLat     , & 
  GeoInXLon     , &
  Std_parallel

! 2D variables 
REAL, DIMENSION(:,:), ALLOCATABLE :: &
  psl_in        , &
  t2_in         , &
  q2_in         , &
  cldfra_inv    , &
  u10_in        , &
  v10_in        , &
  cape          , &
  cin           , &
  li            , &
  lcl           , &
  lfc           , &
  prw           , &
  clwvi         , &
  clivi         , &
  clgvi         , &
  clhvi         , &
  sinalpha_in   , &
  cosalpha_in   , &
  psfc_in       , &
  rainc_max_in  , &
  rainnc_max_in , &
  landmask_in   , &
  var2d_in      , &
  var2d_out     , &
  var2d_u       , &
  var2d_v       , &
  hgt_in        , &       
  cd_in         , &
  tmp_2d        , &
  tmp2_2d       

! 3D variables 
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
  qg_in       , &
  qh_in       , &
  qs_in       , &
  theta_in    , &
  t_in        , &
  ph_fl       , &
  p_in        , &
  u_in        , &
  v_in        , &
  w_in        , &
  var3d_in    , &
  var3d_in_u  , &
  var3d_in_v  , &
  potevp_in   , &
  rainnc_in   , &
  rainc_in    , &
  rad_in      , &
  hflux_in    , &
  t_p         , &
  acsnom_in   , &
  GeoInLonLat , &
  sfcevp_in   , &
  sfroff_in   , &
  udroff_in   , &
  swdown_in   , &
  tmp_3d      , &
  pl_in       , &
  pl_in_u     , &
  pl_in_v     , &
  smois_in    , & 
  sh2o_in    , & 
  tslb_in  


! 4D variables
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: &
  cldfra_in

! Soil data dimension
REAL, DIMENSION(:), ALLOCATABLE :: DZS ! soil layer thickness may vary from simulation to simulation
REAL, DIMENSION(4) :: DZShc = (/0.1, 0.3, 0.6, 1.0/)

! pressure levels [Pa]; Shorter set: REAL, DIMENSION(6) :: pout = (/100000.,92500.,85000.,70000.,50000.,20000./)
INTEGER, PARAMETER :: num_levels = 6
REAL, DIMENSION(num_levels) :: pout = (/1000.,925.,850.,700.,500.,200./)

! bucket system
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
  i_rainnc_in, i_rainc_in, i_rad_in

!-------------------------------------------------------------------------------
! time and date handling
REAL :: dtHours
REAL(KIND=8) :: tsact_singlenumber, teact_singlenumber
REAL, DIMENSION(:), ALLOCATABLE :: InDateTimeCombined
CHARACTER (LEN=4) :: InDateTimeYearStr, tsactYearStr, teactYearStr
CHARACTER (LEN=3), DIMENSION(:), ALLOCATABLE :: frequency 
CHARACTER (LEN=2) :: InDateTimeMonthStr, FirstHourStr, FirstMinuteStr, LastDayStr, &
  LastHourStr, LastMinuteStr, tsactMonthStr, tsactDayStr, tsactHourStr, tsactMinuteStr, &
  teactMonthStr, teactDayStr, teactHourStr, teactMinuteStr
CHARACTER (LEN=12) :: FileNameStartDateTime, FileNameEndDateTime
CHARACTER (LEN=100), DIMENSION(nvars) :: cell_methods
INTEGER :: InDateTimeYearPrev = 0, InDateTimeMonthPrev = 0
INTEGER :: tsactYear, tsactMonth, tsactDay, tsactHour, tsactMinute, tsactSecond, & 
  teactYear, teactMonth, teactDay, teactHour, teactMinute, teactSecond
INTEGER, DIMENSION(:), ALLOCATABLE :: InDateTimeYear, InDateTimeMonth, &
  InDateTimeDay, InDateTimeHour, InDateTimeMinute, InDateTimeSecond 
LOGICAL, DIMENSION(nvars) :: procflag

! some special switch for input handling
INTEGER :: i, j, k, sts, ivar, ifrq, ifl, it, counter, np, nl, ii, ivarnml, prevpass = 0, iteration
LOGICAL :: FileExists, newpass, time_match, calc = .TRUE., inputtimesteptruncate = .FALSE.
REAL :: cpuTs, cpuTe
INTEGER, ALLOCATABLE :: ipos(:)

#ifndef SERIAL
CHARACTER (LEN = 3) :: FileNrStr
#endif
!-------------------------------------------------------------------------------
! system calls for tracking ID and date
CHARACTER (LEN=*), PARAMETER :: cmdUUID = "uuidgen -t > tmpfileUUID"
CHARACTER (LEN=37) :: trackingID

CHARACTER (LEN=*), PARAMETER :: cmdDate = "date -u +%Y-%m-%dT%H:%M:%SZ > tmpfileDate"
CHARACTER (LEN=20) :: creationDate

!===============================================================================

! MPI
#ifndef SERIAL
INTEGER, DIMENSION(:), ALLOCATABLE :: ivar_list
INTEGER :: ierr, rank, numtasks

! Initialise MPI evironment
  CALL MPI_INIT(ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_INIT"

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_COMM_RANK"

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "MPI_COMM_WORLD"
#endif

PRINT *, "============================================================"
PRINT *, "WRF RCM CMORizer"
PRINT *, "============================================================"

OPEN(2,FILE=fnVersion,STATUS='old',ACTION='read',IOSTAT=sts)
  IF (sts /= 0) STOP "*** Could not open VERSION.txt, stopping***"
  READ(UNIT=2,FMT='(A)',IOSTAT=sts) cmorizer_version
CLOSE(2)
IF (sts /= 0) STOP "*** Could not read VERSION.txt, stopping***"
cmorizer_version = TRIM(ADJUSTL(cmorizer_version))

PRINT *, "============================================================"
PRINT *, "*** CENTRAL NAMELIST READING ***"
PRINT *, fnNMLexp

OPEN(2,FILE=fnNMLexp)
READ(UNIT=2,NML=globalvars)
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

! this is a dirty hack for the mrsol, mrsfl and tsl variables
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
ALLOCATE( landmask_in(xfocus, yfocus) )

#ifndef SERIAL
sts = NF90_OPEN(TRIM(PnFnGeo), IOR(NF90_NOWRITE, NF90_MPIIO), ncidin, &
   comm = MPI_COMM_WORLD, info = MPI_INFO_NULL )
#else
sts = NF90_OPEN(TRIM(PnFnGeo), NF90_NOWRITE, ncidin)
#endif

sts = NF90_INQ_VARID(ncidin, "XLONG_M", varid)
sts = NF90_GET_VAR(ncidin, varid, GeoInLonLat(:, :, 1), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /))

sts = NF90_INQ_VARID(ncidin, "XLAT_M", varid)
sts = NF90_GET_VAR(ncidin, varid, GeoInLonLat(:, :, 2), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /))



IF ( projection == "LCC" ) THEN
  ALLOCATE( Std_parallel(2) )
  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "TRUELAT1", Std_parallel(1))
  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "TRUELAT2", Std_parallel(2))
  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "CEN_LON", GeoCenLon)
  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "CEN_LAT", GeoCenLat)

  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "DX", dx_distance)
  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "DY", dy_distance)  
     
  ALLOCATE( GeoInYLat(yfocus) )
  ALLOCATE( GeoInXLon(xfocus) )
  DO i = 1,xfocus
    GeoInXLon(i) = i*dx_distance - dx_distance/2 - (xfocus)*dx_distance/2
  END DO
  DO j = 1,yfocus
    GeoInYLat(j) = j*dy_distance - dy_distance/2 - (yfocus)*dy_distance/2
  END DO

ELSE
  ALLOCATE( GeoInRLat(yfocus) )
  ALLOCATE( GeoInRLon(xfocus) )

  sts = NF90_INQ_VARID(ncidin, "CLONG", varid)
  sts = NF90_GET_VAR(ncidin, varid, GeoInRLon(:), &
    START = (/ xoffset, 1, 1 /), COUNT = (/ xfocus, 1, 1 /))

  sts = NF90_INQ_VARID(ncidin, "CLAT", varid)
  sts = NF90_GET_VAR(ncidin, varid, GeoInRLat(:), &
    START = (/ 1, yoffset, 1 /), COUNT = (/ 1, yfocus, 1 /))

  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "POLE_LAT", GeoNPLat)
  sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "POLE_LON", GeoNPLon)
  IF ( GeoNPLon > 0.0 ) THEN 
    GeoNPLon = GeoNPlon -180.0
  END IF
  
END IF

IF (.not. ALLOCATED(sinalpha_in)) ALLOCATE( sinalpha_in( xfocus, yfocus ), STAT=sts )
sts = NF90_INQ_VARID(ncidin, "SINALPHA", sinalpha_varid)
sts = NF90_GET_VAR(ncidin, sinalpha_varid, sinalpha_in(:,:), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )

IF (.not. ALLOCATED(cosalpha_in)) ALLOCATE( cosalpha_in( xfocus, yfocus ), STAT=sts )
sts = NF90_INQ_VARID(ncidin, "COSALPHA", cosalpha_varid)
sts = NF90_GET_VAR(ncidin, cosalpha_varid, cosalpha_in(:,:), &
  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )

sts = NF90_CLOSE(ncidin)

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

ALLOCATE ( fnNMLvar(1) )
fnNMLvar(1) = "runctrl.vars.nml" 

!-------------------------------------------------------------------------------
 DO ifrq = 1, 1, 1 ! 1hr

  PRINT *, "============================================================"
  PRINT *, "freq = ", frequency(ifrq)
  
  PRINT *, "============================================================"
  PRINT *, "*** FILELIST CREATION ***"
  
  tmpfileFL_std = "tmpfileFLstd" // TRIM(domain)
  tmpfileFL_3d = "tmpfileFL3d" // TRIM(domain)
  tmpfileFL_xtrm = "tmpfileFLxtrm" // TRIM(domain)

  ! creates a year range, can be expanded also to use months
  IF ( (ts == "0000") .AND. (te == "0000") ) THEN
    fl_filter = ""
  ELSE
    fl_filter = ts
  END IF
  PRINT *, "fl_filter from runctrl nml: ", fl_filter

#ifndef SERIAL
  IF ( rank == 0 ) THEN
    CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfout*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*' | sort > " // tmpfileFL_std)
  END IF
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
#else
  CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfout*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*' | sort > " // tmpfileFL_std)
#endif
  ft = 0 ! file type
  CALL GenerateFilelist

#ifndef SERIAL
  IF ( rank == 0 ) THEN 
    CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfxtrm*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*' | sort > " // tmpfileFL_xtrm)
  END IF
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
#else
  CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfxtrm*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*' | sort > " // tmpfileFL_xtrm)
#endif
  ft = 1
  CALL GenerateFilelist

#ifndef SERIAL
  IF ( rank == 0) THEN
    CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfpress*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*' | sort > " // tmpfileFL_3d)
  END IF
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
#else
  CALL SYSTEM("find " // TRIM(DirInputSimResRoot) // "/ -name 'wrfpress*" // TRIM(domain) // "*" // TRIM(fl_filter) // "*' | sort > " // tmpfileFL_3d)
#endif
  ft = 2
  CALL GenerateFilelist
 
  !DO i=1,SIZE(fl_wrfout(:)),1
    !PRINT '(100A)', fl_wrfout(i)
  !END DO
  !DO i=1,SIZE(fl_wrfxtr(:)),1
    !PRINT '(100A)', fl_wrfxtr(i)
  !END DO
  !DO i=1,SIZE(fl_wrfpres(:)),1
    !PRINT '(100A)', fl_wrfpres(i)
  !END DO

  PRINT *, "============================================================"
  PRINT *, "*** TIME REFERENCE ARRAY ***"
  
  CALL CreateRefTimeArray( frequency(ifrq) )
  
  PRINT *, "size & shape of the TimRefArray = ", &
    SIZE(TimeRefArray), &
    SHAPE(TimeRefArray)
  PRINT *, "SIZE(TimeRefArray,1): ", SIZE(TimeRefArray,1) 

!-------------------------------------------------------------------------------
   DO ivarnml = 1, 1, 1  
    PRINT *, "============================================================"
    PRINT *, "var. namelist nr. and name: ", ivarnml, TRIM(fnNMLvar(ivarnml))

    ! read the very specific namelist
    OPEN(2,FILE=TRIM(fnNMLvar(ivarnml)))
      READ(UNIT=2,NML=vars)
    CLOSE(2)

    ! must have read the namelist already
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
    CASE ('fx')
      cell_methods(:) = cmfx(:)
      procflag(:) = timefx(:)
      dtHours = 1.
    CASE DEFAULT
      PRINT *, "invalid time interval specified"
      STOP
    END SELECT

    ! nvar_nml might be determined through the namelist itself
    ! sequence here does not matter
    SELECT CASE (TRIM(fnNMLvar(ivarnml)))
    CASE ("runctrl.vars.nml") ! OK
      nvar_nml = nvar
    END SELECT
 
    PRINT *, "number of vars inside current namelist: nvar_nml = ", nvar_nml
  
  !-------------------------------------------------------------------------------
  ! loop over all vars in the individual namelist
#ifndef SERIAL
    ALLOCATE( ivar_list(nvar_nml) )
    DO ivar = 1, nvar_nml, 1
      ivar_list(ivar) = ivar
    END DO

    CALL MPI_SCATTER( ivar_list, 1, MPI_INT, ivar, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    IF ( ierr /= MPI_SUCCESS ) STOP "MPI_SCATTER"
#else
    DO ivar = 1, nvar_nml, 1
#endif
  
      PRINT *,"============================================================"
      PRINT *, "*** ", TRIM(var_cmip(ivar)), procflag(ivar), " ***"
 
      IF (procflag(ivar)) THEN
 
      !-------------------------------------------------------------------------------
      ! loop over the filelists

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

        IF ( frequency(ifrq) == "fx" ) THEN
          nfiles = 1                ! Since the loop over files does not take the last file into account
        ELSE IF (aggregation_yearly) THEN
          nfiles = SIZE(fl_input)-1 ! To avoid creating output for the next year with only 1 day
        ELSE
          nfiles = SIZE(fl_input)
        END IF

        DO ifl = 1, nfiles, 1 ! operational: loop over complete filelist

          PRINT *,"============================================================"
          PRINT *, "filelist filetype = ", filetype(ivar) 
          PRINT *, "#files to process = ", SIZE(fl_input)

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

#ifndef SERIAL
        sts = NF90_OPEN(iflWRFin, IOR(NF90_NOWRITE, NF90_MPIIO), ncid_in, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
#else
	sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncid_in)
#endif
	
        sts = NF90_INQUIRE(ncid_in, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
        sts = NF90_INQ_VARID(ncid_in, "Times", InVarIdRec)
        sts = NF90_INQUIRE_DIMENSION(ncid_in, unlimdimid_in, &
          NAME = InDimNameRec, LEN = InDimLenRec)

        ! For the fixed variable, read only 2 timesteps and not 1 to make sure 
        ! that it does not take 1st timestep, as it can be read from the forcing files and filled with NaNs
	IF ( frequency(ifrq) == "fx" ) THEN
           InDimLenRec=1  ! 
	END IF

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

        ALLOCATE(InVarDataRec(InDimLenRec))
        sts = NF90_GET_VAR(ncid_in, InVarIdRec, InVarDataRec)
        sts = NF90_CLOSE(ncid_in)
        PRINT *, "RCM input file time coverage:"
        PRINT *, InVarDataRec(1), " to ", InVarDataRec(InDimLenRec)
 
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
        DO i = 1, SIZE(InVarDataRec), 1
        
          READ( InVarDataRec(i), '(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)' ) &
            InDateTimeYear(i), &
            InDateTimeMonth(i), &
            InDateTimeDay(i), &
            InDateTimeHour(i), &
            InDateTimeMinute(i), &
            InDateTimeSecond(i)
            
          ! 2009062008 use for later-on comparisons
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
            
            !PRINT *,'SIZE(TimeRefArray, 1)', SIZE(TimeRefArray, 1)
            !PRINT *,'SHAPE(TimeRefArray, 1)', SHAPE(TimeRefArray, 1)
            !PRINT *,'TimeRefArray(1,2)', TimeRefArray(1,2)
            !PRINT *,'TimeRefArray(744,2)', TimeRefArray(744,2)
            !PRINT *,'InDateTimeYear(it)', InDateTimeYear(it)

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
              PRINT *,tsactYear, tsactMonth, tsactDay, tsactHour
              PRINT *,teactYear, teactMonth, teactDay, teactHour

              tsact_singlenumber = REAL(tsactYear,8)*1000000._8 + REAL(tsactMonth,8)*10000._8 + REAL(tsactDay,8)*100._8 + REAL(tsactHour,8)
              teact_singlenumber = REAL(teactYear,8)*1000000._8 + REAL(teactMonth,8)*10000._8 + REAL(teactDay,8)*100._8 + REAL(teactHour,8)
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
              !IF ( ( cell_methods(ivar) == "mean" ) .OR. &
              !     ( cell_methods(ivar) == "sum" ) .OR. &
              !     ( cell_methods(ivar) == "minimum" ) .OR. &
              !     ( cell_methods(ivar) == "maximum" ) ) THEN
              !  counter = counter - 1
              !  ipos = ipos(1:counter)
              !  PRINT *, "size ipos", SIZE(ipos)
              !  PRINT *, "counter = ", counter
              !END IF     

            ELSE IF (aggregation_monthly) THEN
              PRINT *, "aggregation_monthly"
              DO i = 1, SIZE(TimeRefArray, 1), 1
                IF (( TimeRefArray(i,2) == InDateTimeYear(it)) .AND. &
                   ( TimeRefArray(i,3) == InDateTimeMonth(it))) THEN
                  counter = counter + 1
                  ipos = [ipos, i]
                END IF
              END DO
 
            ELSE IF (aggregation_yearly) THEN ! default, CORDEX annual files
              PRINT *, "aggregation_annually"
              PRINT *, InDateTimeYear(it)
              DO i = 1, SIZE(TimeRefArray, 1), 1
                IF ( TimeRefArray(i,2) == InDateTimeYear(it)) THEN
                  counter = counter + 1 
                  ipos = [ipos, i]
                END IF
              END DO
            END IF

            IF ( frequency(ifrq) == "fx" ) THEN
              counter = 1
            ELSE
              counter = counter
            END IF

            PRINT *, "timesteps in the time ref. subset = ", counter

            ALLOCATE( TimeRefArraySubset( counter, 6 ) ) ! index, y, m, d, h
            ALLOCATE( TimeRefArraySubsetMean( counter ) )
            ALLOCATE( Time_bnds( 2, counter ) )
            DO i = 1, counter, 1
              j = ipos(i)
              TimeRefArraySubset(i,1:6) = TimeRefArray(j,1:6)
            END DO
            !PRINT '(F9.3,1X,F5.0,1X,F3.0,1X,F3.0,1X,F3.0)', TRANSPOSE( TimeRefArraySubset(:,:) )
            !PRINT *, TimeRefArraySubset(0,:)
            
            PRINT *, "Start_year=",TimeRefArraySubset(1,2) ! index, y, m, d, h
            PRINT *, "End_date=",TimeRefArraySubset(counter,1:6) ! index, y, m, d, h

            DEALLOCATE(ipos)

            !-------------------------------------------------------------------------------
            ! create path- and filenames according to the ruleset of the CORDEX data 
            ! protocol: 2 main options here:
            ! Non standard:
            ! TODO: not all is implemented as the averaging is not yet implemented

            PRINT *, "generate path and filename"

            ! only a single file is created and filled successively
            ! use the information from the namelist
            ! 2009-06-20_00:00:00, may be full days or end e.g. at 23UTC
            ! arbitrary timespan
            IF (aggregation_individually) THEN

              WRITE (tsactYearStr,'(I4.4)') INT(TimeRefArraySubset(1,2))
              WRITE (tsactMonthStr,'(I2.2)') INT(TimeRefArraySubset(1,3))    
              WRITE (tsactDayStr,'(I2.2)') INT(TimeRefArraySubset(1,4))   
              WRITE (tsactHourStr,'(I2.2)') INT(tsactHour)
              !WRITE (tsactMinuteStr,'(I2.2)') INT(TimeRefArraySubset(1,6))
              WRITE (teactYearStr,'(I4.4)') INT(TimeRefArraySubset(counter,2))
              WRITE (teactMonthStr,'(I2.2)') INT(TimeRefArraySubset(counter,3))
              WRITE (teactDayStr,'(I2.2)') INT(TimeRefArraySubset(counter,4))   
              WRITE (teactHourStr,'(I2.2)') INT(teactHour) 
              !WRITE (teactHourStr,'(I2.2)') INT(TimeRefArraySubset(counter,5))
              !WRITE (teactMinuteStr,'(I2.2)') INT(TimeRefArraySubset(counter,6))

              IF ((frequency(ifrq) == '1hr') .OR. &
                  (frequency(ifrq) == '3hr') .OR. &
                  (frequency(ifrq) == '6hr')) THEN
                  
              IF ( (cell_methods(ivar) == "mean") .OR. (cell_methods(ivar) == "sum") ) THEN 
                !WRITE (tsactHourStr,'(I2.2)') INT( FLOOR( ((dtHours/2.)*60.) / 60. ) )
                WRITE (tsactMinuteStr,'(I2.2)') INT( MOD( (dtHours/2.)*60., 60. ) )
                !WRITE (teactHourStr,'(I2.2)') INT( FLOOR( ( (24.*60.) - (dtHours/2.)*60.)  / 60. ) )
                WRITE (teactMinuteStr,'(I2.2)') INT( MOD( (24.*60.) - (dtHours/2.)*60., 60. ) )
              ELSE
                !WRITE (tsactHourStr,'(I2.2)') INT(TimeRefArraySubset(1,5))
                WRITE (tsactMinuteStr,'(I2.2)') INT(0)
                !WRITE (teactHourStr,'(I2.2)') INT(TimeRefArraySubset(counter,5))
                WRITE (teactMinuteStr,'(I2.2)') INT(0)
              END IF
              
                FileNameStartDateTime = tsactYearStr//tsactMonthStr//tsactDayStr//tsactHourStr//tsactMinuteStr
                FileNameEndDateTime = teactYearStr//teactMonthStr//teactDayStr//teactHourStr//teactMinuteStr
                PRINT *, "DATE_START == ", FileNameStartDateTime
                PRINT *, "DATE_END == ", FileNameEndDateTime

              ELSE IF (frequency(ifrq) == 'day') THEN
  
                FileNameStartDateTime = tsactYearStr//tsactMonthStr//tsactDayStr
                FileNameEndDateTime = teactYearStr//teactMonthStr//teactDayStr

              ELSE IF ((frequency(ifrq) == 'mon') .OR. &
                       (frequency(ifrq) == 'sem')) THEN

                FileNameStartDateTime = tsactYearStr//tsactMonthStr
                FileNameEndDateTime = teactYearStr//teactMonthStr

              ELSE IF (frequency(ifrq) == 'fx') THEN 
                FileNameStartDateTime = ''
                FileNameEndDateTime = ''

              END IF

            ! determined by the date information in the wrf outputs, automatic
            ! monthly aggregation is also special
            ELSE IF ( (aggregation_monthly) .OR. (aggregation_yearly) ) THEN

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
                !Josipa: add minutes to the hourly files to avoid differences in file nameing between cumulative and non-cumulative variabels. 
                FirstHourStr = "00"
                WRITE (LastHourStr,'(I2.2)') 24-INT(dtHours)
		FirstMinuteStr = "00"
		LastMinuteStr = "00"
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
                    FileNameStartDateTime = InDateTimeYearStr//InDateTimeMonthStr//"01"//FirstHourStr//FirstMinuteStr
                    FileNameEndDateTime = InDateTimeYearStr//InDateTimeMonthStr//LastDayStr//LastHourStr//LastMinuteStr
                  END IF  

                ELSE IF (frequency(ifrq) == 'day') THEN
  
                  FileNameStartDateTime = InDateTimeYearStr//InDateTimeMonthStr//"01"
                  FileNameEndDateTime = InDateTimeYearStr//InDateTimeMonthStr//LastDayStr

                ELSE IF (frequency(ifrq) == 'fx') THEN 
                  FileNameStartDateTime = ''
                  FileNameEndDateTime = ''

                END IF

              ELSE IF (aggregation_yearly) THEN 

                IF ((frequency(ifrq) == '1hr') .OR. &
                    (frequency(ifrq) == '3hr') .OR. &
                    (frequency(ifrq) == '6hr')) THEN
                  IF ( (cell_methods(ivar) == "mean") .OR. (cell_methods(ivar) == "sum") ) THEN
                    FileNameStartDateTime = InDateTimeYearStr//"0101"//FirstHourStr//FirstMinuteStr
                    FileNameEndDateTime = InDateTimeYearStr//"1231"//LastHourStr//LastMinuteStr
                  ELSE
                    FileNameStartDateTime = InDateTimeYearStr//"0101"//FirstHourStr//FirstMinuteStr
                    FileNameEndDateTime = InDateTimeYearStr//"1231"//LastHourStr//LastMinuteStr
                  END IF

                ! TODO: 5 yearly
                ELSE IF (frequency(ifrq) == 'day') THEN  
                  FileNameStartDateTime = InDateTimeYearStr//"0101"
                  FileNameEndDateTime = InDateTimeYearStr//"1231"

                ! TODO: 10 yearly
                ELSE IF ((frequency(ifrq) == 'mon') .OR. &
                         (frequency(ifrq) == 'sem')) THEN
                  STOP "not yet implemented functionality, monthly/seasonally"
                  !FileNameStartDateTime = tsactYearStr//tsactMonthStr
                  !FileNameEndDateTime = teactYearStr//teactMonthStr

                ELSE IF (frequency(ifrq) == 'fx') THEN 
                  FileNameStartDateTime = ''
                  FileNameEndDateTime = ''

                END IF
              END IF
            END IF

            ! /hpc/shared/int/eva/ramod_WRF_CRPGL/WRFrv021rXXrcc3CpCdx/postpro/
            ! EUR-44/CRPGL/ECMWF-ERAINT/evaluation/r1i1p1/CRPGL-WRFARW331/v1
            pn_out = TRIM(project_id)                   // "/" // &
                     TRIM(mip_era)                      // "/" // &
                     TRIM(activity_id)                 	// "/" // &
                     TRIM(domain_id)                  	// "/" // &
                     TRIM(institution_id)              	// "/" // &
                     TRIM(driving_source_id)        	// "/" // &
                     TRIM(driving_experiment_id)    	// "/" // &
                     TRIM(driving_variant_label)        // "/" // &
                     TRIM(source_id)                	// "/" // &
                     TRIM(version_realization)      	// "/" // &
                     TRIM(frequency(ifrq))             	// "/" // &
                     TRIM(var_cmip(ivar))		// "/" // &
                     TRIM(version)

	    IF (frequency(ifrq) == 'fx') THEN
            	fn_out = TRIM(var_cmip(ivar))            // "_" // &
                     TRIM(domain_id)                 	 // "_" // &
                     TRIM(driving_source_id)             // "_" // &
                     TRIM(driving_experiment_id)         // "_" // &
                     TRIM(driving_variant_label)         // "_" // &
                     TRIM(institution_id)                // "_" // &
                     TRIM(source_id)                     // "_" // &
                     TRIM(version_realization)           // "_" // &
                     TRIM(frequency(ifrq))		// &
                     ".nc"
            ELSE

            	fn_out = TRIM(var_cmip(ivar))            // "_" // &
                     TRIM(domain_id)                 	 // "_" // &
                     TRIM(driving_source_id)             // "_" // &
                     TRIM(driving_experiment_id)         // "_" // &
                     TRIM(driving_variant_label)         // "_" // &
                     TRIM(institution_id)                // "_" // &
                     TRIM(source_id)                     // "_" // &
                     TRIM(version_realization)           // "_" // &
                     TRIM(frequency(ifrq))               // "_" // &
                     TRIM(FileNameStartDateTime) // "-" // TRIM(FileNameEndDateTime) // &
                     ".nc"
            END IF
  
            PRINT *, "CORDEX compliant pathname, pn_out = ", TRIM(pn_out)
            PRINT *, "CORDEX compliant filename, fn_out = ", TRIM(fn_out)

            PRINT *, "-----------------------------------------------------------------------------"
            PRINT *, "*** CHECK FOR FILE EXISTANCE / CREATE FILE WITH BASIC STRUCTURE + METADATA***"

            INQUIRE( FILE=TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) // "/" // TRIM(fn_out), EXIST=FileExists )

            IF ( FileExists ) THEN
              PRINT *, "++++ path and file exist, continue filling"
            ELSE
              PRINT *, "++++ path and file do not yet exist, create path and netCDF file first"
              PRINT '(150A)', "path = ", TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out)

              CALL SYSTEM("mkdir -p " // TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) )
 
#ifndef SERIAL
              WRITE(FileNrStr,'(i3)') rank
              CALL SYSTEM("uuidgen -t > tmpfileUUID"//TRIM(domain)//TRIM(fnNMLvar(ivarnml))//TRIM(ADJUSTL(FileNrStr)))
              OPEN(1,FILE="tmpfileUUID"//TRIM(domain)//TRIM(fnNMLvar(ivarnml))//TRIM(ADJUSTL(FileNrStr)),STATUS='old')
              READ(1,*) trackingID
              CLOSE(1)
              PRINT *, "uuidgen externally generated trackingID = ", trackingID

              CALL SYSTEM("date -u +%Y-%m-%dT%H:%M:%SZ > tmpfileDate"//TRIM(domain)//TRIM(fnNMLvar(ivarnml))//TRIM(ADJUSTL(FileNrStr)))
              OPEN(1,FILE="tmpfileDate"//TRIM(domain)//TRIM(fnNMLvar(ivarnml))//TRIM(ADJUSTL(FileNrStr)),STATUS='old')
              READ(1,*) creationDate
              CLOSE(1)
              PRINT *, "date externally generated creation date = ", creationDate

              CALL mpi_barrier(MPI_COMM_WORLD, ierr)
              IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
#else

              CALL SYSTEM("uuidgen -t > tmpfileUUID"//TRIM(domain)//TRIM(fnNMLvar(ivarnml)))
              OPEN(1,FILE="tmpfileUUID"//TRIM(domain)//TRIM(fnNMLvar(ivarnml)),STATUS='old')
              READ(1,*) trackingID
              CLOSE(1)
              PRINT *, "uuidgen externally generated trackingID = ", trackingID

              CALL SYSTEM("date -u +%Y-%m-%dT%H:%M:%SZ > tmpfileDate"//TRIM(domain)//TRIM(fnNMLvar(ivarnml)))
              OPEN(1,FILE="tmpfileDate"//TRIM(domain)//TRIM(fnNMLvar(ivarnml)),STATUS='old')
              READ(1,*) creationDate
              CLOSE(1)
              PRINT *, "date externally generated creation date = ", creationDate
#endif

              PRINT *, "create netCDF file"
              sts = NF90_CREATE(TRIM(DirOutputPostProRoot) // "/" // &
                TRIM(pn_out) // "/" // TRIM(fn_out), IOR(NF90_NETCDF4, &
                NF90_CLASSIC_MODEL), ncid)

              !-----------------------------------------------------------------
              ! define time dimension
              IF (frequency(ifrq) /= 'fx') THEN
                sts = NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, rec_dimid)
                IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                     ( cell_methods(ivar) == "sum" ) .OR. &
                     ( cell_methods(ivar) == "minimum" ) .OR. &
                     ( cell_methods(ivar) == "maximum" ) ) THEN
                  sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
                ENDIF
              END IF

              ! define horizontal dimensions depending on the projection      
              IF ( projection == "LCC" ) THEN 
                sts = NF90_DEF_DIM(ncid, "x", xfocus, x_dimid)
                sts = NF90_DEF_DIM(ncid, "y", yfocus, y_dimid)
              ELSE
                sts = NF90_DEF_DIM(ncid, "rlon", xfocus, lon_dimid)
                sts = NF90_DEF_DIM(ncid, "rlat", yfocus, lat_dimid) 
              END IF
              
              ! define vertical dimension for near-surface variables      
              IF ( height(ivar) /= -999 ) THEN
                sts = NF90_DEF_DIM(ncid, "height", 1, height_dimid)
              END IF
              
              ! define vertical dimension for variable on pressure levels    
              IF ( ( plevel(ivar) /= -999 ) ) THEN
                sts = NF90_DEF_DIM(ncid, "plev", 1, plev_dimid)
              END IF

              ! define vertical dimension for variable on pressure levels    
              IF (TRIM(var_cmip(ivar)) == 'mrsos' ) THEN
                sts = NF90_DEF_DIM(ncid, "sdepth", 1, depth_dimid)
                sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
              END IF

              ! define special vertical dimension dimension 4D vars
              IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. &
                  (TRIM(var_cmip(ivar)) == 'mrsfl') .OR. &
                  (TRIM(var_cmip(ivar)) == 'tsl')) THEN
                sts = NF90_DEF_DIM(ncid, "sdepth", SIZE(DZShc), depth_dimid)
                sts = NF90_DEF_DIM(ncid, "bnds", 2, nb2_dimid)
              ENDIF


              ! Add coordinates for LCC and rotated grid projections

              IF ( projection == "LCC" ) THEN
              	! included for for lamber conformal projection 

                sts = nf90_def_var(ncid, "x", NF90_DOUBLE, (/ x_dimid /), xlon_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, xlon_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, xlon_varid, "standard_name", "projection_x_coordinate")
                sts = nf90_put_att(ncid, xlon_varid, "long_name", "X Coordinate Of Projection")
                sts = nf90_put_att(ncid, xlon_varid, "units", "m")
                sts = nf90_put_att(ncid, xlon_varid, "axis", "X")
  
                sts = nf90_def_var(ncid, "y", NF90_DOUBLE, (/ y_dimid /), ylat_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, ylat_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, ylat_varid, "standard_name", "projection_y_coordinate")
                sts = nf90_put_att(ncid, ylat_varid, "long_name", "Y Coordinate Of Projection")
                sts = nf90_put_att(ncid, ylat_varid, "units", "m")
                sts = nf90_put_att(ncid, ylat_varid, "axis", "Y")
 
                sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ x_dimid, y_dimid /), lon_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, lon_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
                sts = nf90_put_att(ncid, lon_varid, "long_name", "Longitude")
                sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")
                sts = nf90_put_att(ncid, lon_varid, "_CoordinateAxisType", "Lon") ! special addon, not needed, but allowed

                sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ x_dimid, y_dimid /), lat_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, lat_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
                sts = nf90_put_att(ncid, lat_varid, "long_name", "Latitude")
                sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north")
                sts = nf90_put_att(ncid, lat_varid, "_CoordinateAxisType", "Lat") ! special addon, not needed, but allowed

                sts = nf90_def_var(ncid, "crs", NF90_CHAR, lambert_varid)
                sts = nf90_put_att(ncid, lambert_varid, "grid_mapping_name", "lambert_conformal_conic")
                sts = nf90_put_att(ncid, lambert_varid, "standard_parallel", Std_parallel)
                sts = nf90_put_att(ncid, lambert_varid, "longitude_of_central_meridian", GeoCenLon)
                sts = nf90_put_att(ncid, lambert_varid, "latitude_of_projection_origin", GeoCenLat)
                sts = nf90_put_att(ncid, lambert_varid, "false_easting", "0.")
                sts = nf90_put_att(ncid, lambert_varid, "false_northing", "0.")
                sts = nf90_put_att(ncid, lambert_varid, "earth_radius", erad)
                
              ELSE              
              	! included for for rotated projection  
                sts = nf90_def_var(ncid, "rlon", NF90_DOUBLE, (/ lon_dimid /), rlon_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, rlon_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, rlon_varid, "standard_name", "grid_longitude")
                sts = nf90_put_att(ncid, rlon_varid, "long_name", "Longitude in rotated pole grid")
                sts = nf90_put_att(ncid, rlon_varid, "units", "degrees")
                sts = nf90_put_att(ncid, rlon_varid, "axis", "X")
  
                sts = nf90_def_var(ncid, "rlat", NF90_DOUBLE, (/ lat_dimid /), rlat_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, rlat_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, rlat_varid, "standard_name", "grid_latitude")
                sts = nf90_put_att(ncid, rlat_varid, "long_name", "Latitude in rotated pole grid")
                sts = nf90_put_att(ncid, rlat_varid, "units", "degrees")
                sts = nf90_put_att(ncid, rlat_varid, "axis", "Y")

                sts = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lon_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, lon_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, lon_varid, "standard_name", "longitude")
                sts = nf90_put_att(ncid, lon_varid, "long_name", "Longitude")
                sts = nf90_put_att(ncid, lon_varid, "units", "degrees_east")

                sts = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ lon_dimid, lat_dimid /), lat_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, lat_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, lat_varid, "standard_name", "latitude")
                sts = nf90_put_att(ncid, lat_varid, "long_name", "Latitude")
                sts = nf90_put_att(ncid, lat_varid, "units", "degrees_north") 

		            sts = nf90_def_var(ncid, "crs", NF90_CHAR, rotated_pole_varid)
		            !sts = nf90_put_att(ncid, rotated_pole_varid, "long_name", "Coordinates of the rotated North Pole")
		            sts = nf90_put_att(ncid, rotated_pole_varid, "grid_mapping_name", "rotated_latitude_longitude")
		            sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_latitude", GeoNPLat)
		            sts = nf90_put_att(ncid, rotated_pole_varid, "grid_north_pole_longitude", GeoNPLon)
                sts = nf90_put_att(ncid, rotated_pole_varid, "earth_radius", erad)

		! additional and useful
		! important for remapping with cdo (conservative remapping only 
		! possible with this information)
		!        vertices = 4 ;
		!        float lat_vertices(rlat, rlon, vertices) ;
		!                lat_vertices:units = "degrees_north" ;
		!        float lon_vertices(rlat, rlon, vertices) ;
		!                lon_vertices:units = "degrees_east" ;
            
              END IF

              ! included for near surface variables at some height
              IF ( height(ivar) /= -999 ) THEN
                sts = nf90_def_var(ncid, "height", NF90_DOUBLE, (/ height_dimid /), height_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, height_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, height_varid, "standard_name", "height")
                sts = nf90_put_att(ncid, height_varid, "long_name", "height")
                sts = nf90_put_att(ncid, height_varid, "units", "m")
                sts = nf90_put_att(ncid, height_varid, "positive", "up")
                sts = nf90_put_att(ncid, height_varid, "axis", "Z")
              END IF

              ! included for variables on some pressure level
              IF ( plevel(ivar) /= -999 ) THEN
                sts = nf90_def_var(ncid, "plev", NF90_DOUBLE, (/ plev_dimid /), plev_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, plev_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, plev_varid, "standard_name", "air_pressure")
                sts = nf90_put_att(ncid, plev_varid, "long_name", "pressure")
                sts = nf90_put_att(ncid, plev_varid, "units", "Pa")
                sts = nf90_put_att(ncid, plev_varid, "positive", "down")
                sts = nf90_put_att(ncid, plev_varid, "axis", "Z")
                IF ( cell_methods(ivar) == "vmean" ) THEN ! if this is layers over which there has been some everaging
                  sts = nf90_put_att(ncid, plev_varid, "bounds", "plev_bnds")
                END IF
              END IF

              ! included for soil moisture and soil temperature at all levels   
              IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. &
                  (TRIM(var_cmip(ivar)) == 'mrsfl') .OR. &
                  (TRIM(var_cmip(ivar)) == 'tsl')) THEN
                sts = nf90_def_var(ncid, "sdepth", NF90_DOUBLE, (/ depth_dimid /), depth_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, depth_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, depth_varid, "standard_name", "depth")
                sts = nf90_put_att(ncid, depth_varid, "long_name", "Soil layer depth")
                sts = nf90_put_att(ncid, depth_varid, "units", "m")
                sts = nf90_put_att(ncid, depth_varid, "positive", "down")
                sts = nf90_put_att(ncid, depth_varid, "axis", "Z")
                sts = nf90_put_att(ncid, depth_varid, "bounds", "sdepth_bnds")
                sts = nf90_def_var(ncid, "sdepth_bnds", NF90_DOUBLE, (/ nb2_dimid, depth_dimid /), soillayerbnds_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, soillayerbnds_varid, 1, 1, 1)
              ENDIF
 
              ! included variabels averaged between levels  
              IF ( cell_methods(ivar) == "vmean" ) THEN
                sts = nf90_def_var(ncid, "plev_bnds", NF90_DOUBLE, (/ nb2_dimid /), plevbnds_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, plevbnds_varid, 1, 1, 1)
              END IF
  
              ! included always
              IF (frequency(ifrq) /= 'fx') THEN
                sts = nf90_def_var(ncid, "time", NF90_DOUBLE, (/ rec_dimid /), rec_varid, fletcher32 = .true.)
                sts = nf90_def_var_deflate(ncid, rec_varid, 1, 1, 1)
                sts = nf90_put_att(ncid, rec_varid, "standard_name", "time")
                sts = nf90_put_att(ncid, rec_varid, "long_name", "Time")
                sts = nf90_put_att(ncid, rec_varid, "units", "days since " // tstot(1:10) // "T" // tstot(12:19) // "Z" )
                sts = nf90_put_att(ncid, rec_varid, "calendar", calendar)
                sts = nf90_put_att(ncid, rec_varid, "axis", "T")
                IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                     ( cell_methods(ivar) == "sum" ) .OR. &
                     ( cell_methods(ivar) == "minimum" ) .OR. &
                     ( cell_methods(ivar) == "maximum" ) ) THEN
                  sts = nf90_put_att(ncid, rec_varid, "bounds", "time_bnds")
                END IF
  
                IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                     ( cell_methods(ivar) == "sum" ) .OR. &
                     ( cell_methods(ivar) == "minimum" ) .OR. &
                     ( cell_methods(ivar) == "maximum" ) ) THEN
                  sts = nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ nb2_dimid, rec_dimid /), recbnds_varid, fletcher32 = .true.)
                  sts = nf90_def_var_deflate(ncid, recbnds_varid, 1, 1, 1)
                END IF
              END IF        

              !-----------------------------------------------------------------              
              ! global attributes always included
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "WRF_CMORizer_version", cmorizer_version)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "activity_id", activity_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "contact", contact)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "Conventions", Conventions)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "creation_date", creationDate)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "domain", att_domain)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "domain_id", domain_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment", driving_experiment)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_experiment_id", driving_experiment_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_institution_id", driving_institution_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_source_id", driving_source_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "driving_variant_label", driving_variant_label)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "frequency", frequency(ifrq))
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "grid", grid)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution", institution)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "institution_id", institution_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "license", license)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "mip_era", mip_era)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "product", product)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "project_id", project_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "references", references)                    
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source", source)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source_id", source_id)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "source_type", source_type)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "tracking_id","hdl:21.14103/" //trackingID)
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "variable_id", var_cmip(ivar))
              sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "comment", comment)
              !sts = NF90_PUT_ATT(ncid, NF90_GLOBAL, "version_realization", version_realization)

              !-----------------------------------------------------------------
              ! always included -- definition of the individual variable
              ! compression is always on variable level, use compression with 
              ! level = 1 (=> CORDEX protocol),  
              ! compression on/off: 19MB->12MB, 0.1s->0.5s
              ! chunking: needs some careful considerations
              ! nf90_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
             
              ! Define dimensions for the LCC and ROT projections
              IF ( projection == "LCC" ) THEN
                IF ((var_cmip(ivar) == "mrsol") .OR. (var_cmip(ivar) == "mrsfl") .OR. (var_cmip(ivar) == "tsl")) THEN
                  sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ x_dimid, y_dimid, depth_dimid, rec_dimid /), x_varid, fletcher32 = .true.)
                ELSE IF (frequency(ifrq) == 'fx') THEN
                  sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ x_dimid, y_dimid/), x_varid, fletcher32 = .true.)
                ELSE
                  sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ x_dimid, y_dimid, rec_dimid /), x_varid, fletcher32 = .true.)
                END IF                 
              ELSE
                IF ((var_cmip(ivar) == "mrsol") .OR. (var_cmip(ivar) == "mrsfl") .OR. (var_cmip(ivar) == "tsl")) THEN
                  sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, depth_dimid, rec_dimid /), x_varid, fletcher32 = .true.)
                ELSE IF (frequency(ifrq) == 'fx') THEN
                  sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid/), x_varid, fletcher32 = .true.)
                ELSE
                  sts = nf90_def_var(ncid, var_cmip(ivar), NF90_FLOAT, (/ lon_dimid, lat_dimid, rec_dimid /), x_varid, fletcher32 = .true.)
                END IF
              END IF            

              ! Define variable metadata         
              sts = nf90_put_att(ncid, x_varid, "standard_name", standard_name(ivar))
              sts = nf90_put_att(ncid, x_varid, "long_name", long_name(ivar))
              sts = nf90_put_att(ncid, x_varid, "units", units(ivar))
              IF ( positive(ivar) /= '-999' ) THEN
                sts = nf90_put_att(ncid, x_varid, "positive", positive(ivar))
              END IF
              
              IF ( ( var_cmip(ivar) == "mrro" ) .OR. ( var_cmip(ivar) == "mrros" ) ) THEN
                sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar))//" area: "//TRIM(cell_methods(ivar))//" where land")
              ELSE
                sts = nf90_put_att(ncid, x_varid, "cell_methods", "time: "//TRIM(cell_methods(ivar)))
              END IF
              
              IF ( height(ivar) /= -999 ) THEN
                sts = nf90_put_att(ncid, x_varid, "coordinates", "height lat lon")
              ELSE IF ( plevel(ivar) /= -999 ) THEN 
                sts = nf90_put_att(ncid, x_varid, "coordinates", "plev lat lon") 
              ELSE IF ( (TRIM(var_cmip(ivar)) == 'mrsol')   .OR. &
                        (TRIM(var_cmip(ivar)) == 'mrsfl')   .OR. &
                	(TRIM(var_cmip(ivar)) == 'tsl'  ) ) THEN
                sts = nf90_put_att(ncid, x_varid, "coordinates", "sdepth lat lon") 
              ELSE
                sts = nf90_put_att(ncid, x_varid, "coordinates", "lat lon")
              END IF

              sts = nf90_put_att(ncid, x_varid, "grid_mapping", "crs")

              ! Set missing value and fill_value for the variable
              sts = nf90_put_att(ncid, x_varid, "missing_value", mv)
              sts = nf90_put_att(ncid, x_varid, "_FillValue", mv)
              
              !-----------------------------------------------------------------
              sts = nf90_def_var_deflate(ncid, x_varid, 1, 1, 1)
              sts = NF90_ENDDEF(ncid)

              !-----------------------------------------------------------------
              ! Fill time dimension for instantaneous variables
              IF ( cell_methods(ivar) == "point" ) THEN  
                sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySubset(:,1) )
              END IF

              ! Fill time dimension for averaged variables                          
              IF ( ( cell_methods(ivar) == "mean" ) .OR. &
                   ( cell_methods(ivar) == "sum" ) .OR. &
                   ( cell_methods(ivar) == "minimum" ) .OR. &
                   ( cell_methods(ivar) == "maximum" ) ) THEN
                TimeRefArraySubsetMean (:) = TimeRefArraySubset(:,1) + ( 0.5_8 * (1._8 / (24._8/dtHours) ) )
                sts = NF90_PUT_VAR(ncid, rec_varid, TimeRefArraySubsetMean(:) )                 
                Time_bnds(1,:) = TimeRefArraySubset(:,1)
                Time_bnds(2,:) = TimeRefArraySubset(:,1) + ( 1._8 * (1._8 / (24._8/dtHours) )) 
                sts = NF90_PUT_VAR(ncid, recbnds_varid, Time_bnds(:,:), START = (/ 1, 1 /) , COUNT = (/ 2, SIZE(Time_bnds(1,:)) /) )
              END IF
              

              ! TODO
              ! define plev_bnds, needed mainly for clh clm cli

              !-----------------------------------------------------------------
              ! write coordinates and height/level information
              
              sts = NF90_PUT_VAR(ncid, lon_varid, GeoInLonLat(:,:,1), &
                START = (/ 1, 1, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
              sts = NF90_PUT_VAR(ncid, lat_varid, GeoInLonLat(:,:,2), &
                START = (/ 1, 1, 2 /), COUNT = (/ xfocus, yfocus, 1 /) )


              IF ( projection == "LCC" ) THEN     
                sts = NF90_PUT_VAR(ncid, xlon_varid, GeoInXLon )
                sts = NF90_PUT_VAR(ncid, ylat_varid, GeoInYLat )
              ELSE 
                sts = NF90_PUT_VAR(ncid, rlon_varid, GeoInRLon )
                sts = NF90_PUT_VAR(ncid, rlat_varid, GeoInRLat )
              END IF  

              IF ( height(ivar) /= -999 ) THEN
                sts = NF90_PUT_VAR(ncid, height_varid, height(ivar) )
              END IF

              IF ( plevel(ivar) /= -999 ) THEN
                sts = NF90_PUT_VAR(ncid, plev_varid, plevel(ivar)*100. )
              END IF

              IF (TRIM(var_cmip(ivar)) == 'mrsos') THEN
                sts = NF90_PUT_VAR(ncid, depth_varid, (/ 0.05D0 /) ) 
                sts = NF90_PUT_VAR(ncid, soillayerbnds_varid, RESHAPE((/ 0.0D0, 0.1D0/), (/2,1/)))
              END IF

              IF ((TRIM(var_cmip(ivar)) == 'mrsol') .OR. (TRIM(var_cmip(ivar)) == 'mrsfl') .OR. (TRIM(var_cmip(ivar)) == 'tsl')) THEN
                sts = NF90_PUT_VAR(ncid, depth_varid, (/ 0.05D0, 0.25D0, 0.7D0, 1.5D0 /) ) ! center points of the depth layers, hardcoded
                sts = NF90_PUT_VAR(ncid, soillayerbnds_varid, RESHAPE((/ 0.0D0, 0.1D0, 0.1D0, 0.4D0, 0.4D0, 1.0D0, 1.0D0, 2.0D0 /), (/2,4/))) ! bnds hardcoded
              END IF 

              !-----------------------------------------------------------------
  
              sts = NF90_CLOSE(ncid)

              !-----------------------------------------------------------------
  
            END IF ! file exists y/n

          !-------------------------------------------------------------------------------
  
          ELSE ! checking whether this is the first pass or previously run
            PRINT *, "This is a subsequent pass, continue filling: read, procecess and sort data"
          END IF
  
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
  
            PRINT *, "index, where the WRF data will be located in the nc file = ", counter
  
          !-------------------------------------------------------------------------------
          ! read orig WRF outputs
          ! there is always a corresponding time-slot in the NC file
          ! extracted time from above
          ! "it" controls it all: timestep in the individual WRF file
          ! there is only one variable at a time under processing

  
	       ! Open the file
            sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin)

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            PRINT *, "*** READING OF VARIABLES ***"
            PRINT *, "variable to work on = ",TRIM(var_cmip(ivar))

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
                .OR.  (var_cmip(ivar) == "prw") &
                .OR.  (var_cmip(ivar) == "cin") &
                .OR.  (var_cmip(ivar) == "clivi") &
                .OR.  (var_cmip(ivar) == "clgvi") &
                .OR.  (var_cmip(ivar) == "clhvi") &
                .OR.  (var_cmip(ivar) == "clwvi") &
                .OR.  (var_cmip(ivar) == "cape") &
                .OR.  (var_cmip(ivar) == "cin") &
                .OR.  (var_cmip(ivar) == "li") &
                .OR.  ( (plevel(ivar) /= -999) .AND. (filetype(ivar) == "s") ) ) THEN

              PRINT *,'read 3D vars'

              !-------------------------------------------------------------------
              ! internal vars needed in the pressure level / 3D processing section

              ! always needed
              
              IF ( (plevel(ivar) /= -999) .AND. ( filetype(ivar) == "s") ) THEN
                PRINT *, "prep. int. 3D pres. level vars"
                IF (.not. ALLOCATED(var3d_in)) ALLOCATE( var3d_in( xfocus, yfocus, nz ), STAT=sts )
                IF ( ( INDEX(var_cmip(ivar),"ua") == 1 ) .OR. ( INDEX(var_cmip(ivar),"va") == 1 ) ) THEN
                  IF (.not. ALLOCATED(var3d_in_u)) ALLOCATE( var3d_in_u( xfocus, yfocus, nz ), STAT=sts )
                  IF (.not. ALLOCATED(var3d_in_v)) ALLOCATE( var3d_in_v( xfocus, yfocus, nz ), STAT=sts )
                END IF
              END IF
              
              PRINT *, "allocate p_in, t_in, ph_fl" 
              IF (.not. ALLOCATED(p_in))  ALLOCATE( p_in( xfocus, yfocus, nz ), STAT=sts ) 
              IF (.not. ALLOCATED(t_in))  ALLOCATE( t_in( xfocus, yfocus, nz ), STAT=sts )
              IF (.not. ALLOCATED(ph_fl)) ALLOCATE( ph_fl( xfocus, yfocus, nz ), STAT=sts )

              IF ( var_cmip(ivar) == "prw" ) THEN
                IF (.not. ALLOCATED(prw)) ALLOCATE( prw( xfocus, yfocus ), STAT=sts )
              END IF

              IF ( var_cmip(ivar) == "psl" ) THEN
                IF (.not. ALLOCATED(psl_in)) ALLOCATE( psl_in ( xfocus, yfocus ), STAT=sts )
              END IF

              IF ( (var_cmip(ivar) == "cape" ) .OR. (var_cmip(ivar) == "cin" ) .OR. &
                   (var_cmip(ivar) == "li" ) ) THEN
                IF (.not. ALLOCATED(t_p))  ALLOCATE( t_p( xfocus, yfocus, nz ), STAT=sts )
                IF (.not. ALLOCATED(qvs))  ALLOCATE( qvs( xfocus, yfocus, nz ), STAT=sts )
                IF (.not. ALLOCATED(cape)) ALLOCATE( cape( xfocus, yfocus ), STAT=sts )
                IF (.not. ALLOCATED(cin))  ALLOCATE( cin( xfocus, yfocus ), STAT=sts )
                IF (.not. ALLOCATED(lcl))  ALLOCATE( lcl( xfocus, yfocus ), STAT=sts )
                IF (.not. ALLOCATED(lfc))  ALLOCATE( lfc( xfocus, yfocus ), STAT=sts )
                IF (.not. ALLOCATED(li))   ALLOCATE( li( xfocus, yfocus ), STAT=sts )
              END IF

              IF ( var_cmip(ivar) == "clwvi" ) THEN
                IF (.not. ALLOCATED(clwvi)) ALLOCATE( clwvi( xfocus, yfocus ), STAT=sts )
              ENDIF

              IF ( var_cmip(ivar) == "clivi" ) THEN
                IF (.not. ALLOCATED(clivi)) ALLOCATE( clivi( xfocus, yfocus ), STAT=sts )
              END IF

              IF ( var_cmip(ivar) == "clgvi" ) THEN
                IF (.not. ALLOCATED(clgvi)) ALLOCATE( clgvi( xfocus, yfocus ), STAT=sts )
              END IF

              IF ( var_cmip(ivar) == "clhvi" ) THEN
                IF (.not. ALLOCATED(clhvi)) ALLOCATE( clhvi( xfocus, yfocus ), STAT=sts )
              END IF

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
            !-------------------------------------------------------------------

              IF ( ( var_cmip(ivar) == "prw" ) &
                 .OR. ( var_cmip(ivar) == "psl" ) &
                 .OR. ( INDEX(var_cmip(ivar),"hus") == 1 ) &
                 .OR. ( var_cmip(ivar) == "cape" ) &
                 .OR. ( var_cmip(ivar) == "cin" )  &
                 .OR. ( var_cmip(ivar) == "li") ) THEN
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
  
              IF ( ( var_cmip(ivar) == "clwvi" ) &
                 .OR. ( var_cmip(ivar) == "clivi" ) ) THEN
                PRINT *, "read QSNOW"
                IF (.not. ALLOCATED(qs_in)) ALLOCATE( qs_in( xfocus, yfocus, nz  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "QSNOW", qs_varid)
                sts = NF90_GET_VAR(ncidin, qs_varid, qs_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
              END IF

              IF ( var_cmip(ivar) == "clgvi" ) THEN
                PRINT *, "read QGRAUP"
                IF (.not. ALLOCATED(qg_in)) ALLOCATE( qg_in( xfocus, yfocus, nz  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "QGRAUP", qg_varid)
                sts = NF90_GET_VAR(ncidin, qg_varid, qg_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
              END IF

              IF ( var_cmip(ivar) == "clhvi" ) THEN
                PRINT *, "read QHAIL"
                IF (.not. ALLOCATED(qh_in)) ALLOCATE( qh_in( xfocus, yfocus, nz  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "QHAIL", qh_varid)
                sts = NF90_GET_VAR(ncidin, qh_varid, qh_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )
              END IF

              IF ( (INDEX(var_cmip(ivar),"ua") == 1 ) .OR. ( INDEX(var_cmip(ivar),"va") == 1 ) ) THEN
                PRINT *, "read U"
                IF (.not. ALLOCATED(u_in)) ALLOCATE( u_in( xfocus+1, yfocus, nz ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "U", u_varid)
                sts = NF90_GET_VAR(ncidin, u_varid, u_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus+1, yfocus, nz, 1 /) )
  
                PRINT *, "read V"
                IF (.not. ALLOCATED(v_in)) ALLOCATE( v_in( xfocus, yfocus+1, nz ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "V", v_varid)
                sts = NF90_GET_VAR(ncidin, v_varid, v_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus+1, nz, 1 /) )
              END IF
            
              IF ( (INDEX(var_cmip(ivar),"wa") == 1) ) THEN
                PRINT *, "read W"
                IF (.not. ALLOCATED(w_in)) ALLOCATE( w_in( xfocus, yfocus, nz+1 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "W", w_varid)
                sts = NF90_GET_VAR(ncidin, w_varid, w_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz+1, 1 /) )
              END IF
  
              IF ( var_cmip(ivar) == "psl" ) THEN
                PRINT *, "read PSFC"
                IF (.not. ALLOCATED(psfc_in)) ALLOCATE( psfc_in( xfocus, yfocus ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "PSFC", psfc_varid)
                sts = NF90_GET_VAR(ncidin, psfc_varid, psfc_in(:,:), &
                  START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
              END IF
              
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Reading data from wrfpress files
            ELSE IF ( filetype(ivar) == "p" )  THEN
              sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "P_LEV_MISSING", p_miss)
              IF ( INDEX(var_cmip(ivar),"hus") == 1 ) THEN
                PRINT *, "read Q_PL from wrfpress files"
                IF (.not. ALLOCATED(pl_in)) ALLOCATE( pl_in( xfocus, yfocus, npl  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "Q_PL", pl_varid)
                sts = NF90_GET_VAR(ncidin, pl_varid, pl_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, npl, 1 /) )
                  WHERE (pl_in < -900.) pl_in = mv
              END IF

              IF ( INDEX(var_cmip(ivar),"ta") == 1 ) THEN
                PRINT *, "read T_PL from wrfpress files"
                IF (.not. ALLOCATED(pl_in)) ALLOCATE( pl_in( xfocus, yfocus, npl  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "T_PL", pl_varid)
                sts = NF90_GET_VAR(ncidin, pl_varid, pl_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, npl, 1 /) )
                  WHERE  (pl_in < -900.) pl_in(:,:,:) = mv
              END IF

              IF ( INDEX(var_cmip(ivar),"zg") == 1 ) THEN
                PRINT *, "read GHT_PL from wrfpress files"
                IF (.not. ALLOCATED(pl_in)) ALLOCATE( pl_in( xfocus, yfocus, npl  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "GHT_PL", pl_varid)
                sts = NF90_GET_VAR(ncidin, pl_varid, pl_in(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, npl, 1 /) )
                  WHERE (pl_in < -900.) pl_in = mv
              END IF

              IF ( (INDEX(var_cmip(ivar),"ua") == 1 ) .OR. ( INDEX(var_cmip(ivar),"va") == 1) ) THEN
                PRINT *, "read U_PL from wrfpress files"
                IF (.not. ALLOCATED(pl_in_u)) ALLOCATE( pl_in_u( xfocus, yfocus, npl  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "U_PL", pl_varid_u)
                sts = NF90_GET_VAR(ncidin, pl_varid_u, pl_in_u(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, npl, 1 /) )
                  WHERE (pl_in_u < -900.) pl_in_u = mv
                PRINT *, "read V_PL from wrfpress files"
                IF (.not. ALLOCATED(pl_in_v)) ALLOCATE( pl_in_v( xfocus, yfocus, npl  ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin, "V_PL", pl_varid_v)
                sts = NF90_GET_VAR(ncidin, pl_varid_v, pl_in_v(:,:,:), &
                  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, npl, 1 /) )
                  WHERE (pl_in_v < -900.) pl_in_v = mv
              END IF

              IF ( var_cmip(ivar) == "od550aer" ) THEN
                PRINT *, "variable to read/write with no additional processing = ", var_wrf(ivar)
                sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
                sts = NF90_GET_VAR(ncidin, varid, data_in(:,:), &
                  START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
              END IF

              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! Josipa: ua and va at height > 10m [ms-1] wind components at heights higher then 10m 

              ELSE IF ( ( height(ivar) > 10. ) .AND. &
			( (INDEX(var_cmip(ivar),"ua") == 1 ) .OR.  (INDEX(var_cmip(ivar),"va") == 1 )  .OR. &
			  (INDEX(var_cmip(ivar),"ta") == 1 ) .OR.  (INDEX(var_cmip(ivar),"hus") == 1 ) ) )THEN

		PRINT *, "alocating 2D and 3D variables"
		IF (.not. ALLOCATED(hgt_in)) ALLOCATE( hgt_in( xfocus, yfocus ), STAT=sts )
		IF (.not. ALLOCATED(ph_in))  ALLOCATE( ph_in( xfocus, yfocus, nz_lowest + 1 ), STAT=sts )
		IF (.not. ALLOCATED(phb_in)) ALLOCATE( phb_in( xfocus, yfocus, nz_lowest + 1 ), STAT=sts )
		IF (.not. ALLOCATED(ph_fl))  ALLOCATE( ph_fl( xfocus, yfocus, nz_lowest ), STAT=sts )

		PRINT *, "read HGT"
		sts = NF90_INQ_VARID(ncidin, "HGT", hgt_varid)
		sts = NF90_GET_VAR(ncidin, hgt_varid, hgt_in(:,:), &
		  START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
		PRINT *, "read PH"
		sts = NF90_INQ_VARID(ncidin, "PH", ph_varid)
		sts = NF90_GET_VAR(ncidin, ph_varid, ph_in(:,:,:), &
		  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz_lowest+1, 1 /) )
		PRINT *, "read PHB"
		sts = NF90_INQ_VARID(ncidin, "PHB", phb_varid)
		sts = NF90_GET_VAR(ncidin, phb_varid, phb_in(:,:,:), &
		  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz_lowest+1, 1 /) )
		  
		  
		IF ((INDEX(var_cmip(ivar),"ua") == 1 ) .OR.  (INDEX(var_cmip(ivar),"va") == 1 )) THEN		  
		  
		  IF (.not. ALLOCATED(u_in))   ALLOCATE( u_in( xfocus+1, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(v_in))   ALLOCATE( v_in( xfocus, yfocus+1, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(u10_in)) ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts )
		  IF (.not. ALLOCATED(v10_in)) ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts ) 
		  IF (.not. ALLOCATED(var3d_in_u))  ALLOCATE( var3d_in_u ( xfocus, yfocus, nz_lowest), STAT=sts )
		  IF (.not. ALLOCATED(var3d_in_v))  ALLOCATE( var3d_in_v ( xfocus, yfocus, nz_lowest), STAT=sts )
		  IF (.not. ALLOCATED(sinalpha_in)) ALLOCATE( sinalpha_in( xfocus, yfocus ), STAT=sts )
		  IF (.not. ALLOCATED(cosalpha_in)) ALLOCATE( cosalpha_in( xfocus, yfocus ), STAT=sts )
		  IF (.not. ALLOCATED(var2d_u)) ALLOCATE( var2d_u( xfocus, yfocus ), STAT=sts )
		  IF (.not. ALLOCATED(var2d_v)) ALLOCATE( var2d_v( xfocus, yfocus ), STAT=sts )
		
            	  PRINT *, "read U"
            	  sts = NF90_INQ_VARID(ncidin, "U", u_varid)
            	  sts = NF90_GET_VAR(ncidin, u_varid, u_in(:,:,:), &
              	    START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus+1, yfocus, nz_lowest, 1 /) )
            	  PRINT *, "read V"
            	  sts = NF90_INQ_VARID(ncidin, "V", v_varid)
            	  sts = NF90_GET_VAR(ncidin, v_varid, v_in(:,:,:), &
              	    START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus+1, nz_lowest, 1 /) )
            	  PRINT *, "read U10"
            	  sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
            	  sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            	  PRINT *, "read V10"
            	  sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)  
            	  sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
                    
                ELSE IF (INDEX(var_cmip(ivar),"ta") == 1 ) THEN
                
		  IF (.not. ALLOCATED(var3d_in)) ALLOCATE( var3d_in( xfocus, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(theta_in)) ALLOCATE( theta_in( xfocus, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(pp_in)) ALLOCATE( pp_in( xfocus, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(pb_in)) ALLOCATE( pb_in( xfocus, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(p_in)) ALLOCATE( p_in( xfocus, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(var2d_in)) ALLOCATE( var2d_in( xfocus, yfocus), STAT=sts )
		  IF (.not. ALLOCATED(var2d_out)) ALLOCATE( var2d_out( xfocus, yfocus), STAT=sts )
             	 
             	  PRINT *, "read P"

                  sts = NF90_INQ_VARID(ncidin, "P", pp_varid)
                  sts = NF90_GET_VAR(ncidin, pp_varid, pp_in(:,:,:), &
                    START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz_lowest, 1 /) )

                  PRINT *, "read PB"
                  IF (.not. ALLOCATED(pb_in)) ALLOCATE( pb_in( xfocus, yfocus, nz_lowest ), STAT=sts )
                  sts = NF90_INQ_VARID(ncidin, "PB", pb_varid)
                  sts = NF90_GET_VAR(ncidin, pb_varid, pb_in(:,:,:), &
                    START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz_lowest, 1 /) )
		  
                  PRINT *, "read T"
            	  sts = NF90_INQ_VARID(ncidin, "T", theta_varid)
            	  sts = NF90_GET_VAR(ncidin, theta_varid, theta_in(:,:,:), &
              	    START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz_lowest, 1 /) )
              	    
                  PRINT *, "read T2"
            	  sts = NF90_INQ_VARID(ncidin, "T2", t2_varid)
            	  sts = NF90_GET_VAR(ncidin, t2_varid, var2d_in(:,:), &
              	    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
                                               				
                                               				
                ELSE IF (INDEX(var_cmip(ivar),"hus") == 1 ) THEN
                
		  IF (.not. ALLOCATED(var3d_in))   ALLOCATE( var3d_in( xfocus, yfocus, nz_lowest ), STAT=sts )
		  IF (.not. ALLOCATED(var2d_in))   ALLOCATE( var2d_in( xfocus, yfocus), STAT=sts )
		  IF (.not. ALLOCATED(var2d_out)) ALLOCATE( var2d_out( xfocus, yfocus), STAT=sts )
		  
                  PRINT *, "read QVAPOR"
            	  sts = NF90_INQ_VARID(ncidin, "QVAPOR", qv_varid)
            	  sts = NF90_GET_VAR(ncidin, qv_varid, var3d_in(:,:,:), &
              	    START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz_lowest, 1 /) )
              	    
                  PRINT *, "read T2"
            	  sts = NF90_INQ_VARID(ncidin, "Q2", q2_varid)
            	  sts = NF90_GET_VAR(ncidin, q2_varid, var2d_in(:,:), &
              	    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
              	    
              	END IF
      				
	    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! cll, clm, clh, clt [%] Cloud Fractions

            ELSE IF ( (var_cmip(ivar) == "cll") .OR. &
                  (var_cmip(ivar) == "clm") .OR.  &
                  (var_cmip(ivar) == "clh") .OR.  &
                  (var_cmip(ivar) == "clt") ) THEN
              
              IF (.not. ALLOCATED(pp_in)) ALLOCATE( pp_in( xfocus, yfocus, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "P", pp_varid)
              sts = NF90_GET_VAR(ncidin, pp_varid, pp_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )

              PRINT *, "read PB"
              IF (.not. ALLOCATED(pb_in)) ALLOCATE( pb_in( xfocus, yfocus, nz ), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "PB", pb_varid)
              sts = NF90_GET_VAR(ncidin, pb_varid, pb_in(:,:,:), &
                START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, nz, 1 /) )

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
                iflWRFin = fl_input(ifl+1) ! if not the first set to the previous wrfoutfile 

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
   
            ELSE IF ( (var_cmip(ivar) == "pr") .OR. ( (var_cmip(ivar) == "prhmax") .AND. (filetype(ivar) == "s" ) ) ) THEN 
             
              PRINT *, "read iflWRFin " , iflWRFin
              PRINT *, "it = ", it
              PRINT *, "inDimLenRec, number of output times in input = ", InDimLenRec

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
                  START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )

                sts = NF90_GET_ATT(ncidin, NF90_GLOBAL, "BUCKET_MM", bucket_mm)
                IF ( bucket_mm > 0. ) THEN
                  IF (.not. ALLOCATED(i_rainnc_in)) ALLOCATE( i_rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
                  sts = NF90_INQ_VARID(ncidin, "I_RAINNC", varid)
                  sts = NF90_GET_VAR(ncidin, varid, i_rainnc_in(:,:,1), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )   
                  IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
                  sts = NF90_INQ_VARID(ncidin, "I_RAINC", varid)
                  sts = NF90_GET_VAR(ncidin, varid, i_rainc_in(:,:,1), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
                END IF

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
		  IF (.not. ALLOCATED(i_rainnc_in)) ALLOCATE( i_rainnc_in ( xfocus, yfocus, 2 ), STAT=sts )
                  sts = NF90_INQ_VARID(ncidin0, "I_RAINNC", varid)
                  sts = NF90_GET_VAR(ncidin0, varid, i_rainnc_in(:,:,2), &
                    START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) ) 
                  IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )  
                  sts = NF90_INQ_VARID(ncidin0, "I_RAINC", varid)
                  sts = NF90_GET_VAR(ncidin0, varid, i_rainc_in(:,:,2), &
                    START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
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
		  IF (.not. ALLOCATED(i_rainc_in)) ALLOCATE( i_rainc_in ( xfocus, yfocus, 2 ), STAT=sts )
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
            !  prsn [kg m-2 s-1] ?>a Snowfall Flux or snm [kg m-2 s-1] a Surface Snow Melt
            ! [kg m-2 s-1]
            
            ELSE IF ( (var_cmip(ivar) == "prsn") .OR. (var_cmip(ivar) == "snm") ) THEN
          
	      PRINT *, "read XLAND"
              IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in( xfocus, yfocus), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "LANDMASK", mask_varid)
              sts = NF90_GET_VAR(ncidin, mask_varid, landmask_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
  
              IF (.not. ALLOCATED(acsnom_in)) ALLOCATE( acsnom_in ( xfocus, yfocus, 2 ), STAT=sts )
              IF (it /= InDimLenRec) THEN  
                sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), acsnom_varid)
                sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,:), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 2 /) )  
  
              ELSE IF ( (it == InDimLenRec) .and. (ifl /= nfiles) ) THEN  
                sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), acsnom_varid)
                sts = NF90_GET_VAR(ncidin, acsnom_varid, acsnom_in(:,:,1), &
                  START = (/ xoffset, yoffset, it /), &
                  COUNT = (/ xfocus, yfocus, 1 /))
  
                iflWRFin = fl_input(ifl+1)  
                sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)              
                sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), acsnom_varid)
                sts = NF90_GET_VAR(ncidin0, acsnom_varid, acsnom_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), &
                  COUNT =(/xfocus, yfocus, 1 /) )
            
                ! Close the next file
                sts = NF90_CLOSE(ncidin0)  
                ! Open the current file
                iflWRFin = fl_input(ifl)
              
              ELSE
                PRINT *, "no data available for average calculation any more"
                calc = .FALSE.  
  
              END IF
              
             !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             ! snc snw snd [kg m-2 s-1] a Surface Snow Melt
  
            ELSE IF ( (var_cmip(ivar) == "snc") .OR. &
           	 (var_cmip(ivar) == "snw") .OR. &
            	 (var_cmip(ivar) == "snd") )THEN
            
	      PRINT *, "read XLAND"
              IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in( xfocus, yfocus), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "LANDMASK", mask_varid)
              sts = NF90_GET_VAR(ncidin, mask_varid, landmask_in(:,:), &
               START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )


              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin, varid, data_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
  
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! evspsbl [kg m-2 s-1] a Evaporation  
          ELSE IF (var_cmip(ivar) == "evspsbl") THEN
 
            IF (.not. ALLOCATED(sfcevp_in)) ALLOCATE( sfcevp_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN  
              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), sfcevp_varid)
              sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )            
            ELSE IF ( (it == InDimLenRec) .AND. (ifl /= SIZE(fl_input)) ) THEN  
              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), sfcevp_varid)
              sts = NF90_GET_VAR(ncidin, sfcevp_varid, sfcevp_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 1/)) 
 
              iflWRFin = fl_input(ifl+1)   
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)  
              sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), sfcevp_varid)
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
          
            PRINT *, "read XLAND"
            IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in( xfocus, yfocus), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "LANDMASK", mask_varid)
              sts = NF90_GET_VAR(ncidin, mask_varid, landmask_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
  
            IF (.not. ALLOCATED(sfroff_in)) ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 2 /) )
  
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus,1/))
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
                sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), sfroff_varid)
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
          
            PRINT *, "read XLAND"
            IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in( xfocus, yfocus), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "LANDMASK", mask_varid)
              sts = NF90_GET_VAR(ncidin, mask_varid, landmask_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
  
            IF (.not. ALLOCATED(sfroff_in)) ALLOCATE( sfroff_in ( xfocus, yfocus, 2 ), STAT=sts )
            IF (.not. ALLOCATED(udroff_in)) ALLOCATE( udroff_in ( xfocus, yfocus, 2 ), STAT=sts )

            IF (it /= InDimLenRec) THEN
  
              sts = NF90_INQ_VARID(ncidin, "RUNSF", sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
  
              sts = NF90_INQ_VARID(ncidin, "RUNSB", udroff_varid)
              sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,:), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus, 2 /) )
  
            ELSE IF ( (it == InDimLenRec) .and. (ifl /= SIZE(fl_input)) ) THEN
  
              sts = NF90_INQ_VARID(ncidin, "RUNSF", sfroff_varid)
              sts = NF90_GET_VAR(ncidin, sfroff_varid, sfroff_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT = (/ xfocus, yfocus,1/))
  
              sts = NF90_INQ_VARID(ncidin, "RUNSB", udroff_varid)
              sts = NF90_GET_VAR(ncidin, udroff_varid, udroff_in(:,:,1), &
                START = (/ xoffset, yoffset, it /), &
                COUNT =(/xfocus,yfocus,1 /) )
  
              iflWRFin = fl_input(ifl+1)
  
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
  
              sts = NF90_INQ_VARID(ncidin0, "RUNSF", sfroff_varid)  
              sts = NF90_GET_VAR(ncidin0, sfroff_varid, sfroff_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), &
                COUNT =(/xfocus, yfocus, 1 /) )
 
              sts = NF90_INQ_VARID(ncidin0, "RUNSB", udroff_varid)
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
! rsdsdir [W m-2] a Surface Direct Downwelling Shortwave Radiation
! rsdscs [W m-2] a Surface Downwelling Clear-Sky Shortwave Radiation
! rldscs [W m-2] a Surface Downwelling Clear-Sky Longwave Radiation
! rsuscs [W m-2] a Surface Upwelling Clear-Sky Shortwave Radiation
! rluscs [W m-2] a Surface Upwelling Clear-Sky Longwave Radiation
! rsutcs [W m-2] a TOA Outgoing Clear-Sky Shortwave Radiation
! rlutcs [W m-2] a TOA Outgoing Clear-Sky Longwave Radiation
! only method in use finally: accumulated values with bucket
! restriction: assume bucket is used
            
          ELSE IF ( (var_cmip(ivar) == "rsds") &
            .OR. (var_cmip(ivar) == "rlds") &
            .OR. (var_cmip(ivar) == "rsus") &
            .OR. (var_cmip(ivar) == "rlus") &
            .OR. (var_cmip(ivar) == "rlut") &
            .OR. (var_cmip(ivar) == "rsdt") &
            .OR. (var_cmip(ivar) == "rsut") &
            .OR. (var_cmip(ivar) == "rsdsdir") &
            .OR. (var_cmip(ivar) == "rsdscs") &
            .OR. (var_cmip(ivar) == "rldscs") &
            .OR. (var_cmip(ivar) == "rsuscs") &
            .OR. (var_cmip(ivar) == "rluscs") &
            .OR. (var_cmip(ivar) == "rsutcs") &
            .OR. (var_cmip(ivar) == "rlutcs") ) THEN
            
            

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
                  START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
              END IF

              iflWRFin = fl_input(ifl+1)
              sts = NF90_OPEN(iflWRFin, NF90_NOWRITE, ncidin0)
              sts = NF90_INQ_VARID(ncidin0, TRIM(var_wrf(ivar)), varid)
              sts = NF90_GET_VAR(ncidin0, varid, rad_in(:,:,2), &
                START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1/) )
              sts = NF90_GET_ATT(ncidin0, NF90_GLOBAL, "BUCKET_J", bucket_J)
              IF ( bucket_J > 0. ) THEN
                PRINT *, "BUCKET_J = ", bucket_J
                IF (.not. ALLOCATED(i_rad_in)) ALLOCATE( i_rad_in ( xfocus, yfocus, 2 ), STAT=sts )
                sts = NF90_INQ_VARID(ncidin0, TRIM('I_' // var_wrf(ivar)), varid)
                sts = NF90_GET_VAR(ncidin0, varid, i_rad_in(:,:,2), &
                  START = (/ xoffset, yoffset, 1 /), COUNT = (/ xfocus, yfocus, 1 /) )
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
            .OR. (var_cmip(ivar) == "mrsos") &
            .OR. (var_cmip(ivar) == "mrsol") &
            .OR. (var_cmip(ivar) == "tsl")) THEN

	     PRINT *, "read XLAND"
             IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in( xfocus, yfocus), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "LANDMASK", mask_varid)
              sts = NF90_GET_VAR(ncidin, mask_varid, landmask_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
 
	     IF ( ( var_cmip(ivar) == "mrso" ) .OR. ( var_cmip(ivar) == "mrsos" ) .OR. &
	        ( var_cmip(ivar) == "mrsol" ) )THEN
 	    	PRINT *, "read SMOIS"
            	IF (.not. ALLOCATED(smois_in)) ALLOCATE( smois_in( xfocus, yfocus, 4), STAT=sts )           
            	sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
            	sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:), &
		  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, 4, 1 /) ) 
             	IF (.not. ALLOCATED(DZS)) ALLOCATE( DZS( SIZE(DZShc) ), STAT=sts )
                DZS(:) = DZShc(:)
                PRINT*,'DZS ', DZS(:)
             END IF  
             
	     IF ( ( var_cmip(ivar) == "tsl" ) )THEN
 	    	PRINT *, "read TSLB"
            	IF (.not. ALLOCATED(tslb_in)) ALLOCATE( tslb_in( xfocus, yfocus, 4), STAT=sts )           
            	sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
            	sts = NF90_GET_VAR(ncidin, varid, tslb_in(:,:,:), &
		  START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, 4, 1 /) ) 
             	IF (.not. ALLOCATED(DZS)) ALLOCATE( DZS( SIZE(DZShc) ), STAT=sts )
                DZS(:) = DZShc(:)
                PRINT*,'DZS ', DZS(:)
             END IF   
             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrfso, mrfsos, mrsfl [kg m-2] 
  
          ELSE IF ((var_cmip(ivar) == "mrfso") &
            .OR. (var_cmip(ivar) == "mrfsos") &
            .OR. (var_cmip(ivar) == "mrsfl")) THEN

	     PRINT *, "read XLAND"
             IF (.not. ALLOCATED(landmask_in)) ALLOCATE( landmask_in( xfocus, yfocus), STAT=sts )
              sts = NF90_INQ_VARID(ncidin, "LANDMASK", mask_varid)
              sts = NF90_GET_VAR(ncidin, mask_varid, landmask_in(:,:), &
                START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
 
 	     PRINT *, "read SMOIS"
             IF (.not. ALLOCATED(smois_in)) ALLOCATE( smois_in( xfocus, yfocus, 4), STAT=sts )           
             sts = NF90_INQ_VARID(ncidin, "SMOIS", varid)
             sts = NF90_GET_VAR(ncidin, varid, smois_in(:,:,:), &
               START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, 4, 1 /) ) 

 	     PRINT *, "read SH20"
             IF (.not. ALLOCATED(sh2o_in)) ALLOCATE( sh2o_in( xfocus, yfocus, 4), STAT=sts )           
             sts = NF90_INQ_VARID(ncidin, "SH2O", sh2o_varid)
             sts = NF90_GET_VAR(ncidin, sh2o_varid, sh2o_in(:,:,:), &
               START = (/ xoffset, yoffset, 1, it /), COUNT = (/ xfocus, yfocus, 4, 1 /) ) 
               
             IF (.not. ALLOCATED(DZS)) ALLOCATE( DZS( SIZE(DZShc) ), STAT=sts )
             DZS(:) = DZShc(:)
             PRINT*,'DZS ', DZS(:)
  
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
 
          ELSE IF ( (var_cmip(ivar) == "sfcWind") .OR. &
            ( (var_cmip(ivar) == "sfcWindmax") .AND. (filetype(ivar) == "s" ) ) ) THEN
  
            IF (.not. ALLOCATED(u10_in)) ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts ) 
            IF (.not. ALLOCATED(v10_in)) ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts )
  
            sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
            sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)
  
            sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
  
           !- - - - - - - - - - - -/ - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           ! uas [m s-1] i Eastward Near-Surface Wind
           ! vas [m s-1] i Northward Near-Surface Wind

          ELSE IF ((var_cmip(ivar) == "uas") .OR. (var_cmip(ivar) == "vas")) THEN
  
            IF (.not. ALLOCATED(u10_in)) ALLOCATE( u10_in ( xfocus, yfocus ), STAT=sts )
            IF (.not. ALLOCATED(v10_in)) ALLOCATE( v10_in ( xfocus, yfocus ), STAT=sts )
  
            sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
            sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)
  
            sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
              
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! tauu [Pa] Surface Downward Eastward Wind Stress
          ! tauv [Pa] Surface Downward Northward Wind Stress

          ELSE IF ( (var_cmip(ivar) == "tauu") .OR. (var_cmip(ivar) == "tauv") ) THEN
        
            IF (.not. ALLOCATED(cd_in)) ALLOCATE(cd_in ( xfocus, yfocus ), STAT=sts )
	    IF (.not. ALLOCATED(u10_in)) ALLOCATE(u10_in ( xfocus, yfocus ), STAT=sts ) 
	    IF (.not. ALLOCATED(v10_in)) ALLOCATE(v10_in ( xfocus, yfocus ), STAT=sts ) 
	    IF (.not. ALLOCATED(t2_in)) ALLOCATE(t2_in ( xfocus, yfocus ), STAT=sts ) 
	    IF (.not. ALLOCATED(psfc_in)) ALLOCATE(psfc_in ( xfocus, yfocus ), STAT=sts ) 
	    	
	    sts = NF90_INQ_VARID(ncidin, "CD", cd_varid)  
	    sts = NF90_GET_VAR(ncidin, cd_varid, cd_in(:,:), &
	      START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )  
	    sts = NF90_INQ_VARID(ncidin, "U10", u10_varid)
	    sts = NF90_GET_VAR(ncidin, u10_varid, u10_in(:,:), &
	      START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            sts = NF90_INQ_VARID(ncidin, "V10", v10_varid)  
            sts = NF90_GET_VAR(ncidin, v10_varid, v10_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )   		            
            sts = NF90_INQ_VARID(ncidin, "T2", t2_varid)  
            sts = NF90_GET_VAR(ncidin, t2_varid, t2_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )    
            sts = NF90_INQ_VARID(ncidin, "PSFC", psfc_varid)  
            sts = NF90_GET_VAR(ncidin, psfc_varid, psfc_in(:,:), &
              START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )       

          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! prhmax [kg m-2 s-1] m Daily Maximum Hourly Precipitation Rate (using wrfxtrm)
          ! reuse variables form rainc and rainnc
          ! special offset, as needed with minimum and maximum per file 
          ! as 'it' counts per input file this is OK

          ELSE IF ( (var_cmip(ivar) == "prhmax") .AND. (filetype(ivar) == "x" ) ) THEN
 
            IF (.not. ALLOCATED(rainc_max_in)) ALLOCATE( rainc_max_in ( xfocus, yfocus ), STAT=sts ) 
            IF (.not. ALLOCATED(rainnc_max_in)) ALLOCATE( rainnc_max_in ( xfocus, yfocus ), STAT=sts )
  
            sts = NF90_INQ_VARID(ncidin, "RAINCVMAX", rainc_varid)
            sts = NF90_INQ_VARID(ncidin, "RAINNCVMAX", rainnc_varid)
  
            sts = NF90_GET_VAR(ncidin, rainc_varid, rainc_max_in(:,:), &
              START = (/ xoffset, yoffset, it+1 /), COUNT = (/ xfocus, yfocus, 1 /) )
            sts = NF90_GET_VAR(ncidin, rainnc_varid, rainnc_max_in(:,:), &
              START = (/ xoffset, yoffset, it+1 /), COUNT = (/ xfocus, yfocus, 1 /) )

          ELSE

            PRINT *, "variable to read/write with no additional processing = ", var_wrf(ivar)  
            sts = NF90_INQ_VARID(ncidin, TRIM(var_wrf(ivar)), varid)
  
            IF ( ( cell_methods(ivar) == "minimum" ) .OR. &
                 ( cell_methods(ivar) == "maximum" ) ) THEN
              sts = NF90_GET_VAR(ncidin, varid, data_in(:,:), &
                    START = (/ xoffset, yoffset, it+1 /), COUNT = (/ xfocus, yfocus, 1 /) )
            ELSE 
              sts = NF90_GET_VAR(ncidin, varid, data_in(:,:), &
                    START = (/ xoffset, yoffset, it /), COUNT = (/ xfocus, yfocus, 1 /) )
            END IF

          END IF
  
          sts = NF90_CLOSE(ncidin)


!-------------------------------------------------------------------------------
!--------------------------------/-----------------------------------------------
! processing the data 

  
          PRINT *, "*** PROCESSING OF VARIABLES ***"

          IF ( (var_cmip(ivar) == "psl") .OR. &
               (var_cmip(ivar) == "prw") .OR. &
               (var_cmip(ivar) == "clwvi") .OR. &
               (var_cmip(ivar) == "clivi") .OR. &
               (var_cmip(ivar) == "clgvi") .OR. &
               (var_cmip(ivar) == "clhvi") .OR. &
               (var_cmip(ivar) == "cape") .OR. &
               (var_cmip(ivar) == "cin") .OR. &
               (var_cmip(ivar) == "li") .OR. &
               ( (plevel(ivar) /= -999) .AND. (filetype(ivar) == "s" ) ) ) THEN


            PRINT *, var_cmip(ivar),plevel(ivar)

            ! needs this initialisation, might be possible that for very low
            ! high levels calculation cannot be done, then it is set to official
            ! missing value
            data_in(:,:) = mv

            ! total pressure [Pa]
            ! perturbation pressure + base state pressure
            p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)

            ! PH and PHB are on nz+1 levels
            ! vertically destaggering
            DO nl = 1,nz 
              ph_fl(:,:,nl) = ((ph_in(:,:,nl)+phb_in(:,:,nl))+ &
                              (ph_in(:,:,nl+1)+phb_in(:,:,nl+1)))/2./gr
            END DO

            ! transfer theta-t0 to total potential temperature [K]
            ! and then convert potential temperature theta to absolute temperature
            t_in(:,:,:) = ( theta_in(:,:,:) + T00) * &
                          ( p_in(:,:,:) / P00 )**(R/cp)
                         
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! Set level on which the value is calculated to the plevel [hPa], 
            ! and convert it to [Pa]           
            pout(np)=plevel(ivar)*100.           
 
            ! vars needed: 3x int.: t_in, p_in, ph_fl 
            !              -> need: ph_in, phb_in, theta_in, pp_in, pb_in
            ! var_pl, var3d_in
            ! pout, slope, zg_pout
            ! comparison with wrfpress, same time: distribution OK, but ours is
            ! much ciolder, old vs new scheme about identical, ta500
            IF (INDEX(var_cmip(ivar),"ta") == 1) THEN
              PRINT *, "calc ta..."
              var3d_in(:,:,:) = t_in(:,:,:)
              ! linear in p
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN                    
                      slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                              (p_in(i,j,nl)-p_in(i,j,nl+1))
                      data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO

            ELSE IF (INDEX(var_cmip(ivar),"hus") == 1) THEN
              PRINT *, "interpolating hus..."
              var3d_in(:,:,:) = qv_in(:,:,:)
              !$OMP PARALLEL DO
              DO i = 1,xfocus
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN
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

            ELSE IF ( (INDEX(var_cmip(ivar),"ua") == 1 ) .AND. (filetype(ivar) == "s" ) ) THEN
              ! rotate to earth grid and destagger
              PRINT *, "calc ua..."
              DO i = 1,xfocus
                 var3d_in_u(i,:,:) = 0.5*(u_in(i,:,:)+u_in(i+1,:,:))     
              END DO
              DO j = 1,yfocus
                 var3d_in_v(:,j,:) = 0.5*(v_in(:,j,:)+v_in(:,j+1,:))
              END DO

              ! linear in p
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  var3d_in(i,j,:) = var3d_in_u(i,j,:)*cosalpha_in(i,j)- var3d_in_v(i,j,:)*sinalpha_in(i,j)
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN
                      slope = (var3d_in(i,j,nl)-var3d_in(i,j,nl+1))/&
                              (p_in(i,j,nl)-p_in(i,j,nl+1))
                      data_in(i,j) = var3d_in(i,j,nl+1) + &
                                    slope*(pout(np)-p_in(i,j,nl+1))
                    END IF
                  END DO
                END DO
              END DO
              !$OMP END PARALLEL DO
              
            ! comparison with wrfpress, same time: simple scheme about the same results
            ELSE IF ( (INDEX(var_cmip(ivar),"va") == 1) .AND. (filetype(ivar) == "s" ) ) THEN
              ! rotate to earth grid and destagger
              DO i = 1,xfocus
                 var3d_in_u(i,:,:) = 0.5*(u_in(i,:,:)+u_in(i+1,:,:)) 
              END DO
              DO j = 1,yfocus
                 var3d_in_v(:,j,:) = 0.5*(v_in(:,j,:)+v_in(:,j+1,:))  
              END DO

              ! linear in p
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  var3d_in(i,j,:) = var3d_in_u(i,j,:)*sinalpha_in(i,j) + var3d_in_v(i,j,:)*cosalpha_in(i,j)
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN                
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
            ELSE IF (INDEX(var_cmip(ivar),"wa") == 1) THEN
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

            ELSE IF ( INDEX(var_cmip(ivar),"zg") == 1 ) THEN
              PRINT *, "calc zg..."
              var3d_in(:,:,:) = ph_fl(:,:,:)
              ! either linear in p or log(p)
              !$OMP PARALLEL DO
              DO i = 1,xfocus 
                DO j = 1,yfocus
                  DO nl = 1,nz - 1
                    IF (pout(np).LE.p_in(i,j,nl) .AND. pout(np).GT.p_in(i,j,nl+1)) THEN               
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

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
                CALL calcslp(psl_in,p_in,qv_in,theta_in,ph_fl,nz,yfocus,xfocus,T00)
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
                CALL calcslptwo(psl_in, p_in, psfc_in, ph_in(:,:,1)+phb_in(:,:,1), t_in(:,:,1), yfocus, xfocus)

              CASE DEFAULT
                PRINT *, 'CAUTION: unknown setting for calcslp, proceed with default method'
                STOP 'calc_slp_type not properly set'
              END SELECT
              data_in(:,:) = psl_in(:,:)
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prw [kg m-2] i Water Vapor Path

            IF ( (var_cmip(ivar) == "prw") ) THEN  
              prw(:,:) = 0.
              DO nl = 1, nz
                prw(:,:) = prw(:,:) + &
                  (qv_in(:,:,nl) * p_in(:,:,nl)/(R*t_in(:,:,nl)) * &
                  ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+ &
                  phb_in(:,:,nl)))/gr)
              END DO
              data_in(:,:) = prw(:,:)
              WHERE (data_in < 0.) data_in = 0.
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! clwvi [kg m-2] i Condensed Water Path  

            IF ( (var_cmip(ivar) == "clwvi") ) THEN

              clwvi(:,:) = 0.
              DO nl = 1,nz - 1
                clwvi(:,:) = clwvi(:,:) + (qc_in(:,:,nl) + qi_in(:,:,nl) + &
                             qr_in(:,:,nl) + qs_in(:,:,nl) ) * p_in(:,:,nl)/ &
                             (R*t_in(:,:,nl)) * ((ph_in(:,:,nl+1)+ &
                             phb_in(:,:,nl+1)) - (ph_in(:,:,nl)+ &
                             phb_in(:,:,nl)))/gr
              END DO
              data_in(:,:) = clwvi(:,:)
              WHERE (data_in < 0.) data_in = 0.
            END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! clivi  [kg m-2] i Ice Water Path

            IF ( (var_cmip(ivar) == "clivi")) THEN
    
              !t_in(:,:,:) = (theta_in(:,:,:)+T00)*((pp_in(:,:,:)+pb_in(:,:,:))/P00)**(R/cp)
              !p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)
    
              clivi(:,:) = 0.
              DO nl = 1,nz - 1
                clivi(:,:) = clivi(:,:) + (qi_in(:,:,nl) + qs_in(:,:,nl)) * &
                             p_in(:,:,nl)/(R*t_in(:,:,nl)) * &
                             ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - &
                             (ph_in(:,:,nl)+phb_in(:,:,nl)))/gr
              END DO
              data_in(:,:) = clivi(:,:)
              WHERE (data_in < 0.) data_in = 0.
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! clgvi  [kg m-2] i Groupel Water Path

            IF ( (var_cmip(ivar) == "clgvi")) THEN

              !t_in(:,:,:) = (theta_in(:,:,:)+T00)*((pp_in(:,:,:)+pb_in(:,:,:))/P00)**(R/cp)
              !p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)

              clgvi(:,:) = 0.
              DO nl = 1,nz - 1
                clgvi(:,:) = clgvi(:,:) + (qg_in(:,:,nl)) * &
                             p_in(:,:,nl)/(R*t_in(:,:,nl)) * &
                             ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - &
                             (ph_in(:,:,nl)+phb_in(:,:,nl)))/gr
              END DO              
              data_in(:,:) = clgvi(:,:)
              WHERE (data_in < 0.) data_in = 0.
            END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! clhvi  [kg m-2] i Ice Water Path

            IF ( (var_cmip(ivar) == "clhvi")) THEN

              !t_in(:,:,:) = (theta_in(:,:,:)+T00)*((pp_in(:,:,:)+pb_in(:,:,:))/P00)**(R/cp)
              !p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)

              clhvi(:,:) = 0.
              DO nl = 1,nz - 1
                clhvi(:,:) = clhvi(:,:) + (qh_in(:,:,nl)) * &
                             p_in(:,:,nl)/(R*t_in(:,:,nl)) * &
                             ((ph_in(:,:,nl+1)+phb_in(:,:,nl+1)) - &
                             (ph_in(:,:,nl)+phb_in(:,:,nl)))/gr
              END DO
              data_in(:,:) = clhvi(:,:)
              WHERE (data_in < 0.) data_in = 0.
            END IF

  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! cape [J kg-1] i 2-D Maximum Convective Available Potential Energy
! cin [J kg-1] i 2-D Maximum Convective Inhibition

            IF ( (var_cmip(ivar) == "cape") .OR. (var_cmip(ivar) == "cin") .OR. &
                 (var_cmip(ivar) == "li") ) THEN        

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
                      t_p(i,j,nl+1) = (theta_in(i,j,1)+T00)*(p_in(i,j,nl+1)/P00)**(R/cp)
    
                    ELSE IF (qvs(i,j,nl) .lt. qv_in(i,j,1)) THEN ! moist adiabatic ascent    
                      IF (lcl(i,j) .eq. -999) THEN ! lifting condensation level
                        lcl(i,j) = p_in(i,j,nl)
                      END IF    
                      t_ii = t_p(i,j,nl)
                      t_tmp = 0.
                      DO ii = 1,100 ! solve iteratively 
                        qvs(i,j,nl+1) = 0.622*a*exp(b*(t_ii-c)/(t_ii-d))/p_in(i,j,nl+1)    
                        t_ii = t_ii - (t_ii*(P00/p_in(i,j,nl+1))**(R/cp)*exp(L*qvs(i,j,nl+1)/(cp*t_ii)) - &
                               (t_p(i,j,nl)*(P00/p_in(i,j,nl))**(R/cp)*exp(L*qvs(i,j,nl)/(cp*t_p(i,j,nl))))) / &
                               ( (P00/p_in(i,j,nl+1))**(R/cp)*exp(n/(p_in(i,j,nl+1)*t_ii)*exp(b*(t_ii-c)/(t_ii-d))) * &
                               (1 - (n/p_in(i,j,nl+1)*exp(b*(t_ii-c)/(t_ii-d))*(t_ii*(t_ii-b*c)+(b-2)*d*t_ii+d**2))/(t_ii*(d-t_ii)**2)) )  
                        IF ( ABS(t_ii - t_tmp) .le. 0.01 ) THEN
                          exit
                        ELSE
                          t_tmp = t_ii
                        END IF
                      END DO  
                      t_p(i,j,nl+1) = t_ii   
                    END IF                 
  
                    IF (t_p(i,j,nl) .gt. t_in(i,j,nl)) THEN                   
                      IF (lfc(i,j) .eq. -999) THEN ! level of free convection
                        lfc(i,j) = p_in(i,j,nl)
                      END IF
                    END IF
                  END DO
                 
                  DO nl = 1,nz-1
                    IF (lfc(i,j) .gt. 0.) THEN
                      IF ( (t_p(i,j,nl) .gt. t_in(i,j,nl)) .AND. (p_in(i,j,nl) .lt. lfc(i,j)) ) THEN   
                        cape(i,j) = cape(i,j) + (t_p(i,j,nl) - t_in(i,j,nl)) / t_in(i,j,nl) * ((phb_in(i,j,nl)+ph_in(i,j,nl))-(phb_in(i,j,nl-1)+ph_in(i,j,nl-1)))      
                      ELSE IF ( (t_p(i,j,nl) .lt. t_in(i,j,nl)) .AND. (p_in(i,j,nl) .ge. lfc(i,j)) ) THEN ! convective inhibition 
                        cin(i,j) = cin(i,j) + (t_in(i,j,nl) - t_p(i,j,nl)) / t_in(i,j,nl) * ((phb_in(i,j,nl)+ph_in(i,j,nl))-(phb_in(i,j,nl-1)+ph_in(i,j,nl-1)))
                      END IF
                    END IF
                    
                    IF ( (p_in(i,j,nl) .lt. 50000.) ) THEN
                      li(i,j) = t_in(i,j,nl) - t_p(i,j,nl)
                      exit
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
              IF ( var_cmip(ivar) == "li") THEN
                data_in(:,:) = li(:,:)
              END IF
             
            END IF
                  
             
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END IF ! pressure levels processing
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Processing variables on pressure levels        

          IF ( ( filetype(ivar) == "p" ) .AND. ( plevel(ivar) .NE. -999 ) ) THEN
            PRINT *, "Working on the pressure level", plevel(ivar)            
            ! Set levels            
            DO i = 1, num_levels
              IF (plevel(ivar) == pout(i)) THEN
                np = i
              EXIT
              END IF
            END DO

            !PRINT *, "extracting/calculating", var_cmip(ivar) 
            PRINT *, "New np level is", np
            !If wind components, derotate, first u component, then v component.
            !sinalpha and cosalha extracted fromt the geo_em file.
            IF ( (INDEX(var_cmip(ivar),"ua") == 1) ) THEN
              data_in(:,:) = pl_in_u(:,:,np)*cosalpha_in(:,:) - pl_in_v(:,:,np)*sinalpha_in(:,:)
            ELSE IF ( (INDEX(var_cmip(ivar),"va") == 1) ) THEN
              data_in(:,:) = pl_in_u(:,:,np)*sinalpha_in(:,:) + pl_in_v(:,:,np)*cosalpha_in(:,:)
            ELSE
              !for all other variables such as hus, ta and zg just read the datafrom the coresponing level. 
              data_in(:,:) = pl_in(:,:,np) 
            END IF

            WHERE (abs(data_in(:,:)) > 100000.) data_in(:,:) = mv
   
          END IF 
 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ua and va at height > 10m [m s-1] (e.g. ua100m, va100m) 
! The code implemented from CORDEX-WRF v1.3 module https://gmd.copernicus.org/articles/12/1029/2019/

          IF ( ( height(ivar) > 10. ) .AND. &
             ( (INDEX(var_cmip(ivar),"ua") == 1 ) .OR.  (INDEX(var_cmip(ivar),"va") == 1 ) ) ) THEN

            !calculating height
            DO nl = 1,nz_lowest 
               ph_fl(:,:,nl) = (((ph_in(:,:,nl)+phb_in(:,:,nl))+ &
                               (ph_in(:,:,nl+1)+phb_in(:,:,nl+1)))/2./gr) - hgt_in(:,:)
            END DO 
            
            !horizontal unstaggering     
            DO i = 1,xfocus 
               var3d_in_u(i,:,:) = 0.5*(u_in(i,:,:) + u_in(i+1,:,:))
            END DO
            
            DO j = 1,yfocus 
               var3d_in_v(:,j,:) = 0.5*(v_in(:,j,:) + v_in(:,j+1,:))
            END DO

            DO i = 1,xfocus 
              DO j = 1,yfocus
                CALL var_zwind( nz_lowest, var3d_in_u(i,j,:), var3d_in_v(i,j,:), ph_fl(i,j,:), u10_in(i,j), v10_in(i,j), &
                  sinalpha_in(i,j), cosalpha_in(i,j), real(height(ivar)), var2d_u(i,j), var2d_v(i,j) )
              END DO
            END DO 

            PRINT *, ph_fl(391,197,:), var3d_in_u(391,197,:), var3d_in_v(391,197,:)
            PRINT *, u10_in(391,197), v10_in(391,197), var2d_u(391,197), var2d_v(391,197)
              
            IF ( (INDEX(var_cmip(ivar),"ua") == 1 ) .AND. ( height(ivar) > 10 ) ) THEN
            	data_in(:,:) = var2d_u(:,:)
            ELSE IF ( (INDEX(var_cmip(ivar),"va") == 1 ) .AND. ( height(ivar) > 10 ) ) THEN
            	data_in(:,:) = var2d_v(:,:)
            END IF   
   
          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ta and huss at height > 10m [m s-1] (e.g. ta50m, hus50m) 

          IF ( ( height(ivar) > 10. ) .AND. &
             ( (INDEX(var_cmip(ivar),"ta") == 1 ) .OR.  (INDEX(var_cmip(ivar),"hus") == 1 ) ) ) THEN             
             
          
            IF ((INDEX(var_cmip(ivar),"ta") == 1)) THEN      
              p_in(:,:,:) = pp_in(:,:,:) + pb_in(:,:,:)                
              var3d_in(:,:,:) = ( theta_in(:,:,:) + T00) * ( p_in(:,:,:) / P00 )**(R/cp)
            END IF

            !calculating height
            DO nl = 1,nz_lowest 
               ph_fl(:,:,nl) = (((ph_in(:,:,nl)+phb_in(:,:,nl))+ &
                               (ph_in(:,:,nl+1)+phb_in(:,:,nl+1)))/2./gr) - hgt_in(:,:)
            END DO 
      

            DO i = 1,xfocus 
              DO j = 1,yfocus
                CALL linear_int(nz_lowest, var3d_in(i,j,:), ph_fl(i,j,:), var2d_in(i,j), real(height(ivar)), var2d_out(i,j))
              END DO
            END DO 
            
            data_in(:,:) = var2d_out(:,:)
   
          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface downward wind stress
! Following Trenberth et al.(1990; doi: https://doi.org/10.1175/1520-0485(1990)020%3C1742:TMACIG%3E2.0.CO;2)

          IF ( (var_cmip(ivar) == "tauu") .OR. (var_cmip(ivar) == "tauv") ) THEN

            IF (var_cmip(ivar) == "tauu") THEN
  
              data_in(:,:) = cd_in*(psfc_in/(R*t2_in)) * &
                             ( (u10_in(:,:)**2 + v10_in(:,:)**2)**0.5 ) * &
                             ( u10_in(:,:)*cosalpha_in(:,:) - v10_in(:,:)*sinalpha_in(:,:) )    

            ELSE IF (var_cmip(ivar) == "tauv")  THEN   

              data_in(:,:) = cd_in*(psfc_in/(R*t2_in)) * &
                             ( (u10_in(:,:)**2 + v10_in(:,:)**2)**0.5 ) * &
                             ( u10_in(:,:)*sinalpha_in(:,:) + v10_in(:,:)*cosalpha_in(:,:) )

            END IF

          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! cll, clm, clh [%] a Total Cloud Fraction
! not instantaneous, but average, have to recode here

          IF ( (var_cmip(ivar) == "cll") .OR. &
               (var_cmip(ivar) == "clm") .OR. &
               (var_cmip(ivar) == "clh") .OR. &
               (var_cmip(ivar) == "clt") ) THEN
            
            p_in = ( pp_in(:,:,:) + pb_in(:,:,:) ) / 100.

            SELECT CASE (var_cmip(ivar))
            CASE('clh')
              low_lev=440.
              high_lev=minval(p_in(:,:,:))
            CASE('clm')
              low_lev=680.
              high_lev=440.
            CASE('cll')
              low_lev=maxval(p_in(:,:,:))
              high_lev=680.
            CASE('clt')
              low_lev=maxval(p_in(:,:,:))
              high_lev=minval(p_in(:,:,:))
            CASE DEFAULT
              PRINT *, "Invalid time interval specified"
              STOP
            END SELECT

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
                        IF ( p_in(i,j,nl) < low_lev .AND. p_in(i,j,nl) >= high_lev) THEN
                          IF (( i == 150 ) .AND. ( j == 150)) THEN
                            PRINT *, p_in(i,j,nl)
                          END IF
                          cldfra_inv(i,j) = cldfra_inv(i,j) * &
                            (1- max(cldfra_in(i,j,nl,k),cldfra_in(i,j,nl-1,k)) / &
                            (1-cldfra_in(i,j,nl-1,k)))
                        END IF
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
! pr [kg m-2 s-1] Precipitation; [mm dtHours-1] -> [kg m-2 s-1]
! accumulated quantity (using a bucket)

          IF ( (var_cmip(ivar) == "pr") .OR. &
             ( (var_cmip(ivar) == "prhmax") .AND. (filetype(ivar) == "s" ) ) ) THEN
            IF (calc) THEN
  
              IF ( bucket_mm > 0. ) THEN 
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
            ELSE
              data_in(:,:) = mv
            END IF
            calc = .TRUE.
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! prsn [kg m-2 s-1] ? Snowfall Flux or snm [kg m-2 s-1] a Surface Snow Melt
! unit [mm/3hr] to [kg m-2 s-1]
  
          IF ( (var_cmip(ivar) == "prsn") .OR. (var_cmip(ivar) == "snm") ) THEN  
            IF (calc) THEN            
              data_in(:,:) = (acsnom_in(:,:,2) - acsnom_in(:,:,1)) / &
                             (dtHours*3600.)  
            ELSE
              data_in(:,:) = mv
            END IF
            WHERE (landmask_in == 0) data_in = 0.
            calc = .TRUE.
         END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! evspsbl [kg m-2 s-1] a Evaporation
! unit [kg m-2 /3hr] to [kg m-2 s-1]
! TO CHECK STILL, in principle OK
  
          IF (var_cmip(ivar) == "evspsbl") THEN
              IF (calc) THEN            
              data_in(:,:) = ( sfcevp_in(:,:,1) + sfcevp_in(:,:,2) ) / 2.  
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
              data_in(:,:) = (sfroff_in(:,:,2) + sfroff_in(:,:,1)) / &
                             (2.*dtHours) !(dtHours*3600)
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
              data_in(:,:) = ( ( sfroff_in(:,:,2) + sfroff_in(:,:,1) )/2. + &
                               ( udroff_in(:,:,2) + udroff_in(:,:,1) )/2. ) / &
                               (dtHours) !(dtHours*3600) 
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
! rsdsdir [W m-2] a Surface Direct Downwelling Shortwave Radiation
! rsdscs [W m-2] a Surface Downwelling Clear-Sky Shortwave Radiation
! rldscs [W m-2] a Surface Downwelling Clear-Sky Longwave Radiation
! rsuscs [W m-2] a Surface Upwelling Clear-Sky Shortwave Radiation
! rluscs [W m-2] a Surface Upwelling Clear-Sky Longwave Radiation
! rsutcs [W m-2] a TOA Outgoing Clear-Sky Shortwave Radiation
! rlutcs [W m-2] a TOA Outgoing Clear-Sky Longwave Radiation
  
          IF ( (var_cmip(ivar) == "rsds") &
            .OR. (var_cmip(ivar) == "rlds") &
            .OR. (var_cmip(ivar) == "rsus") &
            .OR. (var_cmip(ivar) == "rlus") &
            .OR. (var_cmip(ivar) == "rlut") &
            .OR. (var_cmip(ivar) == "rsdt") &
            .OR. (var_cmip(ivar) == "rsut") &
            .OR. (var_cmip(ivar) == "rsdsdir") &
            .OR. (var_cmip(ivar) == "rsdscs") &
            .OR. (var_cmip(ivar) == "rldscs") &
            .OR. (var_cmip(ivar) == "rsuscs") &
            .OR. (var_cmip(ivar) == "rluscs") &
            .OR. (var_cmip(ivar) == "rsutcs") &
            .OR. (var_cmip(ivar) == "rlutcs") ) THEN
    
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
            data_in(:,:) =  (smois_in(:,:,1)*DZS(1) + &
                             smois_in(:,:,2)*DZS(2) + &
                             smois_in(:,:,3)*DZS(3) + &
                             smois_in(:,:,4)*DZS(4) ) * 1000.
            WHERE (landmask_in == 0) data_in = mv
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrsos [kg m-2] Soil Moisture Content at the 1st model layer
! see comment in the variable aquisition section
! m3 m-3 -> kg m-2
! SMOIS_i [m3 m-3] * DZS_1 [m] * 1000 [kg m-3] = [kg m-2]
! assume standard area 1m2, assume 1l water = 1kg, 1mm m-2 water depth = 1l
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrsos") THEN    
            data_in(:,:) =  (smois_in(:,:,1)*DZS(1)) * 1000.
            WHERE (landmask_in == 0) data_in = mv
          END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrsol
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrsol") THEN
	    DO i = 1, 4
	      data_in_3D(:,:,i) = smois_in(:,:,i) * DZS(i) * 1000.
              WHERE (landmask_in == 0) data_in_3D(:,:,i) = mv 
            END DO
          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrfso [kg m-2] Total soil Frozen Water Content
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrfso") THEN   
            data_in(:,:) =  ( ( smois_in(:,:,1) - sh2o_in(:,:,1) )*DZS(1) + &
                              ( smois_in(:,:,2) - sh2o_in(:,:,2) )*DZS(2) + &
                              ( smois_in(:,:,3) - sh2o_in(:,:,3) )*DZS(3) + &
                              ( smois_in(:,:,4) - sh2o_in(:,:,4) )*DZS(4) ) * 1000.
            WHERE (landmask_in == 0) data_in = mv
          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrfsos [kg m-2] Soil Frozen Water Content at the 1st model layer
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrfsos") THEN    
            data_in(:,:) =  ( (smois_in(:,:,1) - sh2o_in(:,:,1))*DZS(1)) * 1000.
            WHERE (landmask_in == 0) data_in = mv
          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mrsfl [kg m-2] Soil Frozen Water Content per layer
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "mrsfl") THEN
	    DO i = 1, 4
	      data_in_3D(:,:,i) = (smois_in(:,:,i) - sh2o_in(:,:,i)) * DZS(i) * 1000.
              WHERE (landmask_in == 0) data_in_3D(:,:,i) = mv 
            END DO
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! tsl [K] soil temperature
! masking over ocean, with mv initialized

          IF (var_cmip(ivar) == "tsl") THEN
	    DO i = 1, 4
	      data_in_3D(:,:,i) = tslb_in(:,:,i)
	      WHERE (landmask_in == 0) data_in_3D(:,:,i) = mv
            END DO
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! snc [%] i Snow Area Fraction
! unit [] to [%]
! masking over ocean, with mv initialized

          IF ( (var_cmip(ivar) == "snc") ) THEN
            data_in(:,:) = data_in(:,:)*100. 
            WHERE (landmask_in == 0) data_in = mv
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! snd [m] i Snow Depth
! snw [kg m-2] i Surface Snow Amount
! masking over ocean, with mv initialized
! variables are just passed through and need some masking

          IF ( (var_cmip(ivar) == "snd") &
            .OR. (var_cmip(ivar) == "snw") ) THEN 
            WHERE (landmask_in == 0) data_in = mv  
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sic [%] ?>i Sea Ice Area Fraction
! no temporal aggregation defined, treat as tier-2 instantaneous data, most
! reasonable; also offers the possibility of some proper sea-ice treatment

          IF ( (var_cmip(ivar) == "sic") ) THEN
              data_in(:,:) = data_in(:,:) * 100.   
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sftlaf [%] where LU_INDEX=21 set to 100, otherwise set to 0.

          IF ( (var_cmip(ivar) == "sftlaf") ) THEN
              data_in(:,:) = data_in(:,:) * 100.
          END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sfturf [%] Urban fraction
          IF ( (var_cmip(ivar) == "sfturf") ) THEN
              WHERE (data_in == 13 .OR. data_in > 30 ) data_in = 100.
              WHERE (data_in < 100.) data_in = 0.
          END IF
          
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sftlf [%] land-sea mask
          IF ( (var_cmip(ivar) == "sftlf") ) THEN
              WHERE (data_in == 2 ) data_in = 0.
              data_in(:,:) = data_in(:,:) * 100.
          END IF
  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! [ms-1] i Near-Surface Wind Speed
  
          IF ( (var_cmip(ivar) == "sfcWind") .OR. &
             ( (var_cmip(ivar) == "sfcWindmax") .AND. (filetype(ivar) == "s" ) ) )THEN  
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

          IF ( (var_cmip(ivar) == "prhmax") .AND. (filetype(ivar) == "x" ) ) THEN
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
            ALLOCATE(tmp_2d(xfocus,yfocus),tmp1_2d(xfocus,yfocus),tmp2_2d(xfocus,yfocus),tmp3_2d(xfocus,yfocus),tmp4_2d(xfocus,yfocus))
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
          END IF

!-------------------------------------------------------------------------------
! write data to netCDF file
  
          sts = NF90_OPEN( TRIM(DirOutputPostProRoot) // "/" // TRIM(pn_out) &
            // "/" // TRIM(fn_out), NF90_WRITE, ncid )

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

          IF ( (var_cmip(ivar) == "mrsol") &
             .OR. (var_cmip(ivar) == "mrsfl") &
             .OR. (var_cmip(ivar) == "tsl")) THEN
              sts = NF90_PUT_VAR( ncid, x_varid, data_in_3D(:,:,:), &
                START=(/ 1, 1, 1, counter /), COUNT = (/ xfocus, yfocus, 4, 1 /) )
          ELSE IF (frequency(ifrq) == 'fx')  THEN
            sts = NF90_PUT_VAR( ncid, x_varid, data_in(:,:),  &
              START=(/ 1, 1 /), COUNT = (/ xfocus, yfocus /) )
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
          PRINT *, "NF90_CLOSE", sts

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

#ifndef SERIAL
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF ( ierr /= MPI_SUCCESS ) STOP "Problem with MPI_BARRIER"
    CALL mpi_finalize(ierr)
    IF ( ierr /= MPI_SUCCESS ) stop "Problem with MPI_FINALIZE"
    DEALLOCATE( ivar_list )
#else 
    END DO ! ivarnml - namelist loop     
#endif
END DO ! ifrq - different temporal aggregations

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

SUBROUTINE calcslptwo(slp, PP, P_s, PHI_s, T_L, ns, ew)

IMPLICIT NONE

real, DIMENSION(:,:),   INTENT(OUT) :: slp        !(output) sea level pressure
real, DIMENSION(:,:,:), INTENT(IN)  :: PP         !(input) 3D pressure
real, DIMENSION(:,:),   INTENT(IN)  :: P_s        !(input) pressure at surface
real, DIMENSION(:,:),   INTENT(IN)  :: PHI_s      !(input) 2D geopotential of the surface
real, DIMENSION(:,:),   INTENT(IN)  :: T_L        !(input) temperature at lowest level
INTEGER,                INTENT(IN)  :: ns, ew     !(input) dimensions: north-south, east-west

real, DIMENSION(ew,ns)              :: P_L        !(calculated) pressure at lowest level
real                                :: T_surf     !(calculated) surface temperature
real                                :: gamma_mod  !(calculated) modified lapse rate that is actually used for calculations
real                                :: x          !(calculated) expansion coefficient of (eq. 9)
real                                :: T_0        !(calculated) auxiliary variable

real, PARAMETER                     :: Rd    = 287.04 ![J kg-1 K-1] (const.) dry air constant
real, PARAMETER                     :: g     = 9.81   ![m s-2] (const.) acceleration due to gravity
real, PARAMETER                     :: gamma = 0.0065 ![K m-1] (const.) lapse rate at const. 0.0065 K/m, also denoted as (dT/dz)_st
integer                             :: j,k            !loop parameters

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

!*********************************************************************************************
! Josipa: Add subroutine to extrapolate the wind at a given height following 
! the 'power law' methodology, adapted from the CORDEX-WRF v1.3 (Fita et al. 2018; 
! https://doi.org/10.5194/gmd-12-1029-2019)
!    wss[newz] = wss[z1]*(newz/z1)**alpha
!    alpha = (ln(wss[z2])-ln(wss[z1]))/(ln(z2)-ln(z1))
! Original source is Phd Thesis: Benedicte Jourdier. Ressource eolienne en France metropolitaine: 
! methodes d’evaluation du potentiel, variabilite et tendances. Climatologie. 
! Ecole Doctorale Polytechnique, 2015. French
!*********************************************************************************************
SUBROUTINE var_zwind(nz, u, v, z, u10, v10, sa, ca, newz, unewz, vnewz)

IMPLICIT NONE

INTEGER, INTENT(in)                 :: nz
REAL, DIMENSION(nz), INTENT(in)     :: u,v,z
REAL, INTENT(in)                    :: u10, v10, sa, ca, newz
REAL, INTENT(out)                   :: unewz, vnewz 

! Local
INTEGER                             :: inear
REAL                                :: zaground
REAL, DIMENSION(2)                  :: v1, v2, zz, alpha, uvnewz
REAL                                :: min_value=1.0e-7


!!!!!!! Variables
! nz: number of vertical levels
! u,v: vertical wind components [ms-1]
! z: height above surface [m]
! u10,v10: 10-m wind components [ms-1]
! sa, ca: local sine and cosine of map rotation [1.]
! newz: height to which extrapolate
! unewz,vnewz: Wind compoonents at the given height [ms-1]

! Looking for the level  below desired height
IF (z(1) < newz ) THEN
  DO inear = 1,nz-2
    ! L. Fita, CIMA. Feb. 2018
    !! Choose between extra/inter-polate. Maybe better interpolate?
    ! Here we extrapolate from two closest lower levels
    ! zaground = z(inear+2)
    ! Here we interpolate between levels
    zaground = z(inear+1)
    IF ( zaground >= newz) EXIT
  END DO
ELSE
  inear = nz - 2
END IF

IF (inear == nz-2) THEN
! No vertical pair of levels is below newz, using 10m wind as first value
! and the first level as the second
   v1(1) = sign(MAX(ABS(u10),min_value),u10) 
   v1(2) = sign(MAX(ABS(v10),min_value),v10)  
   v2(1) = sign(MAX(ABS(u(1)),min_value),u(1)) 
   v2(2) = sign(MAX(ABS(v(1)),min_value),v(1)) 
   zz(1) = 10.
   zz(2) = z(1)
ELSE
   v1(1) = sign(MAX(ABS(u(inear)),min_value),u(inear))
   v1(2) = sign(MAX(ABS(v(inear)),min_value),v(inear))
   v2(1) = sign(MAX(ABS(u(inear+1)),min_value),u(inear+1))
   v2(2) = sign(MAX(ABS(v(inear+1)),min_value),v(inear+1))
   zz(1) = z(inear)
   zz(2) = z(inear+1)
END IF

! Computing for each component
alpha = (LOG(ABS(v2))-LOG(ABS(v1)))/(LOG(zz(2))-LOG(zz(1)))
uvnewz = v1*(newz/zz(1))**alpha

! Earth-rotation
unewz = uvnewz(1)*ca - uvnewz(2)*sa
vnewz = uvnewz(1)*sa + uvnewz(2)*ca

RETURN

END SUBROUTINE var_zwind

!===============================================================================

SUBROUTINE linear_int(nz, var3d, z, var2d, newz, varout)
!!!!!!! Variables
! nz: number of vertical levels
! var3d: 3d variable
! z: height above surface [m]
! var2d: 2d variable (near-surface variable)
! newz: height to which extrapolate
! unewz,vnewz: Wind compoonents at the given height [ms-1]

IMPLICIT NONE

INTEGER, INTENT(in)                 :: nz
REAL, DIMENSION(nz), INTENT(in)     :: var3d, z
REAL, INTENT(in)                    :: var2d, newz
REAL, INTENT(out)                   :: varout 

! Local
INTEGER                             :: inear
REAL                                :: zaground, slope
REAL, DIMENSION(2)                  :: var, zz

IF (z(1) < newz ) THEN
  DO inear = 1,nz-2
    zaground = z(inear+1)
    IF ( zaground >= newz) EXIT
  END DO
ELSE
  inear = nz - 2
END IF

IF (inear == nz-2) THEN
! No vertical pair of levels is below newz, using 2m variable as first value
! and the first level as the second
   var(1) = var2d 
   var(2) = var3d(1)  
   zz(1)  = 2.
   zz(2)  = z(1)
ELSE
   var(1) = var3d(inear)
   var(2) = var3d(inear+1)
   zz(1) = z(inear)
   zz(2) = z(inear+1)
END IF

slope = (var(2)-var(1))/(zz(2)-zz(1))
varout = var(1) + (newz - zz(1))*slope

RETURN

END SUBROUTINE linear_int

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
CASE ('fx') 
  dtDecDay = 1.0 / 24.0_8
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
  ELSE IF ( CheckForLeapyear( i ) == 360 ) THEN 
    ndpm = (/30,30,30,30,30,30,30,30,30,30,30,30/)
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

    IF ( TRIM(calendar) == "noleap" ) THEN
      CheckForLeapyear = 365
    ELSE IF ( TRIM(calendar) == "360_days" ) THEN
      CheckForLeapyear = 360
    ELSE
      IF ( ( MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0 ) .OR. ( MOD(year, 400) == 0 ) ) THEN
        CheckForLeapyear = 366
      ELSE
        CheckForLeapyear = 365
      END IF
    END IF

  END FUNCTION CheckForLeapyear

!-------------------------------------------------------------------------------

END SUBROUTINE CreateRefTimeArray
