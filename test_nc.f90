PROGRAM test_nc

! source loadenv.JURECA-DC_2020_Intel-PSMPI.ini
! get_unpack_compile_cdi_lib.sh
! make -f Makefile_test_nc clean
! make -f Makefile_test_nc
! export LD_LIBRARY_PATH=/p/scratch/cjjsc39/jjsc3900/sim/tmp2/tools/pCMORizer.f90/extLib/cdi210/lib:$LD_LIBRARY_PATH

USE netcdf

IMPLICIT NONE

INCLUDE 'cdi.inc'

INTEGER :: ncid_in, varid
INTEGER :: sts
REAL, DIMENSION(1600,1552) :: data

INTEGER :: i
INTEGER :: streamID, vlistID, taxisID, tsID, vdate, vtime, nts
PARAMETER (nts = 5)

sts = NF90_OPEN('2020/lffd2020011320.nc', NF90_NOWRITE, ncid_in)
sts = NF90_INQ_VARID(ncid_in, "T_2M", varid)
sts = NF90_GET_VAR(ncid_in, varid, data(:,:))
sts = NF90_CLOSE(ncid_in)

print *, SHAPE(data)
print *, data(500:510,500)

! does only work with classic netCDF files
! but there is a check for 
streamID = streamOpenRead('2020/lffd2020011320.nc')
vlistID = streamInqVlist(streamID)
taxisID = vlistInqTaxis(vlistID)
DO tsID = 0, nts-1 
  sts = streamInqTimestep(streamID, tsID) 
  vdate = taxisInqVdate(taxisID) 
  vtime = taxisInqVtime(taxisID) 
  WRITE(0, *)"read timestep:",sts,tsID+1,"date=",vdate,"time=",vtime 
END DO
call streamClose(streamID)

END PROGRAM test_nc
