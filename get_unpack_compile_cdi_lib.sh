#!/bin/bash

source loadenv.JURECA-DC_2022_GCC-OpenMPI.ini
#wget https://code.mpimet.mpg.de/attachments/download/27452/cdi-2.1.0.tar.gz
#tar -xvzf cdi-2.1.0.tar.gz
InstallPath=$(pwd)/extLib/cdi210
cd cdi-2.1.0
make clean
make distclean
./configure --prefix=$InstallPath --with-netcdf=$EBROOTNETCDFMINFORTRAN
make
make check
make install

# risk is serialism
# export LD_LIBRARY_PATH=/p/project/cjjsc39/jjsc3900/test_cdi_fortran/test_gcc_latest_3/lib:$LD_LIBRARY_PATH
# gfortran -I/p/project/cjjsc39/jjsc3900/test_cdi_fortran/test_gcc_latest_3/include -I$EBROOTNETCDFMINFORTRAN/include test_cdi_cos.f90 -L/p/project/cjjsc39/jjsc3900/test_cdi_fortran/test_gcc_latest_3/lib -lcdi -L$EBROOTNETCDFMINFORTRAN/lib -lnetcdf
