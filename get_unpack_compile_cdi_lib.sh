#!/bin/bash

#wget https://code.mpimet.mpg.de/attachments/download/27452/cdi-2.1.0.tar.gz
#tar -xvzf cdi-2.1.0.tar.gz
InstallPath=$(pwd)/extLib/cdi210

#source loadenv.JURECA-DC_2022_GCC-OpenMPI.ini
source loadenv.JURECA-DC_2020_Intel-PSMPI.ini

cd cdi-2.1.0
make clean
make distclean
#./configure --prefix=$InstallPath --with-netcdf=$EBROOTNETCDFMINFORTRAN
./configure --prefix=$InstallPath --with-netcdf=$EBROOTNETCDFMINFORTRAN FC=mpifort
make
make check
make install

# export LD_LIBRARY_PATH=/p/scratch/cjjsc39/jjsc3900/sim/tmp2/tools/pCMORizer.f90/extLib/cdi210/lib:$LD_LIBRARY_PATH
