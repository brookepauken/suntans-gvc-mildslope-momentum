#!/bin/sh
export PATH=/opt/homebrew/bin:$PATH

SUNTANSHOME=../../main/
MATLABPATH=/Applications/MATLAB_R2024b.app/bin/matlab
numprocs=4
datadir=data_turbtest

echo Starting turbulent channel test case.

# $MATLABPATH -nosplash -nodesktop -r "run('./mfiles/quadgrid_periodic.m'); exit;"
# echo generated grid files in MATLAB!

#create data directory and copy over rundata
cp -r rundata/. $datadir
if [ ! -d  $datadir ] ; then
   cp -r rundata/. $datadir
fi
echo copied over rundata

#generate grid                                                                                                                             
mpirun -np $numprocs $SUNTANSHOME/sun -g -vv --datadir=./$datadir
echo generated grid!

#run model
mpirun -np $numprocs $SUNTANSHOME/sun -s -vv --datadir=./$datadir
echo ran test case!
