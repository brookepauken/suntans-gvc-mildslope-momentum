#!/bin/sh
export PATH=/opt/homebrew/bin:$PATH

SUNTANSHOME=../../main/
MATLABPATH=/Applications/MATLAB_R2024b.app/bin/matlab
datadir=data_bl_test


#set number of grid points in x and directory names 
Nx=400
sh_val=2.5 #sets height of pycnocline realtive to surface where sh_val/10*H = depth of pycnocline
nudge=0
NUMPROCS=4

echo running diffusion test case

#generate grid in matlab
$MATLABPATH -nosplash -nodesktop -r "Nx = $Nx; run('./mfiles/quadgrid_periodic.m'); exit;"
echo Generated grid in matlab

#create data directory and copy over rundata and sources.c so we can see later what nudging we were doing
if [ -d $datadir ]; then
   rm -r $datadir
fi
    
cp -r rundata/. $datadir
if [ ! -d  $datadir ] ; then
   cp -r rundata/. $datadir
fi

echo copied over rundata

#generate grid through suntans 
mpirun -np $NUMPROCS $SUNTANSHOME/sun -g -vv --datadir=$datadir
echo generated grid through suntans

#generate wave with DJLES in matlab 
$MATLABPATH -nosplash -nodesktop -r "datadir = '$datadir'; nudge = $nudge; nx = $Nx; sh_val = $sh_val; Atarget = 1.5E8; run('./mfiles/initial_DJL_sherlock_midstretch.m'); exit;"

echo created wave in matlab

#run suntans 
mpirun -np $NUMPROCS $SUNTANSHOME/sun -g -s -vv --datadir=$datadir
echo ran suntans