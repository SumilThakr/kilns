#!/bin/bash -l

#SBATCH --job-name run01
#SBATCH -N 1-1
#SBATCH --time=72:00:00
#SBATCH --ntasks=96
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=sthakrar@umn.edu
#=================================================

#Run number
RN=01

cd sthakrar/runbrick/

# Remove old .o file
rm -rf runFA-$RN*o*

#copy input file
rm ./input.geos
cp ./input.dir/input.geos.$RN ./input.geos

ulimit
date
time

time ./gcclassic

# show when run finished
date

exit 0
