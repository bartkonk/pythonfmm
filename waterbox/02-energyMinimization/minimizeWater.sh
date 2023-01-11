#!/bin/bash
source /usr/local/gromacs/GMXRC2022
source /usr/local/cuda-11.1/CUDA111RC

shopt -s -o nounset

function func.testquit
{
    if [ "$1" = "0" ] ; then
        echo "OK"
    else
        echo "ERROR: exit code of the last command was $1. Exiting."
        exit
    fi
}

# Waters 100, 1,000 10,000 100,000 1,000,000 10,000,000  33,333,333
# Atoms  300, 3,000 30,000 300,000 3,000,000 30,000,000 999,999,999

CURRDIR=$( pwd )
INDIR=../01-buildBoxes

# Use pdb2gmx to produce a topology for the waterbox system
SYSTEM=waters_100
gmx pdb2gmx -f $INDIR/$SYSTEM.pdb -p $SYSTEM.top -water tip3p -ff amber03
func.testquit $?

# Topologies for the larger systems can be build from the small topology
#sed 's/SOL              1000/SOL             10000/' waters_1000.top > waters_10k.top
#sed 's/SOL              1000/SOL            100000/' waters_1000.top > waters_100k.top
#sed 's/SOL              1000/SOL           1000000/' waters_1000.top > waters_1M.top
#sed 's/SOL              1000/SOL          10000000/' waters_1000.top > waters_10M.top
#sed 's/SOL              1000/SOL          33333333/' waters_1000.top > waters_33M.top

for SIZE in 100 ; do
    gmx grompp -f minimize.mdp -c $INDIR/waters_$SIZE.pdb -p waters_$SIZE.top -o equil_$SIZE.tpr
    func.testquit $?
    
    gmx_threads_AVX2_256 mdrun -s equil_$SIZE.tpr -deffnm minimized_$SIZE
    func.testquit $?
done

