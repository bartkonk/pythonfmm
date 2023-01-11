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

CURRDIR=$( pwd )
INDIR=../02-energyMinimization

for SIZE in 100 ; do
    gmx grompp -f equilibration.mdp -c $INDIR/minimized_$SIZE.gro -p $INDIR/waters_$SIZE.top -o equil_$SIZE.tpr
    func.testquit $?
    
    gmx_threads_d_AVX2_256 mdrun -s equil_$SIZE.tpr -deffnm equilibrated_$SIZE -v
    func.testquit $?
done

