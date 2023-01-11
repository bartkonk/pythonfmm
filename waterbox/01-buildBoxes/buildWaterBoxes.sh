#!/bin/bash
source /usr/local/gromacs/GMXRC2022

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


# Waters 100, 1,000 10,000 100,000 1,000,000 10,000,000
# Atoms  300, 3,000 30,000 300,000 3,000,000 30,000,000

# 100 W occupy 100 * 30 *  A^3 = 3000 A^3 = 1.47 nm ^3 (roughly)
LENGTH=1.47
gmx solvate -cs spc216.gro -o waters_100.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 100

# 1,000 W occupy 1000 * 30 *  A^3 = 30000 A^3 = 3.107 nm ^3
LENGTH=3.152
gmx solvate -cs spc216.gro -o waters_1000.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 1000

# 10,000 W occupy 10,000 * 30 *  A^3 = 300000 A^3 = 6.7 nm ^3
LENGTH=6.74
gmx solvate -cs spc216.gro -o waters_10k.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 10000

# 100k waters occupy 100,000 * 30 *  A^3 = 3,000,000 A^3 = 14.4 nm ^3
LENGTH=14.51
gmx solvate -cs spc216.gro -o waters_100k.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 100000

# 1 M waters occupy 1,000,000 * 30 *  A^3 = 30,000,000 A^3 = 31.1 nm ^3
LENGTH=31.25
gmx solvate -cs spc216.gro -o waters_1M.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 1000000

# 10 M waters occupy 10,000,000 * 30 *  A^3 = 300,000,000 A^3 = 66.94 nm ^3
LENGTH=67.5
gmx solvate -cs spc216.gro -o waters_10M.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 10000000

# 33,333,333 waters occupy 33,333,333 * 30 *  A^3 = 1,000,000,000 A^3 = 100 nm ^3
LENGTH=101
gmx solvate -cs spc216.gro -o waters_10M.pdb -box $LENGTH $LENGTH $LENGTH -maxsol 33333333


