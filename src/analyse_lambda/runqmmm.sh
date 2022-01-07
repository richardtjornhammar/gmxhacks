#!/bin/bash
module add easy
module add heimdal
module add i-compilers
module add mpi
module add gaussian

source /afs/nada.kth.se/home/1/u1exw3x1/PDC/projdir/prg/grompqc/bin/GMXRC.bash
export DEVEL_DIR=/afs/nada.kth.se/home/1/u1exw3x1/PDC/projdir/modlinks
export GAUSS_DIR=/pdc/vol/gaussian/G03RevD.01/g03
export GAUSS_EXE=g03
export GAUSS_SCRDIR=.
export QM_EXE=$GAUSS_EXE
export GAUSS_CPUS=8

PRG="$1"
shift
ARGS="$*"

echo $PWD
mpirun --n 1 /afs/nada.kth.se/home/1/u1exw3x1/PDC/projdir/prg/grompqc/bin/$PRG $ARGS
