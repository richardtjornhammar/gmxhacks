  #!/bin/bash
  #PBS -A SNIC022-09-18
  #PBS -o job.out
  #PBS -e job.err
  #PBS -m abe
  #PBS -l nodes=3:ppn=8
  #PBS -l walltime=82:00:00
  #PBS -N QMMM_zinc_6h2O_L_l0.48

  # RUN IN CORRECT DIR
  cd $PBS_O_WORKDIR

  # MAKE SURE BAD MODULES ARE NOT LOADED
  module rm gromacs
  module rm psc
  module rm openmpi/3.3.3/psc

  # LOAD CORRECT MODULES
  module add intel-compiler
  module add openmpi/1.4/intel
  module add libfftw/3
  module add psc
  module add mpqc

  # SETTING UP THE ENVIRONMENT
  source /home/r/richardt/pfs/prg/grompqc/bin/GMXRC.bash
  export QM_EXE=/lap/mpqc/2.3.1/bin/mpqcrun
  export GMX_EXE=/home/r/richardt/pfs/prg/grompqc/bin/mdrun 
  export MPQC_CPUS=8
  export MPQC_NODES=3
  export OPTIONS="--nodefile=$PBS_NODEFILE --messagegrp=mpi --memorygrp=mtmpi"
  export LAUNCH="mpirun [--hostfile %NODEFILE%] -n %NPROC% %MPQCRUNPROC% [-o %OUTPUT%] %INPUT%"

  # GROMACS DEFAULTS
  RUN=wrk-qm
  MAXH=81

  $GMX_EXE -s $RUN.tpr -maxh $MAXH

