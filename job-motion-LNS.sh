#!/bin/bash
#SBATCH -J motion 
#SBATCH -o motion.o%j
#SBATCH -n 96
#SBATCH --ntasks-per-node=24
#SBATCH -p comp
#SBATCH -t 120:00:00

echo $SLURM_JOB_ID
echo $SLURM_PRUEBA
echo $SLURM_JOB_NUM_NODES

source /opt/rh/devtoolset-7/enable
source /software/LNS/intel-parallel-studio-xe-2017/bin/compilervars.sh -arch intel64 -platform linux
source /software/LNS/intel-parallel-studio-xe-2017/impi/2017.4.239/bin64/mpivars.sh
source /software/LNS/intel-parallel-studio-xe-2017/mkl/bin/mklvars.sh intel64

mpicxx  main.cpp SimplexAbstract.cpp SimplexAlpha.cpp LocalSearch.cpp Matrix.cpp  Covering.hpp  RCC.hpp -DMPICH_IGNORE_CXX_SEEK  -std=c++0x -O3 -mavx  -pthread -o absMotion-intel-ignore-crecxxseek

date;

time mpirun -np 96 ./absMotion-intel-ignore-crecxxseek > absMotion-intel-ignore-crecxxseek.out

date;

exit 0;

