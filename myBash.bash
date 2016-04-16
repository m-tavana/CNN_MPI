#!/bin/sh
#BSUB -J tavana
#BSUB -o output_file
#BSUB -e error_file
#BSUB -n 1
#BSUB -q ht-10g
#BSUB cwd /home/khavaritavana.m/cnn_MPI/

work=/home/khavaritavana.m/cnn_MPI/

cd $work
mpirun -np 2 -prot -TCP -lsf cnn
