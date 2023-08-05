#!/bin/bash

#SBATCH -J thesis-small-test
#SBATCH -n 8
#SBATCH -t 00:02:00
#SBATCH -p compute
#SBATCH -A MSCHPC

#SBATCH --output=/home/users/mschpc/2022/meuwissh/Thesis/ResultsLonsdale/%x.out
#SBATCH --error=/home/users/mschpc/2022/meuwissh/Thesis/ResultsLonsdale/%x.error

# go to working directory
cd /home/users/mschpc/2022/meuwissh/Thesis/Parallel

# load up necessary modules
# module load cports gcc/12.1.0-gnu intel openmpi
# module load apps intel-oneapi/2022.1.0 
# module load apps intel-oneapi/2022.1.0 
# module load gcc/12.2.0 
module load apps intel-oneapi/2022.1.0 openmpi 

# launch code
make caa_gen
mpicc -o caa_gen caa_gen.c parse.o graph.o sparse.o tsqr_mpi.o bgs.o matrix.o hess.o -Wextra -Wall -march=opteron \
        -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -liomp5 -lpthread -lm
mpirun -np 1 ./caa_gen -v b.txt -m 400 -z 200 -d 8 -s 2 -t -a
# mpirun -np 2 ./caa_gen -v b.txt -m 400 -z 200 -d 8 -s 2 -t
# mpirun -np 4 ./caa_gen -v b.txt -m 400 -z 200 -d 8 -s 2 -t
# mpirun -np 8 ./caa_gen -v b.txt -m 400 -z 200 -d 8 -s 2 -t
make clean

