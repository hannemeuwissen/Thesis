#!/bin/bash

#SBATCH -J thesis-small-test
#SBATCH -n 64
#SBATCH -t 02:00:00
#SBATCH -p compute
#SBATCH -A MSCHPC

#SBATCH --output=/home/users/mschpc/2022/meuwissh/Thesis/ResultsLonsdale/%x.out
#SBATCH --error=/home/users/mschpc/2022/meuwissh/Thesis/ResultsLonsdale/%x.error

# go to working directory
cd /home/users/mschpc/2022/meuwissh/Thesis/Parallel

# load up necessary modules
module load apps intel-oneapi/2022.1.0
module load gcc/12.2.0 intel/18.0.5 openmpi

# launch code
make -f MakefileLonsdale caa_gen

echo "----------------------------------------------------------"
echo "Blocksize 2"
echo "----------------------------------------------------------"
mpirun -np 1 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
mpirun -np 2 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
mpirun -np 4 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
mpirun -np 8 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
mpirun -np 16 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
mpirun -np 32 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
mpirun -np 64 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 2 -t
echo "----------------------------------------------------------"
echo "----------------------------------------------------------"
echo "Blocksize 4"
echo "----------------------------------------------------------"
mpirun -np 1 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"
mpirun -np 2 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"
mpirun -np 4 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"
mpirun -np 8 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"
mpirun -np 16 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"
mpirun -np 32 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"
mpirun -np 64 ./caa_gen -v b.txt -m 40000 -z 4000 -d 16 -s 4 -t
echo "----------------------------------------------------------"

make clean

