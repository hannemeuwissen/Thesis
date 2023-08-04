mpirun -np 1 ./caa_gen -v b.txt -m 2500 -z 2000 -d 16 -s 2 -t
mpirun -np 2 ./caa_gen -v b.txt -m 5000 -z 2000 -d 16 -s 2 -t
mpirun -np 4 ./caa_gen -v b.txt -m 10000 -z 2000 -d 16 -s 2 -t
mpirun -np 8 ./caa_gen -v b.txt -m 20000 -z 2000 -d 16 -s 2 -t
echo "--------------------------------------------------------------"
mpirun -np 1 ./caa_gen -v b.txt -m 2500 -z 2000 -d 16 -s 4 -t
mpirun -np 2 ./caa_gen -v b.txt -m 5000 -z 2000 -d 16 -s 4 -t
mpirun -np 4 ./caa_gen -v b.txt -m 10000 -z 2000 -d 16 -s 4 -t
mpirun -np 8 ./caa_gen -v b.txt -m 20000 -z 2000 -d 16 -s 4 -t