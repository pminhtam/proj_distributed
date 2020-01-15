# proj_distributed
mpicc main.c -o graph
mpirun -np 2 graph

mpicc floyd_warshall.c -o floyd -lm
mpirun -np 1 floyd < input12x12.txt
mpirun --hostfile host -np 9 floyd < input12x12.txt