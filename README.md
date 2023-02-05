## Wandering Salesman Problem

This repository provides three solutions to the Wandering Salesman Problem using branch and bound algorithm.

The repository contains the following:
1. The brute force algorithm
2. The serial branch and bound algorithm 
3. The parallel branch and bound solution algorithm using MPI

## How to run
1. Install the OpenMPI library on your system here => [MPI](https://www.open-mpi.org/)

### Brute Force
```
mpicxx wsp_bruteforce.cpp -o wsp &&  ./wsp -i input/dist5 
```

### Branch and Bound Serial 
```
mpicxx branch_n_bound_serial.cpp -o main && mpirun -n {num_of_processes} ./main -i input/dist10 {start_path}
``` 

### Branch and Bound Parallel
```
mpicxx branch_n_bound_parallel.cpp -o main && mpirun -n {num_of_processes} ./main -i input/dist10 {start_path}
```

## Submitting jobs to the clustered queue
1. Copy the codes to the cluster
```
scp -r ./* {username}@{cluster_ip}:/home/{username}/
```

2. Compile the code on the cluster using Step 1 above
3. Navigate to either queue17, queue18 or queue19
```
cd queue17
```
4. Submit the job to the queue
```
sh ./submit.sh
```
5. Check the status of the job
```
qstat -u {username}
```
6. Check the output of the job
```
cat {job_name}.o{job_id}
```

## References
1. [MPI](https://www.open-mpi.org/)
2. [Branch and Bound](https://en.wikipedia.org/wiki/Branch_and_bound)
3. [Wandering Salesman Problem](https://en.wikipedia.org/wiki/Wandering_salesman_problem)
4. [Parallel Branch and Bound](https://www.researchgate.net/publication/2202000_Parallel_branch_and_bound_algorithms_for_the_travelling_salesman_problem)