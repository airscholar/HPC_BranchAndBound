#scp branch_n_bound_parallel_v4.cpp s395207@crescent.central.cranfield.ac.uk:HPC/branch_n_bound_parallel_v4.cpp #&& scp -r ./queue s395207@crescent.central.cranfield.ac.uk:HPC/queue
scp branch_n_bound_serial.cpp s395207@crescent.central.cranfield.ac.uk:HPC/queue/branch_n_bound_serial.cpp #&& scp -r ./queue s395207@crescent.central.cranfield.ac.uk:HPC/queue

#ssh s395207@crescent.central.cranfield.ac.uk 'cd HPC && ./script.sh'
#scp s395207@crescent.central.cranfield.ac.uk:HPC/mpi.sub queue/mpi.sub