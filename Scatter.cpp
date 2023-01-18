#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {
    int myrank, npes;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int a[100], b[25];

    for (int i = 0; i < 100; i++) a[i] = i;

    MPI_Scatter(a, 25, MPI_INT, b, 25, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Gather(a, 10, MPI_INT, b, 10, MPI_INT, i, MPI_COMM_WORLD);

    for (int i = 0; i < 25; i++) printf("data scattered %d from process %d\n", b[i], myrank);

    MPI_Finalize();
}
