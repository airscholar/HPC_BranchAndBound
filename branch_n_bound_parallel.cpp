#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <time.h>
#include <vector>

using namespace std;

int myrank, num_procs;
clock_t tStart;
void wsp(int *dist, bool *visited, int *path, int *best_path, int &min_cost, int curr_pos, int n, int cost, int count) {
    if (count == n) {
        if (cost < min_cost) {
            min_cost = cost;
            for (int i = 0; i < n; i++) {
                best_path[i] = path[i];
            }
        }
        return;
    }
    for (int i = 0; i < n; i++) {
        if (!visited[i] && dist[curr_pos * n + i]) {
            //check if the cost is greater than the min_cost
            if (cost + dist[curr_pos * n + i] > min_cost) {
                continue;
            }
            path[count] = i;
            visited[i] = true;
            wsp(dist, visited, path, best_path, min_cost, i, n, cost + dist[curr_pos * n + i], count + 1);
            visited[i] = false;
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Request request;
    int N = 0;
    string filename;

    if(argv[2] == NULL) {
        cout << "Please enter file filename" << endl;
        return 1;
    } else
        filename = argv[2];

    //open filename
    FILE *fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        printf("Error opening file");
        return 1;
    }
    fscanf(fp, "%d", &N);

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int *dist = new int[N * N];
    bool *visited = new bool[N];
    int *path = new int[N];
    int *best_path = new int[N];
    int min_cost = numeric_limits<int>::max();
    if (myrank == 0) {
        //fill the matrix with 0
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                dist[i * N + j] = 0;
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                fscanf(fp, "%d", &dist[i * N + j]);
                dist[j * N + i] = dist[i * N + j];
            }
        }

        tStart = clock();
    }

    MPI_Bcast(dist, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    int chunk_size = N / num_procs;
    int starting_city = myrank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;
    //the last process will take care of the remaining cities
    if (myrank == num_procs - 1) {
        ending_city = N - 1;
    }
    for (int i = starting_city; i <= ending_city; i++) {
        path[0] = i;
        visited[i] = true;
        wsp(dist, visited, path, best_path, min_cost, i, N, 0, 1);
        visited[i] = false;
    }

    int local_min_cost = min_cost;
    int local_best_path[N];
    for (int i = 0; i < N; i++) {
        local_best_path[i] = best_path[i];
    }

    printf("Process %d: min_cost = %d, best_path = ", myrank, local_min_cost);
    for (int i = 0; i < N; i++) {
        printf("%d ", local_best_path[i]);
    }
    printf("\n");

    MPI_Allreduce(&local_min_cost, &min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (local_min_cost == min_cost) {
        //send to rank 0
        printf("Sending from %d to 0\n", myrank);
        MPI_Isend(local_best_path, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
    }

    if (myrank == 0) {
        MPI_Recv(best_path, N, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Best path: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\n");
        printf("Min cost: %d\n", min_cost);
        //time taken formatted in MM:SS format
        int minutes = (int) ((clock() - tStart) / CLOCKS_PER_SEC) / 60;
        int seconds = (int) ((clock() - tStart) / CLOCKS_PER_SEC) % 60;
        printf("Time taken: %d:%d\n", minutes, seconds);
    }

    delete[] dist;
    delete[] visited;
    delete[] path;
    delete[] best_path;
    MPI_Finalize();
    return 0;
}