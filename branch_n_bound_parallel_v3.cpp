#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <time.h>
#include <vector>
#include <cstring>

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

void parallel_tsp(int *dist, int n, int starting_city, int ending_city, int &min_cost, int *best_path) {
    bool *visited = new bool[n];
    int *path = new int[n];

    min_cost = numeric_limits<int>::max();
    for (int i = starting_city; i <= ending_city; i++) {
        path[0] = i;
        memset(visited, false, sizeof(bool) * n);
        visited[i] = true;

        if (min_cost + dist[starting_city * n + i] > min_cost) {
            continue;
        }

        wsp(dist, visited, path, best_path, min_cost, i, n, 0, 1);
    }
    // MPI all reduce on min_cost
    MPI_Allreduce(MPI_IN_PLACE, &min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

//    if (local_min_cost == min_cost) {
//        // MPI all reduce on best_path
//        MPI_Allreduce(MPI_IN_PLACE, best_path, n, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Request request;
    int N = 0;
    string filename;

    if (argv[2] == NULL) {
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
    int min_cost, best_path[N];

    if (myrank == 0) {
        //fill the matrix with 0
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                dist[i * N + j] = 0;
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                if (i == j) {
                    dist[i * N + j] = 0;
                } else {
                    fscanf(fp, "%d", &dist[i * N + j]);
                    dist[j * N + i] = dist[i * N + j];
                }
            }
        }

        fclose(fp);
    }

    MPI_Bcast(&dist[0], N * N, MPI_INT, 0, MPI_COMM_WORLD);

    int chunk_size = N / num_procs;
    int starting_city = myrank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;

    if (myrank == num_procs - 1) {
        ending_city = N - 1;
    }

    printf("Process: %d, Total size: %lu, Chunk size: %d, starting city: %d, ending city: %d\n",myrank, N, chunk_size, starting_city, ending_city);

    tStart = clock();
    parallel_tsp(dist, N, starting_city, ending_city, min_cost, best_path);
    double time_taken = (double) (clock() - tStart) / CLOCKS_PER_SEC;

    if (myrank == 0) {
        printf("Minimum cost : %d\n", min_cost);
        printf("Path : ");
        for (int i = 0; i < N; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\nTime taken: %.2fs\n", time_taken);
    }

    MPI_Finalize();

    return 0;
}