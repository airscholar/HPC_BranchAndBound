#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <time.h>

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
            path[count] = i;
            visited[i] = true;
            //check the cost against minimum cost
            if (cost + dist[curr_pos * n + i] > min_cost) {
                visited[i] = false;
                continue;
            }
            wsp(dist, visited, path, best_path, min_cost, i, n, cost + dist[curr_pos * n + i], count + 1);
//            MPI_Allreduce(&min_cost, &min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            visited[i] = false;
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    int N = 0;
    if (myrank == 0) {
        scanf("%d", &N);
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int *dist = new int[N*N];
    bool *visited = new bool[N];
    int *path = new int[N];
    int *best_path = new int[N];
    int min_cost = numeric_limits<int>::max();

    if (myrank == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                scanf("%d", &dist[i*N + j]);
                dist[j*N + i] = dist[i*N + j];
            }
        }
        tStart = clock();

//        printf("Matrix:\n");
//        for (int i = 0; i < N; i++) {
//            for (int j = 0; j < N; j++) {
//                printf("%d ", dist[i*N + j]);
//            }
//            printf("\n");
//        }
    }
    MPI_Bcast(dist, N*N, MPI_INT, 0, MPI_COMM_WORLD);

    int chunk_size = N / num_procs;
    int starting_city = myrank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;
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

    MPI_Allreduce(&local_min_cost, &min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (myrank == 0) {
        double time_taken = (double) (clock() - tStart) / CLOCKS_PER_SEC;
        printf("Minimum cost : %d ", min_cost);
        printf("Path : ");
        for (int i = 0; i < N; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\n");
        printf("\nTime taken: %f seconds\n", time_taken);

    }
    delete[] dist;
    delete[] visited;
    delete[] path;
    delete[] best_path;
    MPI_Finalize();
    return 0;
}
