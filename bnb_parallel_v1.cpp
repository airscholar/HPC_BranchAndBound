#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>

using namespace std;

const int MAX_N = 100;
int dist[MAX_N][MAX_N];
bool visited[MAX_N];
int path[MAX_N];
int best_path[MAX_N];
int min_cost = numeric_limits<int>::max();
int myrank, num_process;

void wsp(int curr_pos, int n, int cost, int count) {
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
        if (!visited[i] && dist[curr_pos][i]) {
            path[count] = i;
            visited[i] = true;
            wsp(i, n, cost + dist[curr_pos][i], count + 1);
            visited[i] = false;
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    int N = 0;
    if (myrank == 0) {
        scanf("%d", &N);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                scanf("%d", &dist[i][j]);
                dist[j][i] = dist[i][j];
            }
        }
        printf("Matrix:\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%d ", dist[i][j]);
            }
            printf("\n");
        }
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dist[0][0], N * N, MPI_INT, 0, MPI_COMM_WORLD);

    int chunk_size = N / num_process;
    int starting_city = myrank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;
    if (myrank == num_process - 1) {
        ending_city = N - 1;
    }
    for (int i = starting_city; i <= ending_city; i++) {
        path[0] = i;
        visited[i] = true;
        wsp(i, N, 0, 1);
        visited[i] = false;
    }

    int local_min_cost = min_cost;
    int local_best_path[N];
    for (int i = 0; i < N; i++) {
        local_best_path[i] = best_path[i];
    }

    MPI_Allreduce(&local_min_cost, &min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (myrank == 0) {
        printf("Minimum cost : %d ", min_cost);
        printf("Path : ");
        for (int i = 0; i < N; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\n");
    }
    MPI_Finalize();
}
