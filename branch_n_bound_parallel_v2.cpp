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
//        wsp(dist, visited, path, best_path, min_cost, i, n, dist[starting_city * n + i], 1);
        //MPI all reduce on min_cost
        MPI_Allreduce(&min_cost, &min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    }
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
        tStart = clock();
    }

    MPI_Bcast(dist, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    vector<vector<int> > starting_cities;
    //generate path
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            for (int k = 0; k < N; k++) {
                if (k == i || k == j) continue;
                for (int l = 0; l < N; l++) {
                    if (l == i || l == j || l == k) continue;
                    vector<int> temp;
                    temp.push_back(i);
                    temp.push_back(j);
                    temp.push_back(k);
                    temp.push_back(l);
                    starting_cities.push_back(temp);
                }
            }
        }
    }

    // divide the problem into chunks
    int chunk_size = starting_cities.size() / num_procs;
    int starting_city = myrank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;
    if (myrank == num_procs - 1) {
        ending_city = starting_cities.size() - 1;
    }
    printf("Process: %d, Total size: %lu, Chunk size: %d, starting city: %d, ending city: %d\n",myrank, starting_cities.size(), chunk_size, starting_city, ending_city);

//    int chunk_size = starting_cities.size() / num_procs;
//    int starting_city = myrank * chunk_size;
//    int ending_city = starting_city + chunk_size - 1;
//    if (myrank == num_procs - 1) {
//        ending_city = starting_cities.size() - 1;
//    }
//
//    //divide the problem into chunks
//    int chunk_size = N / num_procs;
//    int starting_city = myrank * chunk_size;
//    int ending_city = starting_city + chunk_size - 1;
//    if (myrank == num_procs - 1) {
//        ending_city = N - 1;
//    }
//    each process will handle different starting city
    parallel_tsp(dist, N, starting_city, ending_city, min_cost, best_path);
    //MPI all reduce on min_cost

    printf("Process: %d, min_cost: %d, best path: ", myrank, min_cost);
    for (int i = 0; i < N; i++) {
        printf("%d ", best_path[i]);
    }
    printf("\n");
    //reduce the results from all the processes
    int global_min_cost = numeric_limits<int>::max();
    int global_best_path[N];
    MPI_Allreduce(&min_cost, &global_min_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);


    printf("process %d best path: ", myrank);
    for (int i = 0; i < N; i++) {
        printf("%d ", best_path[i]);
    }
    printf(" cost: %d", min_cost);
    printf(" global min cost: %d", global_min_cost);
    printf("\n");

    //calculate the cost from the path
    int cost = 0;
    for (int i = 0; i < N - 1; i++) {
        cost += dist[best_path[i] * N + best_path[i + 1]];
    }

    if (global_min_cost == min_cost) {
        for (int i = 0; i < N; i++) {
            global_best_path[i] = best_path[i];
        }
        //send to rank 0
        printf("Sending from %d to 0\n", myrank);
        MPI_Isend(global_best_path, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
    }

    if (myrank == 0) {
        MPI_Recv(best_path, N, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        double time_taken = (double) (clock() - tStart) / CLOCKS_PER_SEC;
        printf("Num of procs: %d\n", num_procs);
        printf("Minimum cost: %d\n", global_min_cost);
//        printf("Size: %lu\n", starting_cities.size());
        printf("Best path: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\nTime taken: %f seconds\n", time_taken);
    }

    MPI_Finalize();
    return 0;
}
