#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <vector>
#include <algorithm>

using namespace std;

int my_rank, num_process;
double tStart;
vector<int> best_path;
vector<int> idle_times;

int estimateUpperBound(const int* dist, bool* &visited, int n, int currentCity) {
    int nearestUnvisitedCityDistance = numeric_limits<int>::max();
    for (int i = 0; i < n; i++) {
        if (!(visited[i])) {
            if (dist[currentCity * n + i] < nearestUnvisitedCityDistance) {
                nearestUnvisitedCityDistance = dist[currentCity * n + i];
            }
        }
    }
    return nearestUnvisitedCityDistance;
}

void wsp(int* &dist, bool *visited, vector<int> &path, int &min_cost, int n, int &cost) {
    if (path.size() == n) {
        if (cost < min_cost) {
            min_cost = cost;
            for (int i = 0; i < n; i++) {
                best_path[i] = path[i];
            }
        }
        return;
    }

    //check unvisited nodes only
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            //check if the cost is greater than the min_cost
            if (cost + dist[path.back() * n + i] > min_cost ||
                    cost + estimateUpperBound(dist, visited, n, i) >  min_cost) {
                // cut the branch
                break;
            }
            cost += dist[path.back() * n + i];
            path.push_back(i);
            visited[i] = true;
            wsp(dist, visited, path, min_cost, n, cost);
            visited[i] = false;
            path.pop_back();
            cost -= dist[path.back() * n + i];
        }
    }
}

vector<vector<int> > generate_path(int N, int starting_city) {
    vector<vector<int> > starting_cities;

    for (int j = 0; j < N; j++) {
        if (j == starting_city) continue;
        for (int k = 0; k < N; k++) {
            if (k == starting_city || k == j) continue;
            for (int l = 0; l < N; l++) {
                if (l == starting_city || l == k || l == j) continue;
                vector<int> cities;
                cities.push_back(starting_city);
                cities.push_back(j);
                cities.push_back(k);
                cities.push_back(l);
                starting_cities.push_back(cities);
            }
        }
    }
    return starting_cities;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    int N = 0;
    string filename;

    if (argv[2] == nullptr) {
        cout << "Please enter file filename" << endl;
        return 1;
    } else
        filename = argv[2];

//open filename
    FILE *fp = fopen(filename.c_str(), "r");
    if (fp == nullptr) {
        printf("Error opening file");
        return 1;
    }

    fscanf(fp, "%d", &N);

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int *dist = new int[N * N];
    int min_cost = numeric_limits<int>::max();
    double total_blocked_time;

    if (my_rank == 0) {
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

    vector<vector<int> > starting_cities = generate_path(N, 0);
    tStart = MPI_Wtime();
    best_path.resize(N);

    int chunk_size = starting_cities.size() / num_process;
    int starting_city = my_rank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;

    if (my_rank == num_process - 1) {
        ending_city = starting_cities.size() - 1;
    }

    for (int i = starting_city; i <= ending_city; i++) {
        vector<int> path = starting_cities[i];

        bool visited[N];
        memset(visited, 0, sizeof(visited));

        //visited
        for (int j : path) {
            visited[j] = true;
        }
        //calculate the cost of the path
        int cost = 0;

        for (int j = 0; j < path.size() - 1; j++) {
            cost += dist[path[j] * N + path[j + 1]];
        }

        wsp(dist, visited, path, min_cost, N, cost);

        int *min_cost_array = new int[num_process];
        double start = MPI_Wtime();
        MPI_Allgather(&min_cost, 1, MPI_INT, min_cost_array, 1, MPI_INT, MPI_COMM_WORLD);
        int *best_path_array = new int[num_process * N];
        MPI_Allgather(&best_path[0], N, MPI_INT, best_path_array, N, MPI_INT, MPI_COMM_WORLD);
        double blocked_time = start - MPI_Wtime();

        MPI_Reduce(&blocked_time, &total_blocked_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        // Find the global minimum cost and best path among all processes
        min_cost = *min_element(min_cost_array, min_cost_array + num_process);
        int idx = min_element(min_cost_array, min_cost_array + num_process) - min_cost_array;

        for (int j = 0; j < N; j++) {
            best_path[j] = best_path_array[idx * N + j];
        }

        //broadcast the best path to all processes
        MPI_Bcast(&best_path[0], N, MPI_INT, idx, MPI_COMM_WORLD);

        //clean up memory
        delete[] min_cost_array;
        delete[] best_path_array;

    }

    if (my_rank == 0) {
        printf("Num of processes: %d\n", num_process);
        printf("Global Min cost: %d\n", min_cost);
        printf("Global Best path: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\n");

        double total_time = MPI_Wtime() - tStart;
        printf("\nTime taken: %f seconds\n", total_time);
        printf("Idle time: %f seconds\n", total_blocked_time);
    }

    //clean up
    delete[] dist;

    MPI_Finalize();
    return 0;
}