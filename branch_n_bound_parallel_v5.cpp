#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <time.h>
#include <vector>
#include <cstring>
#include <algorithm>

using namespace std;

int myrank, num_procs;
double tStart;
std::vector<int> best_path;


int estimateUpperBound(int* dist, std::vector<bool> &visited, int n, int currentCity) {
    int upperBound = 0;
    int nearestUnvisitedCity = -1;
    int nearestUnvisitedCityDistance = INT_MAX;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            upperBound += dist[currentCity * n + i];
            if (dist[currentCity * n + i] < nearestUnvisitedCityDistance) {
                nearestUnvisitedCity = i;
                nearestUnvisitedCityDistance = dist[currentCity * n + i];
            }
        }
    }
    upperBound += nearestUnvisitedCityDistance;
    return upperBound;
}


void wsp(int *dist, std::vector<bool> &visited, std::vector<int> &path, int &global_min_cost, int &min_cost, int n, int &cost) {
    if (path.size() == n) {
        if (cost < min_cost) {
            min_cost = cost;
            for (int i = 0; i < n; i++) {
                best_path[i] = path[i];
            }
        }
        return;
    }
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            //check if the cost is greater than the min_cost
            if (cost + dist[path.back() * n + i] > min_cost || cost + dist[path.back() * n + i] > global_min_cost) {
                // cut the branch
                break;
            }
            int upper_bound = estimateUpperBound(dist, visited, n, i);
            if (cost + upper_bound > global_min_cost) {
                break;
            }
            cost += dist[path.back() * n + i];
            path.push_back(i);
            visited[i] = true;
            wsp(dist, visited, path, global_min_cost, min_cost, n, cost);
            visited[i] = false;
            path.pop_back();
            cost -= dist[path.back() * n + i];
        }
    }
}

vector<int> generate_path_nearest_neighbor(int *dist, int N, int starting_city) {
    vector<int> path;
    path.push_back(starting_city);
    vector<bool> visited(N, false);
    visited[starting_city] = true;
    int currentCity = starting_city;
    for (int i = 0; i < N - 1; i++) {
        int nearestCity = -1;
        int nearestDistance = INT_MAX;
        for (int j = 0; j < N; j++) {
            if (!visited[j]) {
                int distance = dist[currentCity * N + j];
                if (distance < nearestDistance) {
                    nearestCity = j;
                    nearestDistance = distance;
                }
            }
        }
        path.push_back(nearestCity);
        visited[nearestCity] = true;
        currentCity = nearestCity;
    }
    return path;
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
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
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
    int min_cost = numeric_limits<int>::max();
    int local_min_cost = numeric_limits<int>::max();

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

    vector<vector<int> > starting_cities = generate_path(N, 0);
    tStart = MPI_Wtime();
    best_path.resize(N);

    int chunk_size = starting_cities.size() / num_procs;
    int starting_city = myrank * chunk_size;
    int ending_city = starting_city + chunk_size - 1;

    if (myrank == num_procs - 1) {
        ending_city = starting_cities.size() - 1;
    }
//    printf("Process %d: size %lu with starting city %d and ending city %d\n", myrank, starting_cities.size(), starting_city, ending_city);

    int check = 0;

    for (int i = starting_city; i <= ending_city; i++) {
        std::vector<int> path = starting_cities[i];

        //initialize visited array
        std::vector<bool> visited(N, false);
        //visited
        for (int j = 0; j < path.size(); j++) {
            visited[path[j]] = true;
        }
        //calculate the cost of the path
        int cost = 0;

        for (int j = 0; j < path.size() - 1; j++) {
            cost += dist[path[j] * N + path[j + 1]];
        }

        wsp(dist, visited, path, min_cost, local_min_cost, N, cost);

        int *min_cost_array = new int[num_procs];
        MPI_Allgather(&local_min_cost, 1, MPI_INT, min_cost_array, 1, MPI_INT, MPI_COMM_WORLD);
        int *best_path_array = new int[num_procs * N];
        MPI_Allgather(&best_path[0], N, MPI_INT, best_path_array, N, MPI_INT, MPI_COMM_WORLD);

        // Find the global minimum cost and best path among all processes
        min_cost = *min_element(min_cost_array, min_cost_array + num_procs);
        int idx = min_element(min_cost_array, min_cost_array + num_procs) - min_cost_array;
        for (int i = 0; i < N; i++) {
            best_path[i] = best_path_array[idx * N + i];
        }
        //broadcast the best path to all processes
        MPI_Bcast(&best_path[0], N, MPI_INT, idx, MPI_COMM_WORLD);

        local_min_cost = min_cost;

        //clean up memory
        delete[] min_cost_array;
        delete[] best_path_array;
    }

    if (myrank == 0) {
        std::cout << "Num of processes: " << num_procs << std::endl;
        std::cout << "Global Min cost: " << min_cost << std::endl;
        std::cout << "Global Best path: ";
        for (int i = 0; i < N; i++) {
            std::cout << best_path[i] << " ";
        }
        std::cout << std::endl;
        printf("\nTime taken: %f seconds\n", MPI_Wtime() - tStart);

    }

    //clean up
    delete[] dist;

    MPI_Finalize();
    return 0;
}