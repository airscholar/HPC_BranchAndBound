#include <iostream>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

int my_rank, num_process; // rank of process and number of process
double tStart; // start time
int N = 0; // number of cities
string filename; // file name
double running_comm_time=0; // running comm time
double total_comm_time=0; // total comm time
int START_PATH = 0; // start path
int upper_bound = numeric_limits<int>::max();
int lower_bound = 0;

// function to find the best path using branch and bound
inline void wsp(int* &dist, bool *visited, vector<int> &path, int* &best_path, int &min_cost, int n, int &cost) {
    if (path.size() == n) {
        if (cost < min_cost) {
            min_cost = cost;
            std::move(path.begin(), path.end(), best_path);
            // update the lower bound to the min_cost
            lower_bound = min_cost;
        }
        return;
    }

    //check unvisited nodes only
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            //check if the cost is greater than the upper bound
            int current_cost = cost + dist[path.back() * n + i];
            if (current_cost >= upper_bound){
                // cut the branch
                break;
            }
            cost += current_cost;
            path.push_back(i);
            visited[i] = true;
            wsp(dist, visited, path, best_path, min_cost, n, cost);
            visited[i] = false;
            path.pop_back();
            cost -= current_cost;
        }
    }
}

// function to generate the starting paths
inline vector<vector<int> > generate_starting_paths(int N, int starting_city) {
    vector<vector<int> > starting_cities;

    // generate all possible starting cities
    for (int j = 0; j < N; j++) {
        // skip the starting city
        if (j == starting_city) continue;
        for (int k = 0; k < N; k++) {
            // skip the starting city
            if (k == starting_city || k == j) continue;
            for (int l = 0; l < N; l++) {
                // skip the starting city
                if (l == starting_city || l == j || l == k) continue;
                vector<int> current_path;
                current_path.push_back(starting_city);
                current_path.push_back(j);
                current_path.push_back(k);
                current_path.push_back(l);
                // add the current path to the list of starting cities
                starting_cities.push_back(current_path);
            }
        }
    }
    return starting_cities;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv); // initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // get the rank of process
    MPI_Comm_size(MPI_COMM_WORLD, &num_process); // get number of processes
    tStart = MPI_Wtime(); // start the timer
    // read the data from file
    if (my_rank == 0) {
        // code to read data from file
        //accept the starting city from the command line
        //open filename
        FILE *fp = fopen(filename.c_str(), "r");
        if (fp == nullptr) {
            printf("Error opening file");
            return 1;
        }
        // read the number of cities
        fscanf(fp, "%d", &N);
        // allocate memory for the distance matrix
        int *dist = new int[N * N];
        // allocate memory for the best path
        int *best_path = new int[N];
        // variable to store the minimum cost
        int min_cost = numeric_limits<int>::max();
        int running_idle_time = 0; // running idle time
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
                    dist[j * N + i] = dist[i * N + j]; // fill the matrix symmetrically
                }
            }
        }

        fclose(fp);
    }

// broadcast the data to all the processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&filename[0], filename.length(), MPI_CHAR, 0, MPI_COMM_WORLD);

// distribute the work among the processes
    int work_per_process = N / num_process;
    int start = my_rank * work_per_process;
    int end = (my_rank + 1) * work_per_process - 1;

// loop through each starting city
    for (int starting_city = start; starting_city <= end; starting_city++) {
        // generate starting paths for the starting city
        vector<vector<int> > starting_cities = generate_starting_paths(N, starting_city);

        // loop through each starting path
        for (vector<int> starting_path : starting_cities) {
            // allocate memory for the best path
            int* best_path = new int[N];
            memset(best_path, 0, N * sizeof(int));
            int min_cost = numeric_limits<int>::max();

            // allocate memory for visited nodes
            bool *visited = new bool[N];
            memset(visited, false, N * sizeof(bool));

            // allocate memory for distances
            int* dist = new int[N * N];
            memset(dist, 0, N * N * sizeof(int));

            // get the distances


            // code to get the distances

            // call the branch and bound algorithm
            vector<int> path;
            int cost = 0;
            path.push_back(starting_path[0]);
            visited[starting_path[0]] = true;
            path.push_back(starting_path[1]);
            visited[starting_path[1]] = true;
            path.push_back(starting_path[2]);
            visited[starting_path[2]] = true;
            path.push_back(starting_path[3]);
            visited[starting_path[3]] = true;
            cost = dist[starting_path[0] * N + starting_path[1]] + dist[starting_path[1] * N + starting_path[2]] + dist[starting_path[2] * N + starting_path[3]];
            wsp(dist, visited, path, best_path, min_cost, N, cost);

            // free the memory
            delete[] visited;
            delete[] dist;
            delete[] best_path;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // wait for all processes to finish
    if(my_rank == 0) {
        double tEnd = MPI_Wtime();
        cout << "Time taken to find the best path: " << tEnd - tStart << " seconds." << endl;
        cout << "Total communication time: " << total_comm_time << " seconds." << endl;
    }

    MPI_Finalize();
    return 0;
}
