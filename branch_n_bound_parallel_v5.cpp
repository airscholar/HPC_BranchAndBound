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
vector<int> best_path; // best path
int N = 0; // number of cities
string filename; // file name
double running_idle_time=0; // running idle time
double total_idle_time=0; // total idle time
const int START_PATH = 0; // start path

// function to find the best path using branch and bound
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
            if (cost + dist[path.back() * n + i] > min_cost) {
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

// function to generate the starting paths
vector<vector<int> > generate_starting_paths(int N, int starting_city) {
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

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

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
    // read the number of cities
    fscanf(fp, "%d", &N);
    // allocate memory for the distance matrix
    best_path.resize(N);
    // start time
    tStart = MPI_Wtime();
    // broadcast the number of cities to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // allocate memory for the distance matrix
    int *dist = new int[N * N];
    // variable to store the minimum cost
    int min_cost = numeric_limits<int>::max();

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

    double start_block_time = MPI_Wtime(); // start timer for the block
    MPI_Bcast(&dist[0], N * N, MPI_INT, 0, MPI_COMM_WORLD);
    double end_block_time = MPI_Wtime(); // end timer for the block
    running_idle_time += (end_block_time - start_block_time); // calculate the running idle time

    // list of starting cities
    vector<vector<int> > starting_cities;
    // generate the starting paths
    starting_cities = generate_starting_paths(N, START_PATH);

    // split the starting cities between the processes
    int chunk_size = starting_cities.size() / num_process;
    // starting city for the current process
    int starting_city = my_rank * chunk_size;
    // ending city for the current process
    int ending_city = starting_city + chunk_size - 1;
    // if the number of starting cities is not divisible by the number of processes
    // then the last process will take the remaining starting cities
    if (my_rank == num_process - 1) {
        //min of the last process
        ending_city = min(ending_city, (int) starting_cities.size() - 1);
    }
    // start timer for the block
    start_block_time = MPI_Wtime();
    // iterate over the cities for the current process
    for (int i = starting_city; i <= ending_city; i++) {
        // variable to store the current stating path
        vector<int> path = starting_cities[i];
        // variable to store the visited cities
        bool visited[N];
        memset(visited, 0, sizeof(visited));
        //visited
        for (int j = 0; j < path.size(); j++) {
            visited[path[j]] = true;
        }
        //calculate the cost of the path
        int cost = 0;
        for (int j = 0; j < path.size()-1; j++) {
            cost += dist[path[j] * N + path[j + 1]];
        }
        // call the function to find the best path recursively
        wsp(dist, visited, path, min_cost, N, cost);

        int *min_cost_array = new int[num_process]; // array to store the min_cost of each process
        int *best_path_array = new int[num_process * N]; // array to store the best_path of each process
        //gather the min_cost of all the processes
        MPI_Allgather(&min_cost, 1, MPI_INT, min_cost_array, 1, MPI_INT, MPI_COMM_WORLD);
        //gather the best_path of all the processes
        MPI_Allgather(&best_path[0], N, MPI_INT, best_path_array, N, MPI_INT, MPI_COMM_WORLD);

        // Find the global minimum cost and best path among all processes
        min_cost = *min_element(min_cost_array, min_cost_array + num_process);
        //get the index of the process with the minimum cost
        int idx = min_element(min_cost_array, min_cost_array + num_process) - min_cost_array;
        //copy the best path of the process with the minimum cost to the best_path vector
        for (int j = 0; j < N; j++) {
            best_path[j] = best_path_array[idx * N + j];
        }
        //broadcast the best path to all processes
        MPI_Bcast(&best_path[0], N, MPI_INT, idx, MPI_COMM_WORLD);

        //clean up memory
        delete[] min_cost_array;
        delete[] best_path_array;
    }
    end_block_time = MPI_Wtime(); //end of the idle time calculation

    //calculate the total idle time
    running_idle_time += (end_block_time - start_block_time);
    //print the process block time
    running_idle_time = MPI_Wtime() - tStart - running_idle_time;

//    printf("Process %d: running idle time = %f\n", my_rank, running_idle_time);

    // sum up the total idle time of all the processes into the final_idle_time variable
    MPI_Reduce(&running_idle_time, &total_idle_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


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
        printf("Idle time: %f seconds\n", total_idle_time);
        printf("Total time: %f seconds\n", total_time + total_idle_time);
        //percentage efficiency
        printf("Efficiency: %.2f%%\n", (total_time / (total_time + total_idle_time)) * 100);
    }

    //clean up
    delete[] dist;

    MPI_Finalize();
    return 0;
}