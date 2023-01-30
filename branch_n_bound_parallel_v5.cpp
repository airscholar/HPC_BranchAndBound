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

// function to find the best path using branch and bound
inline void wsp(int* &dist, bool *visited, vector<int> &path, int* &best_path, int &min_cost, int n, int &cost) {
    if (path.size() == n) {
        if (cost < min_cost) {
            min_cost = cost;
            std::move(path.begin(), path.end(), best_path);
        }
        return;
    }

    //check unvisited nodes only
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            //check if the cost is greater than the min_cost
            if (cost + dist[path.back() * n + i] > min_cost){
                // cut the branch
                break;
            }
            cost += dist[path.back() * n + i];
            path.push_back(i);
            visited[i] = true;
            wsp(dist, visited, path, best_path, min_cost, n, cost);
            visited[i] = false;
            path.pop_back();
            cost -= dist[path.back() * n + i];
        }
    }
}

// function to generate the starting paths
inline vector<vector<int> > generate_starting_paths(int starting_city) {
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

    if (argc != 4) {
        cout << "Usage: " << argv[0] << " -i <filename> <starting city>" << endl;
        exit(1);
    } else {
        //accept the filename from the command line
        filename = argv[2];
        //accept the starting city from the command line
        START_PATH = atoi(argv[3]);
    }

    //accept the starting city from the command line
    //open filename
    FILE *fp = fopen(filename.c_str(), "r");
    if (fp == nullptr) {
        printf("Error opening file\n");
        exit(1);
    }
    // read the number of cities
    fscanf(fp, "%d", &N);

    if(N < 4){
        printf("Number of cities should be greater than 3\n");
        exit(1);
    }
    else if (START_PATH > N-1){
        printf("Starting city should be less than the number of cities\n");
        exit(1);
    }else if(num_process > N){
        printf("Number of processes should be less than the number of cities\n");
        exit(1);
    }

    // allocate memory for the distance matrix
    int *dist = new int[N * N];
    // allocate memory for the best path
    int *best_path = new int[N];
    // variable to store the minimum cost
    int min_cost = numeric_limits<int>::max();
    int running_idle_time = 0; // running idle time
    tStart = MPI_Wtime();
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

    // list of starting cities
    vector<vector<int> > starting_cities;
    // generate the starting paths
    starting_cities = generate_starting_paths(START_PATH);

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

    int iteration_size = 0;
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
        for (int j = 0; j < path.size() - 1; j++) {
            cost += dist[path[j] * N + path[j + 1]];
        }
        // call the function to find the best path recursively
        wsp(dist, visited, path, best_path, min_cost, N, cost);

        int *min_cost_array = new int[num_process]; // array to store the min_cost of each process
        int *best_path_array = new int[num_process * N]; // array to store the best_path of each process
        double start_block_time = MPI_Wtime();
        //gather the min_cost of all the processes
        MPI_Allgather(&min_cost, 1, MPI_INT, min_cost_array, 1, MPI_INT, MPI_COMM_WORLD);
        running_comm_time += MPI_Wtime() - start_block_time; // calculate the communication time

        start_block_time = MPI_Wtime();
        //gather the best_path of all the processes
        MPI_Allgather(&best_path[0], N, MPI_INT, best_path_array, N, MPI_INT, MPI_COMM_WORLD);
        running_comm_time += MPI_Wtime() - start_block_time; // calculate the communication time
        // Find the global minimum cost and best path among all processes
        min_cost = *min_element(min_cost_array, min_cost_array + num_process);
        //get the index of the process with the minimum cost
        int idx = min_element(min_cost_array, min_cost_array + num_process) - min_cost_array;
        //copy the best path of the process with the minimum cost to the best_path vector
        for (int j = 0; j < N; j++) {
            best_path[j] = best_path_array[idx * N + j];
        }
        iteration_size++;
        //clean up memory
        free(min_cost_array);
        free(best_path_array);
    }

    running_comm_time = running_comm_time / iteration_size;
//    printf("Process %d Running comm time: %f size: %d start size: %d\n", my_rank, running_comm_time, iteration_size, starting_cities.size());

    // sum up the total comm time of all the processes into the final_comm_time variable
    MPI_Reduce(&running_comm_time, &total_comm_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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
        printf("Communication time: %f seconds\n", total_comm_time);
//        printf("Idle time: %f seconds\n", total_time - total_comm_time);
        printf("Total time: %f seconds\n", total_time + total_comm_time);
        //percentage efficiency
        printf("Efficiency: %.2f%%\n", (total_time / (total_time + total_comm_time)) * 100);
    }

    //clean up
    free(dist);

    MPI_Finalize();
    return 0;
}