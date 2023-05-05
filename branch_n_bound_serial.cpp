#include <iostream>
#include <cstdio>
#include <limits>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

double tStart; // start time
int N = 0; // number of cities
string filename; // file name
int START_PATH = 0; // start path

// function to find the best path using branch and bound
/// @param dist the distance matrix
/// @param visited the visited array
/// @param path the current path
/// @param best_path the best path
/// @param min_cost the min cost
/// @param n the number of cities
/// @param cost the current cost
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
/// @param starting_city the starting city
/// @return the list of starting cities
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
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " -i <filename> <starting city>" << endl;
        exit(1);
    } else {
        //accept the filename from the command line
        filename = argv[2];
        //accept the starting city from the command line
        START_PATH = atoi(argv[3]);
    }
    //open filename
    FILE *fp = fopen(filename.c_str(), "r");
    if (fp == nullptr) {
        printf("Error opening file");
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
    }

    // allocate memory for the distance matrix
    int *dist = new int[N * N];
    // allocate memory for the best path
    int *best_path = new int[N];
    // variable to store the minimum cost
    int min_cost = numeric_limits<int>::max();
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

    // list of starting cities
    vector<vector<int> > starting_cities;
    // generate the starting paths
    starting_cities = generate_starting_paths(START_PATH);

    // split the starting cities between the processes
    int chunk_size = starting_cities.size();
    // starting city for the current process
    int starting_city = 0;
    // ending city for the current process
    int ending_city = starting_city + chunk_size - 1;

    int *min_cost_array = new int[starting_cities.size()]; // array to store the min_cost of each process
    int *best_path_array = new int[starting_cities.size() * N];

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

        // store the min_cost and the best path for the current process
        min_cost_array[i] = min_cost;
        for (int j = 0; j < N; j++) {
            best_path_array[i * N + j] = best_path[j];
        }

        // Find the global minimum cost and best path among all processes
        min_cost_array[i] = min_cost;
        min_cost = *min_element(min_cost_array, min_cost_array + i);
        //get the index of the process with the minimum cost
        int idx = min_element(min_cost_array, min_cost_array + i) - min_cost_array;
        //copy the best path of the process with the minimum cost to the best_path vector
        for (int j = 0; j < N; j++) {
            best_path[j] = best_path_array[idx * N + j];
        }
    }

    //calculate the total time
    double total_time = (clock() - tStart) / CLOCKS_PER_SEC;

    //print the results
    printf("Num of processes: %d\n", 1);
    printf("Global Min cost: %d\n", min_cost);
    printf("Global Best path: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", best_path[i]);
    }
    printf("\n");
    //percentage efficiency
    printf("Efficiency: %.2f%%\n", (total_time / (total_time + 0)) * 100);

    //clean up
    free(dist);
    free(best_path);
    free(min_cost_array);
    free(best_path_array);

    return 0;
}