#include <iostream>
#include <cstdio>
#include <limits>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

double tStart; // start time
vector<int> best_path; // best path
int N = 0; // number of cities
string filename; // file name
double running_idle_time=0; // running idle time
double total_idle_time=0; // total idle time
const int START_PATH = 2; // start path

// function to calculate the upper bound of the path
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
    tStart = clock();
    // allocate memory for the distance matrix
    int *dist = new int[N * N];
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
    starting_cities = generate_starting_paths(N, START_PATH);

    // iterate over the cities for the current process
    for (int i = 0; i < starting_cities.size(); i++) {
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
        wsp(dist, visited, path, min_cost, N, cost);
    }

    printf("Num of processes: %d\n", 1);
    printf("Global Min cost: %d\n", min_cost);
    printf("Global Best path: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", best_path[i]);
    }
    printf("\n");

    double total_time = (clock() - tStart)/CLOCKS_PER_SEC;
    printf("\nTime taken: %f seconds\n", total_time);
    printf("Idle time: %f seconds\n", total_idle_time);
    printf("Total time: %f seconds\n", total_time + total_idle_time);
    //percentage efficiency
    printf("Efficiency: %.2f%%\n", (total_time / (total_time + total_idle_time)) * 100);

    // free the memory
    delete[] dist;

    return 0;
}