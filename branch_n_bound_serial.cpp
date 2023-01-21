#include <iostream>
#include <cstdio>
#include <limits>
#include <vector>
#include <algorithm>

using namespace std;

double tStart;
std::vector<int> best_path;

int estimateUpperBound(int* dist, int *visited, int n, int currentCity) {
    int nearestUnvisitedCity = -1;
    int nearestUnvisitedCityDistance = numeric_limits<int>::max();
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            if (dist[currentCity * n + i] < nearestUnvisitedCityDistance) {
                nearestUnvisitedCity = i;
                nearestUnvisitedCityDistance = dist[currentCity * n + i];
            }
        }
    }
    return nearestUnvisitedCityDistance;
}

void wsp(int *dist, int *visited, std::vector<int> &path, int &global_min_cost, int &min_cost, int n, int &cost) {
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
            if (cost + dist[path.back() * n + i] > min_cost || cost + dist[path.back() * n + i] > global_min_cost || min_cost > global_min_cost) {
                // cut the branch
                break;
            }

            if (cost + estimateUpperBound(dist, visited, n, i) >  global_min_cost) {
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

    int *dist = new int[N * N];
    int min_cost = numeric_limits<int>::max();
    int local_min_cost = numeric_limits<int>::max();

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

    vector<vector<int> > starting_cities = generate_path(N, 0);
    tStart = clock();
    best_path.resize(N);


    int check = 0;

    for (int i = 0; i <= starting_cities.size(); i++) {
        std::vector<int> path = starting_cities[i];

        //initialize visited array
//        std::vector<bool> visited(N, false);
        int *visited = new int[N];
        for (int i = 0; i < N; i++) {
            visited[i] = false;
        }

        //visited
        for (int j = 0; j < path.size(); j++) {
            visited[path[j]] = true;
        }
        //calculate the cost of the path
        int cost = 0;

        for (int j = 0; j < path.size() - 1; j++) {
            cost += dist[path[j] * N + path[j + 1]];
        }

        if (min_cost > local_min_cost) break; // cut the branch
        wsp(dist, visited, path, min_cost, local_min_cost, N, cost);

        for (int i = 0; i < N; i++) {
            best_path[i] = path[i];
        }

        local_min_cost = min_cost;

    }
    std::cout << "Global Min cost: " << min_cost << std::endl;
    std::cout << "Global Best path: ";
    for (int i = 0; i < N; i++) {
        std::cout << best_path[i] << " ";
    }
    std::cout << std::endl;
    printf("\nTime taken: %f seconds\n", clock() - tStart);


    //clean up
    delete[] dist;

    return 0;
}