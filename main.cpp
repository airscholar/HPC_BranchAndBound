#include <iostream>
#include <algorithm>
#include <vector>
#include <mpi.h>

using namespace std;

const int INF = 1e9;

int n; // Number of cities
int best_cost = INF; // The best solution cost found so far
vector<int> best_path; // The best solution path found so far
vector<vector<int> > dist; // The distance matrix

// Function to calculate the cost of a given path
int calc_cost(const vector<int>& path) {
    int cost = 0;
    for (int i = 0; i < n; i++) {
        cost += dist[path[i]][path[i + 1]];
    }
    cost += dist[path[n]][path[0]];
    return cost;
}

void tsp(vector<int> path, int current_cost) {
    if (path.size() == n) {
        if (current_cost < best_cost) {
            best_cost = current_cost;
            best_path = path;
        }
        return;
    }

    for (int i = 0; i < n; i++) {
        if (find(path.begin(), path.end(), i) != path.end()) continue;
        path.push_back(i);
        current_cost += dist[path[path.size() - 2]][i];
        tsp(path, current_cost);
        current_cost -= dist[path[path.size() - 2]][i];
        path.pop_back();
    }
}

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Input the number of cities and the distance matrix
    if (world_rank == 0) {
        scanf("%d", &n);
        dist.resize(n, vector<int>(n));
        for (int i = 0; i < n; i++) {
            for(int j =0; j< i; j++) {
                scanf("%d", &dist[i][j]);
            }
        }
//
        //fill the matrix symmetrically
        for (int i = 0; i < n; i++) {
            for(int j =0; j< i; j++) {
                dist[j][i] = dist[i][j];
            }
        }

        //print matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%d ", dist[i][j]);
            }
            printf("\n");
        }
    }

    // Broadcast the number of cities to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory for the distance matrix on all processes
    dist.resize(n, vector<int>(n));

    // Broadcast the distance matrix to all processes
    MPI_Bcast(&dist[0][0], n * n, MPI_INT, 0, MPI_COMM_WORLD);

    // Start the tsp algorithm from each city in parallel
    for (int i = world_rank; i < n; i++)    {
        vector<int> path;
        path.push_back(i);
        tsp(path, 0);
    }

//    // Gather all the best solutions found by each process
    MPI_Allreduce(MPI_IN_PLACE, &best_cost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    int best_path_size = best_path.size();
    MPI_Allreduce(MPI_IN_PLACE, &best_path_size, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (best_path.size() > best_path_size) {
        best_path.resize(best_path_size);
    }

    int *temp = &best_path[0];
    MPI_Allreduce(MPI_IN_PLACE, &temp, best_path.size(), MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Print the best solution
    if (world_rank == 0) {
        cout << "Shortest route: ";
        for (int i = 0; i < best_path.size(); i++) {
            cout << best_path[i] << " ";
        }
        cout << endl;
        cout << "Cost: " << best_cost << endl;
    }

    MPI_Finalize();
}

