#include <iostream>
#include <algorithm>
#include <vector>
#include <mpi.h>

using namespace std;

int n; // Number of cities
int best_cost = 1e9; // The best solution cost found so far
vector<int> best_path; // The best solution path found so far
vector<vector<int> > dist; // The distance matrix

void tree_algorithm(vector<int> path, int current_cost) {
    if (path.size() == n) {
        if (current_cost < best_cost) {
            best_cost = current_cost;
            best_path = path;
        }
        return;
    }

    for (int i = 0; i < n; i++) {
        //
        if (find(path.begin(), path.end(), i) != path.end()) continue;

        path.push_back(i);
        current_cost += dist[path[path.size() - 2]][i];
        tree_algorithm(path, current_cost);
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

    // Start the tree_algorithm algorithm from each city in parallel
    for (int i = world_rank; i < n; i++)    {
        vector<int> path;
        path.push_back(i);
        tree_algorithm(path, 0);
    }

    int best_path_size = best_path.size();
    MPI_Allreduce(MPI_IN_PLACE, &best_path_size, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Print the best solution
    if (world_rank == 0) {
        printf("Shortest route: ");
        for (int i = 0; i < best_path.size(); i++) {
            printf("%d ", best_path[i]);
        }

        //calculate the cost using the path
        int sum = 0;
        for (int i = 0; i < best_path.size()-1; i++) {
            int curr = best_path[i];
            int next = best_path[(i+1)];
            sum += dist[curr][next];
        }

        printf("\nCost = %d \n",sum);
    }

    MPI_Finalize();
    return 0;
}

