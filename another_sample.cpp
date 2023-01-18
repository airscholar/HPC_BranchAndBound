#include <iostream>
#include <vector>
#include <algorithm>
//#include <mpi.h>
#include <time.h>

using namespace std;

const int MAX_CITIES = 17;
int n;
int dist[MAX_CITIES][MAX_CITIES];
vector<int> path;
int best_cost = 1e9;
vector<int> best_path;

void tsp(int current_city, int visited_cities, int distance) {
    // base case: all cities have been visited
    if (visited_cities == (1 << n) - 1) {
        // add distance from last city to starting city
        distance += dist[current_city][0];

        if (distance <= best_cost) {
            best_path = path;
//            best_cost = distance;
        }
        return;
    }
    for (int next_city = 1; next_city < n; next_city++) {
        // check if city has not been visited
        if (!(visited_cities & (1 << next_city))) {
            path.push_back(next_city);
            tsp(next_city, visited_cities | (1 << next_city), distance + dist[current_city][next_city]);
            path.pop_back();
        }
    }
}

int main() {
//    int n = 0;
    vector<vector<int> > d;

    scanf("%d", &n);
    d.resize(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            scanf("%d", &dist[i][j]);
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
    printf("\n");

    clock_t tStart = clock();

    for(int i = 0; i < n; i++) {
        path.push_back(i);
        tsp(i, 1, 0);
        path.pop_back();
    }

    printf("\nPath: ");
    for (int i = 0; i < best_path.size(); i++) {
        printf("%d ", best_path[i]);
    }
//    calculate the cost using the path
    int sum = 0;
    for (int i = 0; i < best_path.size() - 1; i++) {
        int curr = best_path[i];
        int next = best_path[(i + 1)];
//        printf("\n%d %d = %d", curr, next, dist[curr][next]);
        sum += dist[curr][next];
    }

    printf("\nCost = %d \n", sum);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

//    printf("\nBest cost: %d\n", best_cost);
    return 0;
}