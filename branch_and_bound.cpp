#include <iostream>
#include <vector>
#include <time.h>

using namespace std;

vector<vector<int> > dist;
vector<int> visited;
vector<int> path;
vector<int> best_path;
int min_cost = INT_MAX;

void wsp(int curr_pos, int n, int cost, int count) {
    if (count == n) {
        if (cost < min_cost) {
            min_cost = cost;
            for (int i = 0; i < n; i++) {
                best_path[i] = path[i];
            }
        }
        return;
    }
    for (int i = 0; i < n; i++) {
        if (!visited[i] && dist[curr_pos][i]) {
            path[count] = i;
            visited[i] = 1;
            wsp(i, n, cost + dist[curr_pos][i], count + 1);
            visited[i] = 0;
        }
    }
}

int main() {
    int N;
    scanf("%d", &N);
    dist.resize(N, vector<int>(N));
    visited.resize(N);
    path.resize(N);
    best_path.resize(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            scanf("%d", &dist[i][j]);
            dist[j][i] = dist[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        visited[i] = 0;
    }
    clock_t tStart = clock();

    for (int i = 0; i < N; i++) {
        path[0] = i;
        visited[i] = 1;
        wsp(i, N, 0, 1);
        visited[i] = 0;
    }

    cout << "Minimum cost : " << min_cost << endl;
    cout << "Path : ";
    for (int i = 0; i < N; i++) {
        printf("%d ", best_path[i]);
    }
    double time_taken = (double) (clock() - tStart) / CLOCKS_PER_SEC;
    printf("\nTime taken: %f seconds\n", time_taken);

    return 0;
}


//#include <iostream>
//#include <vector>
//using namespace std;
//
//const int N = 10;
//int dist[N][N];
//
//int visited[N];
//int path[N];
//int best_path[N];
//int min_cost = INT_MAX;
//
//void tsp(int curr_pos, int n, int cost, int count) {
//    if (count == n && dist[curr_pos][0]) {
//        min_cost = min(min_cost, cost + dist[curr_pos][0]);
//        for (int i = 0; i < n; i++) {
//            best_path[i] = path[i];
//        }
//        return;
//    }
//    for (int i = 0; i < n; i++) {
//        if (!visited[i] && dist[curr_pos][i]) {
//            path[count] = i;
//            visited[i] = 1;
//            tsp(i, n, cost + dist[curr_pos][i], count + 1);
//            visited[i] = 0;
//        }
//    }
//}
//
//int main() {
////    scanf("%d", &N);
//    for (int i = 0; i < N; i++) {
//        for(int j =0; j< i; j++) {
//            scanf("%d", &dist[i][j]);
//            dist[j][i] = dist[i][j];
//        }
//    }
//
//    //print matrix
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            printf("%d ", dist[i][j]);
//        }
//        printf("\n");
//    }
//
//    for (int i = 0; i < N; i++) {
//        visited[i] = 0;
//    }
//    path[0] = 0;
//    visited[0] = 1;
//    tsp(0, N, 0, 1);
//    cout << "Minimum cost : " << min_cost << endl;
//    cout << "Path Taken : ";
//    for (int i = 0; i < N; i++) {
//        cout << best_path[i] << " ";
//    }
//    return 0;
//}
