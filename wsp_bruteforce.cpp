#include <iostream>
#include <vector>
#include <algorithm>
#include <climits> // for INT_MAX

using namespace std;

int main(int argc, char *argv[]) {
    int N = 0; // number of cities
    // Accept the filename and starting city from the command line
    if (argc != 3 || string(argv[1]) != "-i") {
        cout << "Usage: " << argv[0] << " -i <filename>" << endl;
        exit(1);
    }
    const string filename = argv[2];

    // Open the file
    FILE *fp = fopen(filename.c_str(), "r");
    if (fp == nullptr) {
        cout << "Error opening file" << endl;
        exit(1);
    }

    // Read the number of cities
    fscanf(fp, "%d", &N);

    // Allocate memory for the distance matrix
    vector<vector<int> > dist(N, vector<int>(N));
    // Read the distance matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            fscanf(fp, "%d", &dist[i][j]);
            dist[j][i] = dist[i][j];
        }
    }
    fclose(fp);

    // Initialize path permutations
    vector<int> mat_size;
    for (int i = 0; i < N; i++) {
        mat_size.push_back(i);
    }

    int shortest_path_sum = INT_MAX;
    vector<int> shortest_path_vec(N - 1);

    // Perform permutations of the matrix to get all possible paths
    do {
        int running_sum = 0;
        // Calculate the sum of the path
        for (int i = 0; i < mat_size.size() - 1; i++) {
            running_sum += dist[mat_size[i]][mat_size[i + 1]];
        }

        // Update the shortest path and cost if necessary
        if (running_sum < shortest_path_sum) {
            shortest_path_sum = running_sum;
            shortest_path_vec = mat_size;
        }
    } while (next_permutation(mat_size.begin(), mat_size.end()));

    // Print the shortest path
    cout << "Shortest path: " << shortest_path_sum << endl;
    cout << "Path: ";
    for (int i = 0; i < shortest_path_vec.size(); i++) {
        cout << shortest_path_vec[i] << " ";
    }
    return 0;
}