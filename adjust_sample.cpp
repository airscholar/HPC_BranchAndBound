#include <stdio.h>
#include <mpi.h>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {
    int myrank, npes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int a;

    scanf("%d", &a);
    //create matrix
    int matrix[a][a];

    //input matrix
    if (myrank == 0) {
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < a; j++) {
                scanf("%d", &matrix[i][j]);
            }
        }
    }

    //Broadcast matrix to all processes
    MPI_Bcast(matrix, a*a, MPI_INT, 0, MPI_COMM_WORLD);

    // Initialize path permutations
    vector<int> mat_size;
    for (int i = 1; i < a; i++) {
        mat_size.push_back(i);
    }

    // Initialize global best path
    int shortest_path_sum = numeric_limits<int>::max();
    vector<int> shortest_path_vec;
    int local_shortest_path_sum = shortest_path_sum;
    vector<int> local_shortest_path_vec;

    // Number of permutations per process
    int num_perms = 1;
    for (int i = 1; i <= a - myrank - 1; i++) num_perms *= i;

    // Perform permutations
    for (int i = 0; i < num_perms; i++) {
        int running_sum = matrix[0][mat_size[0]];
        for (int j = 1; j < a - 1; j++) {
            running_sum += matrix[mat_size[j - 1]][mat_size[j]];
        }
        running_sum += matrix[mat_size[a - 2]][0];
        if (local_shortest_path_sum > running_sum) {
            local_shortest_path_sum = running_sum;
            local_shortest_path_vec = mat_size;
        }
        next_permutation(mat_size.begin(), mat_size.end());
    }

    // Gather local best path from all processes
    MPI_Allreduce(&local_shortest_path_sum, &shortest_path_sum, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Check for global best path
    if (shortest_path_sum == local_shortest_path_sum) {
        shortest_path_vec = local_shortest_path_vec;
    }

    if (myrank == 0) {
        printf("Sum: %d\n", shortest_path_sum);
        printf("Path: ");
        printf("1, ");
        for (int i = 0; i < shortest_path_vec.size(); i++) {
            printf("%d, ", shortest_path_vec[i] + 1);
        }
        printf("1\n");
    }

    MPI_Finalize();
    return 0;
}
