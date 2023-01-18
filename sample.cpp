
#include<stdio.h>
#include<mpi.h>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
    int myrank, npes;
//    double t1, t2;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    t1 = MPI_Wtime();
    int a;
    int b;

    scanf("%d", &a);
    //create matrix
    int matrix[a][a];

    for (int i = 0; i < a; i++) {
        for(int j =0; j< i; j++) {
            scanf("%d", &b);
            matrix[i][j] = b;
        }
    }
//
    //fill the matrix symmetrically
    for (int i = 0; i < a; i++) {
        for(int j =0; j< i; j++) {
            matrix[j][i] = matrix[i][j];
        }
    }

    vector<int> mat_size;
    for(int i=0; i<a ; i++){
        mat_size.push_back(i);
    }
//
    //perform permutations of the matrix to get all possible paths
    int shortest_path_sum = numeric_limits<int>::max();
    vector<int> shortest_path_vec;

    do {
        printf("Permutation: from process %d: ", myrank);
        for (int i = 0; i < mat_size.size(); i++) {
            printf("%d ", mat_size[i]);
        }
        printf("\n");
        int running_sum = 0;
        for(int i=0; i<a-1; i++){
            if(mat_size[i] < mat_size[i+1]){
                running_sum += matrix[mat_size[i]][mat_size[i+1]];
            } else {
                running_sum += matrix[mat_size[i+1]][mat_size[i]];
            }
        }

        if(shortest_path_sum > running_sum){
            shortest_path_sum = running_sum;
            shortest_path_vec = mat_size;
        }

    } while (next_permutation(mat_size.begin(), mat_size.end()));


    printf("Sum: %d\n", shortest_path_sum);
    printf("Path: ");
    for(int i=0; i<shortest_path_vec.size(); i++){
        printf("%d, ", shortest_path_vec[i]+1);
    }
////    return 0;
//    printf("\n\nHwllo from %d", myrank);
    MPI_Finalize();
}


//
// Created by Yusuf Ganiyu on 1/9/23.
//
