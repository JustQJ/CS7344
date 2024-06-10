#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define min(a, b) ((a) < (b) ? (a) : (b))

#define N 1024

double A[N][N];

void init_array()
{
    int i, j;

    srand(123);
    for(i= 0; i < N; i++) {
        for(j= 0; j < N; j++) {
            double ra = 1.0 + ((double)N * rand() / (RAND_MAX + 1.0));
            A[i][j] = ra;
        }
    }
    for(i= 0; i < N; i++) {
        A[i][i] = 0;
    }
}

void print_array()
{
    int i, j;

    for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
            fprintf(stderr, "%lf ", A[i][j]); 
        }
        fprintf(stderr, "\n");
    }
}




//using row split
int main()
{
    MPI_Init(NULL, NULL);
    int rank;
    int p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int sub_size = N / p;
    if(sub_size*p!=N){
        printf("N is not divisible by p\n");
        MPI_Finalize();
        return -1;
    }
    int i, j, k;
    init_array();

 
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = -MPI_Wtime();

#pragma scop
    //using row split
    int start_row = rank*sub_size;
    int end_row = (rank+1)*sub_size;
    // printf("rank: %d, start_row: %d, end_row: %d\n", rank, start_row, end_row);
    double temp_row[N]; //receive the row from other processes or send the row to other processes
    for(k = 0; k < N; k++) {
        int root = k / sub_size; //the process that has the k-th row
        if(rank==root){
            for(j = 0; j < N; j++) { //copy the k-th row to temp_row
                temp_row[j] = A[k][j]; 
            }
        }
        //broadcast the k-th row to all processes
        MPI_Bcast(temp_row, N, MPI_DOUBLE, root, MPI_COMM_WORLD);
        for(i = start_row; i < end_row; i++) {
            for(j = 0; j < N; j++) {
                A[i][j] = min(A[i][j], A[i][k] + temp_row[j]);
            }
        }

    }


#pragma endscop

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(rank==0){
        printf("Time: %lf\n", elapsed_time);
    }

    //merge all data to rank 0
    if(rank==0){
        for(i = 1; i < p; i++){
            MPI_Recv(A[i*sub_size], sub_size*N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }else{
        MPI_Send(A[start_row], sub_size*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    print_array();

    MPI_Finalize();

    return 0;
}
