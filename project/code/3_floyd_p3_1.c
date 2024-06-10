#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

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




//using row and colomn split
int main()
{
    MPI_Init(NULL, NULL);
    int rank;
    int p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    
    //p should be a perfect square, and N should be divisible by sqrt(p)
    int sqrt_p = (int)sqrt((double)p);
    if(sqrt_p*sqrt_p!=p){
        printf("p is not a perfect square\n");
        MPI_Finalize();
        return -1;
    }
    int sub_size = N / sqrt_p;
    if(sub_size*sqrt_p!=N){
        printf("N is not divisible by sqrt(p)\n");
        MPI_Finalize();
        return -1;
    }
    int i, j, k;
    init_array();

   
   
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = -MPI_Wtime();

#pragma scop
    //using colomn and row split, each process handle the sub_size rows and sub_size columns (n/sqrt(p) * n/sqrt(p))
    int p_x = rank / sqrt_p;
    int p_y = rank % sqrt_p;
    int start_row = p_x*sub_size;
    int end_row = (p_x+1)*sub_size;
    int start_col = p_y*sub_size;
    int end_col = (p_y+1)*sub_size;
    //printf("rank: %d in (%d, %d), start_row: %d, end_row: %d, start_col: %d, end_col: %d\n", rank, p_x, p_y ,start_row, end_row, start_col, end_col);
    //split the commucation world, row_comm for row communication, col_comm for col communication
    MPI_Comm row_comm, col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, p_x, rank, &row_comm); //all processes in the same row p_x are in the same communicator
    MPI_Comm_split(MPI_COMM_WORLD, p_y, rank, &col_comm); //all processes in the same col p_y are in the same communicator

    double temp_row[sub_size]; //receive the row from other processes or send the row to other processes
    double temp_col[sub_size]; //receive the col from other processes or send the col to other processes
    for(k = 0; k < N; k++) {
        int root_row = k / sub_size; //the process that has the k-th row
        int root_col = k / sub_size; //the process that has the k-th col
        if(root_row==p_x){
            for(j = 0; j < sub_size; j++) { //copy the k-th row to temp_row
                temp_row[j] = A[k][j+start_col]; 
            }
        }
        if(root_col==p_y){
            for(i = 0; i < sub_size; i++) { //copy the k-th col to temp_col
                temp_col[i] = A[i+start_row][k]; 
            }
        }
        //broadcast the k-th row to all processes in the same colomn
        MPI_Bcast(temp_row, sub_size, MPI_DOUBLE, root_row, col_comm);
        //broadcast the k-th col to all processes in the same row
        MPI_Bcast(temp_col, sub_size, MPI_DOUBLE, root_col, row_comm);
        for(i = start_row; i < end_row; i++) {
            for(j = start_col; j < end_col; j++) {
                A[i][j] = min(A[i][j], temp_col[i-start_row] + temp_row[j-start_col]);
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
    //define new data type for (n/sqrt(p)) * (n/sqrt(p)) submatrix
    MPI_Datatype submatrix;
    MPI_Type_vector(sub_size, sub_size, N, MPI_DOUBLE, &submatrix);
    MPI_Type_commit(&submatrix);
    if(rank==0){
        for(i = 1; i < p; i++){
            int temp_p_x = i / sqrt_p;
            int temp_p_y = i % sqrt_p;
            int temp_start_row = temp_p_x*sub_size;
            int temp_start_col = temp_p_y*sub_size;
            
            MPI_Recv(&A[temp_start_row][temp_start_col], 1, submatrix, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
    }else{
        MPI_Send(&A[start_row][start_col], 1, submatrix, 0, 0, MPI_COMM_WORLD);
    }

    print_array();
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Finalize();

    return 0;
}
