#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#define N 1024

double A[N][N][N];
// double B[N][N][N];

void init_array()
{
    int i, j, k;
    for (k=0; k<N; k++) {
        for (j=0; j<N; j++) {
            for (i=0; i<N; i++) {
                A[k][j][i] = (1+(i*j+k)%1024)/3.0;
            }
        }
    }
}



void print_array()
{
    int i, j, k;

    for (k=0; k<N; k++) {
        for (j=0; j<N; j++) {
            for (i=0; i<N; i++) {
                fprintf(stderr, "%lf ", A[k][j][i]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    }
}





int main()
{
    
    MPI_Init(NULL, NULL);
    int rank;
    int p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int sub_size = N / p;

    if (sub_size*p != N) {
        printf("N must be divisible by the number of processes\n");
        MPI_Finalize();
        return 1;
        
    }

    //using serial code to assure the correctness of the parallel code
    // if(rank==0){    
    //     serial_code();
    // }


    int i, j, k;

    init_array();
    
    // start time
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = -MPI_Wtime();

   
    //split in x and y dim, because the data is stored in row major order, so the z dim is the most contiguous   
    //there are p loop in k dim, each loop, a process with send and receive data from it neighbors
    //send and receive data to neighbors
    //send to the left ((rank-1+p)%p), receive from the right ((rank+1)%p)
    //x:i, y:j, z:k, for z, no data exchange needed
    int left_p = (rank-1+p)%p; // 0<-3, 1<-0, 2<-1, 3<-2
    int right_p = (rank+1)%p; // 0->1, 1->2, 2->3, 3->0

    //first part
    int start_i = rank * sub_size; 
    int end_i = start_i + sub_size;
    int start_j = 0;
    int end_j = start_j + sub_size;
    for(int t = 0; t < p; t++){
        for(i=start_i; i<end_i; i++){
            for(j=start_j; j<end_j; j++){
                for(k=1; k<N; k++){
                    A[i][j][k] = A[i][j][k] * 0.4 - A[i][j][k-1] * 0.6;
                }
            }
        }
        start_i = (end_i)%N;
        end_i = start_i + sub_size;
        start_j = end_j;
        end_j = start_j + sub_size;
  
    }

    


    //second part iterate over the y dim
    MPI_Request requests[2];
    MPI_Status statuses[2];

    MPI_Datatype subarray_block1;
    MPI_Type_vector(sub_size, N, N*N, MPI_DOUBLE, &subarray_block1);
    MPI_Type_commit(&subarray_block1);

    start_i = rank * sub_size; 
    end_i = start_i + sub_size;
    start_j = 1; //start from 1,0 has no previous element
    end_j = sub_size;
    for(int t = 0; t < p; t++){
        for(i=start_i; i<end_i; i++){
            for(j=start_j; j<end_j; j++){
                for(k=0; k<N; k++){
                    A[i][j][k] = A[i][j][k] * 0.5 - A[i][j-1][k] * 0.5;
                
                }
            }
        }
        if(t==p-1){
            break;
        }

        //send to the left ((rank-1+p)%p), receive from the right ((rank+1)%p)
        //send the A[start_j:end_j][end_k-1][0:N], each time send a array[N] to the left
        
        // for(int ii=0; ii<sub_size; ii++){
        //     MPI_Isend(A+start_location+ii*N*N, N, MPI_DOUBLE, left_p, 0, MPI_COMM_WORLD, &requests[ii]);
        // }
        //printf("start_j %d, start_location:  %d for process %d \n ",start_j, start_location, rank);
        MPI_Isend(A[start_i][end_j-1], 1, subarray_block1, left_p, 0, MPI_COMM_WORLD, &requests[0]);
        //receive from the right
        start_i = (end_i)%N;
        end_i = start_i + sub_size;
        //store data at A[start_j:end_j][end_k-1][0:N]
        
       
        
        MPI_Irecv(A[start_i][end_j-1], 1, subarray_block1, right_p, 0, MPI_COMM_WORLD, &requests[1]);
        
        MPI_Waitall(2, requests, statuses);
        start_j = end_j;
        end_j = start_j + sub_size;
  
    }


    
    //third part, iterate along x
    start_j = ((p-rank)%p) * sub_size; // loop on x dim, the order in y dim should be 0,p-1, p-2, ..., 1, so process i has (p-i)%p -th block, %p is mainly for process 0
    end_j = start_j + sub_size;
    start_i = 1;
    end_i = sub_size;
    for(int t = 0; t < p; t++){
        for(i=start_i; i<end_i; i++){
            for(j=start_j; j<end_j; j++){
                for(k=0; k<N; k++)
                {
                    A[i][j][k] = A[i][j][k] * 0.6 - A[i-1][j][k] * 0.4;
                
                }
            }
        }
        if(t==p-1){
            break;
        }

        //send to the right , receive from the left 
        //send the A[end_i-1][start_j:end_j][0:N], can directly send all the data in one go, sub_size*N
        
        MPI_Isend(A[end_i-1][start_j], sub_size*N, MPI_DOUBLE, right_p, 0, MPI_COMM_WORLD, &requests[0]);
        
        //receive from the right
        start_j = (end_j)%N;
        end_j = start_j + sub_size;
        //store data at A[end_i-1][start_j:end_j][0:N]    
        MPI_Irecv(A[end_i-1][start_j], sub_size*N, MPI_DOUBLE, left_p, 0, MPI_COMM_WORLD, &requests[1]);
        
        MPI_Waitall(2, requests, statuses);
        start_i = end_i;
        end_i = start_i + sub_size;
  
    }



    // end time
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(rank==0){
        printf("Time: %lf s\n", elapsed_time);
    }
    
    //merge all data to the process 0
    start_i = rank * sub_size; //this is a loop in y dim, each process move to next subarray at each iteration t 0->1->2->3
    end_i = start_i + sub_size;
    start_j = 0;
    end_j = start_j + sub_size;
    MPI_Datatype subarray_block;
    MPI_Type_vector(sub_size, N*sub_size, N*N, MPI_DOUBLE, &subarray_block);
    MPI_Type_commit(&subarray_block);
    for(int t=0; t<p; t++){
        if(rank==0){
            for(int i=1; i<p; i++){
                start_i = ((i+t)%p) * sub_size;
                end_i = start_i + sub_size;
                //receive A[start_j:end_j][start_k:end_k][0:N] from process i
                MPI_Recv(A[start_i][start_j], 1, subarray_block, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }else{
            MPI_Send(A[start_i][start_j], 1, subarray_block, 0, 0, MPI_COMM_WORLD);
            start_i = (end_i)%N;
            end_i = start_i + sub_size;
        }
        start_j = end_j;
        end_j = start_j + sub_size;

    }

    print_array();

    

    MPI_Finalize();

    return 0;
}
/*
No protocol specified
Time: 6.952225 s
*/
