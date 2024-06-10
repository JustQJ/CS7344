#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>



// p2p communication to implement MPI_Reduce
//process i send data to process size/2+i, if the size is odd, the first process won't join the communication.
// Don't exclude the final process rather than the first process because we can't gurantee the the final process's rank in the tree reduce 
int main(int argc, char *argv[]){

    // Initialize MPI
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int local_sum;
    int gloabl_sum; //only for rank 0 for MPI_Reduce
    double elapsed_time;

    //warm up
    MPI_Barrier(MPI_COMM_WORLD);
    local_sum = rank + 1;
    MPI_Reduce(&local_sum, &gloabl_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //################# 1. MPI_Reduce #################
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    local_sum = rank + 1;
    MPI_Reduce(&local_sum, &gloabl_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(rank == 0){
        printf("Origin MPI_Reduce: The sum is %d, and cost time is %10.6f \n", gloabl_sum, elapsed_time);
    }   

    //################# 2. Implementation MPI_Reduce by using Simple sending all data to processes 0 #################
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    local_sum = rank + 1;
    if(rank == 0){
        gloabl_sum = local_sum;
        for(int i=1; i<size; i++){
            int temp;
            MPI_Recv(&temp, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            gloabl_sum += temp;
        }
    }else{
         MPI_Send(&local_sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    elapsed_time += MPI_Wtime();
    if(rank == 0){
        printf("Implementation MPI_Reduce by sending all data to processes 0: The sum is %d, and cost time is %10.6f \n", gloabl_sum, elapsed_time);
    }   

    //################# 3. Implementation MPI_Reduce by using P2P communication with Tree Reduce #################
    MPI_Barrier(MPI_COMM_WORLD);
    local_sum = rank + 1;
    elapsed_time = -MPI_Wtime();
    //use MPI_Send and MPI_Recv to implement MPI_Reduce
    //exlude the first process if the size is odd
    int current_size = size;
    while(current_size>1){
        int half_size = current_size>>1;
        int half_size1 = (current_size+1)>>1;
        if(rank<half_size1 && rank+half_size>=half_size1){
            int temp;
            MPI_Recv(&temp, 1, MPI_INT, half_size+rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_sum += temp;
        }else if(rank>=half_size1){ //for the odd number of process, final process will not send data
            MPI_Send(&local_sum, 1, MPI_INT, rank-half_size, 0, MPI_COMM_WORLD);
            break;
        }
        current_size = half_size1;
    }

    elapsed_time += MPI_Wtime();
    if(rank == 0){
        printf("Implementation MPI_Reduce by using P2P communication with Tree Reduce: The sum is %d and the cost time is %10.6f \n", local_sum, elapsed_time);
    }

    

    MPI_Finalize();
    return 0;
}


