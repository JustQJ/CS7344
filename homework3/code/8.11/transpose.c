#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>




//mpirun -np 4 ./transpose 3

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int rank, work_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &work_size);

    int pow = 13;
    if (argc==2){
        pow = atoi(argv[1]);
    }
    int n = 1<<pow;
    int block_size = n/work_size; // assume n is divisible by work_size
    int* local_rows = (int*)malloc(sizeof(int)*block_size*n); //each process has block_size rows, range from rank*block_size to (rank+1)*block_size-1
    if (rank==0){
        printf("n is %d, block_size is %d\n", n, block_size);
    }
    
    //initlize local_rows
    for(int i=0; i<block_size; i++){
        for(int j=0; j<n; j++){
            local_rows[i*n+j] = (rank*block_size+i)*n+j;
        }
    }
    

    //print local_rows
    // printf("rank %d local rows\n ", rank);
    // for(int i=0; i<block_size; i++){
    //     for(int j=0; j<n; j++){
    //         printf("%d ", local_rows[i*n+j]);
    //     }
    //     printf("\n");
    // }

    int* local_columns = (int*)malloc(sizeof(int)*block_size*n); //store the transpose of local colomns
     // redefine the column type such that receive data canbe stored in local_columns
    MPI_Datatype columnType, columnResized;
    MPI_Type_vector(n, 1, n, MPI_INT, &columnType);
    MPI_Type_create_resized(columnType, 0, sizeof(int), &columnResized);
    MPI_Type_commit(&columnResized);

    MPI_Barrier(MPI_COMM_WORLD);
    //record time
    double elapsed_time = -MPI_Wtime();

    //transpose, using sendrecv
    for(int i=0; i<work_size; i++){
        for(int j=0; j<block_size; j++){
            MPI_Sendrecv(&local_rows[j*n+i*block_size], block_size, MPI_INT, i, 0, 
                        &local_columns[i*block_size+j], 1, columnResized, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(rank==0){
        printf("cost time is %10.6f s\n", elapsed_time);
    }

    //print the result
    // printf("rank %d results\n", rank);
    // for(int i=0; i<block_size; i++){
    //     for(int j=0; j<n; j++){
    //         printf("%d ", local_columns[i*n+j]);
    //     }
    //     printf("\n");
    // }

    free(local_rows);
    free(local_columns);
    MPI_Finalize();
    return 0;      

}