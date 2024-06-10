#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW(id,p,n))
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))

int serial_consecutive_prime_number(int n){
    if(n<=3) return 0;
    char *mark =  (char*)malloc(sizeof(char)*n+1); 
    for(int i=0; i<=n; i++) mark[i] = 0;
    mark[0] = 1;
    mark[1] = 1;
    for(int i=2; i<=n; i++){
        if(mark[i]==1) continue;
        for(int j=i+i; j<=n; j+=i){
            mark[j] = 1;
        }
    }
    int count = 0;
    for(int i=3; i<n-1; i++){
        if(mark[i]==0 && mark[i+2]==0) count++;
    }
    free(mark);
    return count;
} // result is 8169


// int main(int argc, char *argv[]){
//     int n = 1000000;
//     // compute time
//     MPI_Init(&argc, &argv);
//     MPI_Barrier(MPI_COMM_WORLD);
//     double elapsed_time = -MPI_Wtime();
//     int real_count = serial_consecutive_prime_number(n);
//     elapsed_time += MPI_Wtime();
//     printf("The number of consecutive prime numbers is %d , computed by serial code with cost time: %10.6f. \n", real_count, elapsed_time);
//     MPI_Finalize();
// }

int main(int argc, char *argv[]){
    int n = 1000000;
    // int real_count = serial_consecutive_prime_number(n);
    // printf("The number of consecutive prime numbers is %d , computed by serial code.\n", real_count);

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    //measuring the time
    double elapsed_time = -MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    
    n = n/2; //only odd numbers
    int low_value, high_value, block_size;
    low_value = 3 + 2*BLOCK_LOW(rank, size, n-1); //start from 3, delete 1 and 2, n-1 because delete 1
    high_value = 3 + 2*BLOCK_HIGH(rank, size, n-1) + 2; //+1 for judege the consecutive number of the last number
    block_size = BLOCK_SIZE(rank, size, n-1)+1;
    if(3+(n-1)/size < (int) sqrt((double)n)){
        if(rank == 0) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }
    
    
    char *marked = (char *)malloc(block_size);
    if(marked == NULL){
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    int i;
    for(i=0; i<block_size; i++){
        marked[i] = 0;
    }
    int index;
    if(rank == 0){
        index = 0;
    }
    int prime = 3; 
    int first;
    do{
        if(prime*prime > low_value){
            first = prime*prime - low_value;
        }else{
            if(low_value%prime == 0){
                first = 0;
            }else{
                first = prime - (low_value%prime);
            }
        }
        if(first%2 == 1) first += prime; //first must be even, because we only consider odd numbers in the block

        for(i=first/2; i<block_size; i+=prime){
            marked[i] = 1;
        }
        if(rank == 0){
            while(marked[++index]);
            prime = 2*index + 3;
        }
        MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }while(prime*prime <= n*2);

    int count = 0;
    int current_number = low_value;
    //the actual start number of the block
    for(i=0; i<block_size-1; i++){ //最后一个数不需要判断，会在下一个block中判断
        if(!marked[i] && !marked[i+1]){ //连续两个odd都是prime
            count++;   
        }
    }
    if(rank==size-1){
        //最后一个block的最后一个数是多余的，所以如何其和前面的odd组成了连续的情况要减去
        if(marked[block_size-2] == 0 && marked[block_size-1] == 0){
            count--;
        }
    }
    int global_count;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(rank == 0){
        printf("The number of consecutive prime numbers is %d , computed by mpi code with cost time: %10.6f.\n", global_count, elapsed_time);
    }

    free(marked);
    MPI_Finalize();


}