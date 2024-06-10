#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// #define VERSION 3

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id) + 1, p, n) - BLOCK_LOW(id, p, n))
//basic version
int main_1(int argc, char **argv) {
    int count;        // Local prime count
    double elapsed_time; // Parallel execution time
    int first;        // Index of the first filter
    int global_count; // Global prime count
    int high_value;   // Highest value on this proc
    int id;           // Process ID number
    int index;        // Index of current prime
    int low_value;    // Lowest value on this proc
    int n;            // Sieving from 2 to 'n'
    int p;            // Number of processes
    int proc0_size;   // Size of proc 0's subarray
    int prime;        // Current prime
    int size;         // Elements in 'marked'

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime(); //record the start time, - for adding the end time
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoi(argv[1]);
    low_value = 2 + BLOCK_LOW(id, p, n-1);
    high_value = 2 + BLOCK_HIGH(id, p, n-1);
    size = BLOCK_SIZE(id, p, n-1);
    proc0_size = (n - 1) / p;

    if ((2 + proc0_size)< (int)sqrt((double)n)) { // let all sieve prime in the first process
        if (!id) printf("Too many processes: %d processes and %d number \n",p, n);
        MPI_Finalize();
        exit(1);
    }

    char *marked = (char *)malloc(size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    if(!id) 
        printf("Basic Version: \n");

    for (int i = 0; i < size; i++) marked[i] = 0;
    if (!id) index = 0;
    prime = 2;
    do {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        for (int i = first; i < size; i += prime) marked[i] = 1;
        if (!id) {
            while (marked[++index]);
            prime = index + 2;
        }
        MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    count = 0;
    for (int i = 0; i < size; i++)
        if (!marked[i]) count++;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %d primes less than or equal to %d\n", global_count-1, n);
        printf("Sieve of Eratosthenes took %10.6f seconds with %d processes. \n", elapsed_time, p);
    }
    MPI_Finalize();
    free(marked);
    return 0;
}

//improvement1: delete even integer
int main_2(int argc, char **argv) {
    int count;        // Local prime count
    double elapsed_time; // Parallel execution time
    int first;        // Index of the first filter
    int global_count; // Global prime count
    int high_value;   // Highest value on this proc
    int id;           // Process ID number
    int index;        // Index of current prime
    int low_value;    // Lowest value on this proc
    int n;            // Sieving from 2 to 'n'
    int p;            // Number of processes
    int proc0_size;   // Size of proc 0's subarray
    int prime;        // Current prime
    int size;         // Elements in 'marked'

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime(); //record the start time, - for adding the end time
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    //only odd number from 3 to n
    n = atoi(argv[1]);
    low_value = 3 + 2*BLOCK_LOW(id, p, n/2-1);
    high_value = 3 + 2*BLOCK_HIGH(id, p, n/2-1);
    size = BLOCK_SIZE(id, p, n/2-1); //size change half
    proc0_size = (n - 1) / p;

    if ((2 + proc0_size)< (int)sqrt((double)n)) { // let all sieve prime in the first process
        if (!id) printf("Too many processes: %d processes and %d number \n",p, n);
        MPI_Finalize();
        exit(1);
    }

    char *marked = (char *)malloc(size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    if(!id) 
        printf("Improvement1: delete even integer based on Basic Version: \n");

    for (int i = 0; i < size; i++) marked[i] = 0;
    if (!id) index = 0;
    prime = 3;
    do {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        if (first % 2 == 1) first += prime; //skip even number, start(low_value) from odd number, first should be even so that start+first is odd
        for (int i = first/2; i < size; i += prime) marked[i] = 1; //add prime actually is add 2*prime
        if (!id) {
            while (marked[++index]);
            prime = 2*index + 3;
        }
        MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    count = 0;
    for (int i = 0; i < size; i++)
        if (!marked[i]) count++;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %d primes between 3 and %d\n", global_count, n);
        printf("Sieve of Eratosthenes took %10.6f seconds with %d processes. \n", elapsed_time, p);
    }
    MPI_Finalize();
    free(marked);
    return 0;
}

//improvement2: eliminate the mpi boadcast for sieve prime
//each process use sequential method to find the prime between 3 and sqrt(n)
int main_3(int argc, char **argv) {
    int count;        // Local prime count
    double elapsed_time; // Parallel execution time
    int first;        // Index of the first filter
    int global_count; // Global prime count
    int high_value;   // Highest value on this proc
    int id;           // Process ID number
    int index;        // Index of current prime
    int low_value;    // Lowest value on this proc
    int n;            // Sieving from 2 to 'n'
    int p;            // Number of processes
    int proc0_size;   // Size of proc 0's subarray
    int prime;        // Current prime
    int size;         // Elements in 'marked'

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime(); //record the start time, - for adding the end time
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    //only odd number from 3 to n
    n = atoi(argv[1]);
    low_value = 3 + 2*BLOCK_LOW(id, p, n/2-1);
    high_value = 3 + 2*BLOCK_HIGH(id, p, n/2-1);
    size = BLOCK_SIZE(id, p, n/2-1); //size change half
    proc0_size = (n - 1) / p;

    if ((2 + proc0_size)< (int)sqrt((double)n)) { // let all sieve prime in the first process
        if (!id) printf("Too many processes: %d processes and %d number \n",p, n);
        MPI_Finalize();
        exit(1);
    }

    char *marked = (char *)malloc(size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    if(!id) 
        printf("Improvement2: eliminate the mpi boadcast for sieve prime based on delete even integer: \n");

    //find the prime between 3 and sqrt(n) , only odd number
    int temp_size = (int)sqrt(n)/2;
    char *marked_sieve_prime = (char *)malloc(temp_size);
    memset(marked_sieve_prime, 0, temp_size);
    for(int i=0; i<temp_size; i++){
        if(marked_sieve_prime[i]) continue;
        int start = 2*i+3;
        for(int j=start*start; j<=temp_size*2; j += 2*start){
            marked_sieve_prime[(j-3)/2] = 1;
        }
    }


    for (int i = 0; i < size; i++) marked[i] = 0;
    
    // prime = 3;
    // index = 0;
    // index = 3;
    for(int i=0; i<temp_size; i++){
        if(marked_sieve_prime[i]) continue;
        prime = 2*i + 3;
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        if (first % 2 == 1) first += prime; //skip even number, start from odd number, first should be even
        for (int j = first/2; j < size; j += prime) marked[j] = 1; //add prime actually is add 2*prime
    }       
        
    
    count = 0;
    for (int i = 0; i < size; i++)
        if (!marked[i]) count++;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %d primes between 3 and %d\n", global_count, n);
        printf("Sieve of Eratosthenes took %10.6f seconds with %d processes. \n", elapsed_time, p);
    }
    MPI_Finalize();
    free(marked);
    return 0;
}

//improvement3: reorganize the loop for a better cache hit ratio
int main_4(int argc, char **argv) {
    int count;        // Local prime count
    double elapsed_time; // Parallel execution time
    int first;        // Index of the first filter
    int global_count; // Global prime count
    int high_value;   // Highest value on this proc
    int id;           // Process ID number
    int index;        // Index of current prime
    int low_value;    // Lowest value on this proc
    int n;            // Sieving from 2 to 'n'
    int p;            // Number of processes
    int proc0_size;   // Size of proc 0's subarray
    int prime;        // Current prime
    int size;         // Elements in 'marked'

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime(); //record the start time, - for adding the end time
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    //only odd number from 3 to n
    n = atoi(argv[1]);
    low_value = 3 + 2*BLOCK_LOW(id, p, n/2-1);
    high_value = 3 + 2*BLOCK_HIGH(id, p, n/2-1);
    size = BLOCK_SIZE(id, p, n/2-1); //size change half
    proc0_size = (n - 1) / p;

    if ((2 + proc0_size)< (int)sqrt((double)n)) { // let all sieve prime in the first process
        if (!id) printf("Too many processes: %d processes and %d number \n",p, n);
        MPI_Finalize();
        exit(1);
    }

    char *marked = (char *)malloc(size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    if(!id) 
        printf("Improvement3: reorganize the loop for a better cache hit ratio based on eliminate the mpi boadcast and delete even integer: \n");

    //find the prime between 3 and sqrt(n) , only odd number
    int temp_size = (int)sqrt(n)/2;
    char *marked_sieve_prime = (char *)malloc(temp_size);
    memset(marked_sieve_prime, 0, temp_size);
    for(int i=0; i<temp_size; i++){
        if(marked_sieve_prime[i]) continue;
        int start = 2*i+3;
        for(int j=start*start; j<=temp_size*2; j += 2*start){
            marked_sieve_prime[(j-3)/2] = 1;
        }
    }


    for (int i = 0; i < size; i++) marked[i] = 0;
    
    //cache block size, 32KB
    int cache_size = 32*1024;
    int cache_block_low_value, cache_block_high_value;

    //firstly loop the range from low to high, each time loop the range of cache size
    for(cache_block_low_value=low_value; cache_block_low_value<=high_value; cache_block_low_value+=cache_size){
        cache_block_high_value = cache_block_low_value + cache_size-1;
        if(cache_block_high_value > high_value) cache_block_high_value = high_value;
        //compute the range of marked array
        int cache_block_low_index = (cache_block_low_value - low_value)/2;
        int cache_block_high_index = (cache_block_high_value - low_value)/2;
        //loop current block, set the sieve prime's multiple to 1
        for(int i=0; i<temp_size; i++){
            if(marked_sieve_prime[i]) continue;
            prime = 2*i + 3;
            if (prime * prime > cache_block_low_value)
                first = prime * prime - cache_block_low_value;
            else {
                if (!(cache_block_low_value % prime)) first = 0;
                else first = prime - (cache_block_low_value % prime);
            }
            if (first % 2 == 1) first += prime; //skip even number, start from odd number, first should be even
            for(int j=first/2; j<=cache_block_high_index-cache_block_low_index; j+=prime) marked[j+cache_block_low_index] = 1;
           
        }       
    }

    count = 0;
    for (int i = 0; i < size; i++)
        if (!marked[i]) count++;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %d primes between 3 and %d.\n", global_count, n);
        printf("Sieve of Eratosthenes took %10.6f seconds with %d processes. \n", elapsed_time, p);
    }
    MPI_Finalize();
    free(marked);
    return 0;
}


int main(int argc, char **argv) {

#if VERSION == 1
    return main_1(argc, argv);
#elif VERSION == 2
    return main_2(argc, argv);
#elif VERSION == 3
    return main_3(argc, argv);
#elif VERSION == 4
    return main_4(argc, argv);
#endif
}