#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>


//function to merge two arrays
void merge(int *arr, int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
    int *L = (int *)malloc(n1 * sizeof(int));
    int *R = (int *)malloc(n2 * sizeof(int));
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j])
            arr[k++] = L[i++];
        else
            arr[k++] = R[j++];
    }
    while (i < n1)
        arr[k++] = L[i++];
    while (j < n2)
        arr[k++] = R[j++];
    free(L);
    free(R);
}

//function to do merge sort
void mergeSort(int *arr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 1000000;
    if (argc > 1) {
        n = atoi(argv[1]);
    }

    if(size&(size-1)!=0){
            printf("Number of processes should be power of 2\n");
            MPI_Finalize();
            return 0;
    }
    if(rank==0){

        printf("n = %d and p = %d \n", n, size);
    }

    //record time
    double elapsed_time = -MPI_Wtime();

    
    int *local_arr;
    int num_elements = 0;


    if(rank==0){
        //generate array with n random numbers
        srand(time(NULL));
        local_arr = (int *)malloc(n * sizeof(int));
        for (int i = 0; i < n; i++) {
            local_arr[i] = rand()%n;
        }
        
        //split the array into p parts
        //each part has n/p+1 elements or n/p elements
        //send each part to a different process
        int num_per_process = n/size;
        int remander = n%size;

        int start = num_per_process;
        if(remander!=0){
            start++;
        }
        for(int i=1; i<size; i++){
            int end = start+num_per_process-1;
            if(i<remander){
                end++;
            }
            int len = end-start+1;
            MPI_Send(local_arr+start, len, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        //sort its own part
        mergeSort(local_arr, 0, start-1);

        num_elements = start;

    }else{
        //receive the part from process 0
        num_elements = n/size;
        int remander = n%size;

        if(rank<=remander){
            num_elements++;
        }
        local_arr = (int *)malloc(num_elements * sizeof(int));
        MPI_Recv(local_arr, num_elements, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //sort the part
        mergeSort(local_arr, 0, num_elements-1);
    }

    //merge the sorted parts with reduce operation
   
    for(int step=1; step<size; step*=2){
        if(rank%(2*step)==0){ //responsible for merging
            if(rank+step<size){
                int num_elements2;
                MPI_Status status;
                //probe to get the number of elements in the part to be received
                //receive the data coming from the process rank+step
                MPI_Probe(rank+step, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_INT, &num_elements2);
                int *temp = (int *)malloc((num_elements+num_elements2) * sizeof(int));
                //store the received data in last part of temp
                MPI_Recv(temp+num_elements, num_elements2, MPI_INT, rank+step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                //merge the two parts
                int i=0, j=num_elements, k=0;
                while(i<num_elements && j<num_elements+num_elements2){
                    if(local_arr[i]<=temp[j]){
                        temp[k++] = local_arr[i++];
                    }else{
                        temp[k++] = temp[j++];
                    }
                }
                while(i<num_elements){
                    temp[k++] = local_arr[i++];
                }

                free(local_arr);
                local_arr = temp;
                num_elements += num_elements2;
            }
        }else{
            //send the part to the process rank-step
            //exclude the sending of the process last step
            if(rank-step>=0){
                MPI_Send(local_arr, num_elements, MPI_INT, rank-step, 0, MPI_COMM_WORLD);
                break;
            }
        }
    }

    //print the sorted array
    elapsed_time += MPI_Wtime();
    // if(rank==0){
    //     for(int i=0; i<n; i++){
    //         printf("%d ", local_arr[i]);
    //     }
    //     printf("\n");
    // }
    if(rank==0){
        printf("Time taken: %f\n", elapsed_time);
    }
    if(local_arr!=NULL){
        free(local_arr);
    }
    MPI_Finalize();
    return 0;
    


}