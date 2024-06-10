#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

//judge if the interval has root
int has_root(double left, double right){
    // f(x) = -2 +sin(x^1) + sin(x^2) + .. + sin(x^100)

    // double f_left = -2;
    // double f_right = -2;
    // double current_left = 1;
    // double current_right = 1;
    // for(int i=1; i<=100; i++){
    //     current_left *= left;
    //     current_right *= right;
    //     f_left += sin(current_left);
    //     f_right += sin(current_right);       
    // }
    // if (f_left*f_right<=0){
    //     return 1;
    // }
    // return 0;

    double f_left = -2;
    double f_right = -2;
    for(int i=1; i<=100; i++){
        f_left += sin(pow(left, i));
        f_right += sin(pow(right, i));
    }
    if (f_left*f_right<=0){
        return 1;
    }
    return 0;
}


int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int rank, work_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &work_size);

    if(work_size==1){
        //only one process, no need to use manager-worker style, use serial version
        printf("Only one process, use serial version with binary search\n");
        //binary search
        //recording time
        double elapsed_time = -MPI_Wtime();
        double interval[2] = {0, 1};
        double eps = 1e-11;
        while(interval[1]-interval[0]>eps){
            double mid = (interval[0]+interval[1])/2;
            if(has_root(interval[0], mid)){
                interval[1] = mid;
            }else{
                interval[0] = mid;
            }
        }
        printf("The root interval: [ %0.12f, %0.12f], the root is %0.12f\n", interval[0], interval[1], (interval[0]+interval[1])/2);
        elapsed_time += MPI_Wtime();
        printf("Elapsed time: %0.6f\n", elapsed_time);
        MPI_Finalize();
        return 0;
    }

    if(work_size<3){
        printf("At least 3 processes are needed for manager-worker style\n");
        MPI_Finalize();
        return 0;
    }

    //recording time
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = -MPI_Wtime();
    

    //1 manager and work_size-1 workers
    if(rank==0){//manager process
        double interval[2] = {0, 1};
        double eps = 1e-11;
        while(interval[1]-interval[0]>eps){ //when the interval is small enough, stop
            double interval_size = (interval[1]-interval[0])/(work_size-1); //split the interval into worker number sub-intervals
            for(int i=0; i<work_size-2; i++){ //assgin the sub-interval to worker processes
                double sub_interval[2] = {interval[0]+interval_size*i, interval[0]+interval_size*(i+1)};
                MPI_Send(sub_interval, 2, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
            }
            double sub_interval[2] = {interval[0]+interval_size*(work_size-2), interval[1]}; //ensure the last worker process the rest, avoid precision problem
            MPI_Send(sub_interval, 2, MPI_DOUBLE, work_size-1, 0, MPI_COMM_WORLD);

            //recevice the interval that has root, only the process that has root will send back the interval
            MPI_Recv(interval, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        printf("The root interval: [ %0.12f, %0.12f], the root is %0.12f\n", interval[0], interval[1], (interval[0]+interval[1])/2);
        elapsed_time += MPI_Wtime();
        printf("Elapsed time: %0.6f\n", elapsed_time);
        
        //send stop signal to all worker processes
        MPI_Abort(MPI_COMM_WORLD, 0);

    }else{
        //worker process
        while (1)
        {
            double sub_interval[2];
            MPI_Recv(sub_interval, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive the interval from manager process
            if(has_root(sub_interval[0], sub_interval[1])){ //judge if the interval has root, if has, send back to manager process
                MPI_Send(sub_interval, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Finalize();

    
    
    return 0;
}


