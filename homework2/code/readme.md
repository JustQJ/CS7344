
## 4.8
代码运行命令：
```shell
mpicc -O -o consecutive_prime_count consecutive_prime_count.c
mpirun -np 2 ./consecutive_prime_count #-np 可以选择不同的进程数量
```

## 5.7
代码运行命令：
```shell
mpicc -DVERSION=3 -o sieve_prime sieve_prime.c -lm #-DVERSION 可以是1，2，3，4；1是基本版本，2去掉了偶数，3去掉了通信，4利用了cache block，提高hit ratio
mpirun -np 6 ./sieve_prime 100000000 #-np 可以选择不同的进程数量
```

##6.9
代码运行命令：
```shell
mpicc -o mpi_reduce mpi_reduce.c
mpirun -np 256 ./mpi_reduce #-np 可以选择不同的进程数量
```