all: 1_adi_p1 3_floyd_p3 3_floyd_p3_1

1_adi_p1: 1_adi_p1.c
	mpicc -o 1_adi_p1 1_adi_p1.c 

3_floyd_p3: 3_floyd_p3.c
	mpicc -o 3_floyd_p3 3_floyd_p3.c

3_floyd_p3_1: 3_floyd_p3_1.c
	mpicc -o 3_floyd_p3_1 3_floyd_p3_1.c -lm

clean:
	rm -f 1_adi_p1 3_floyd_p3 3_floyd_p3_1

.PHONY: all clean