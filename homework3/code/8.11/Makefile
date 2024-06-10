all: transpose

transpose: transpose.c
	mpicc -o transpose transpose.c 

clean:
	rm -f transpose

.PHONY: all clean