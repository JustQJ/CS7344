all: find_root

find_root: find_root.c
	mpicc -o find_root find_root.c -lm

clean:
	rm -f find_root

.PHONY: all clean