GCC=gcc
CCFLAGS=-O3
LDFLAGS=-lm -O3
MPICC=mpicc

all: tsp-mpi

tsp-mpi: tsp-mpi.o
	$(MPICC) -o tsp-mpi tsp-mpi.o $(LDFLAGS)

tsp-mpi.o : tsp-mpi.c
	$(MPICC) -o tsp-mpi.o -c tsp-mpi.c $(CCFLAGS)

cleanall:
	rm tsp-mpi tsp-mpi.o
