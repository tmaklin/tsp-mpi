GCC=gcc
CCFLAGS=
LDFLAGS=-lm
MPICC=mpicc

all: tsp tsp-mpi

tsp: tsp.o
	$(GCC) -o tsp tsp.o $(LDFLAGS)

tsp.o : tsp.c
	$(GCC) -o tsp.o -c tsp.c $(CCFLAGS)

tsp-mpi: tsp-mpi.o
	$(MPICC) -o tsp-mpi tsp-mpi.o $(LDFLAGS)

tsp-mpi.o : tsp-mpi.c
	$(MPICC) -o tsp-mpi.o -c tsp-mpi.c $(CCFLAGS)

cleanall:
	rm tsp tsp.o tsp-mpi tsp-mpi.o
