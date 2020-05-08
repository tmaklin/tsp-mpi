GCC=gcc
CCFLAGS=-O3
LDFLAGS=-O3

tsp: tsp.o
	$(GCC) -o tsp tsp.o $(LDFLAGS)

tsp.o : tsp.c
	$(GCC) -o tsp.o -c tsp.c $(CCFLAGS)

cleanall:
	rm tsp tsp.o
