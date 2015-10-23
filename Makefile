CC = g++
CFLAGS = -O3 -arch sm_20 g
NVCCFLAGS = -g -O3 -arch sm_20
LIBS = -lm

all: simulate_empty_box simulate_lattice_model

simulate_empty_box: empty_GCMC.cc
	$(CC) -o $@ empty_GCMC.cc -std=c++11

simulate_lattice_model: lattice_GCMC.cc
	$(CC) -o $@ lattice_GCMC.cc -std=c++11

clean: 
	rm *.o 
