CC = g++
CFLAGS = -O3 -arch sm_20 g
NVCCFLAGS = -g -O3 -arch sm_20
LIBS = -lm

all: simulate

simulate: empty_GCMC.cc
	$(CC) -o $@ empty_GCMC.cc -std=c++11

clean: 
	rm *.o 
