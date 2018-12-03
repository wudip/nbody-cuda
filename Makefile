CC=nvcc
CFLAG=-g -std=c++11 -arch compute_35

all: nbody generator

nbody: main.cu cell.cu simple_cell.o cell.o particle.o vec3.o
	$(CC) $(CFLAG) main.cu cell.o simple_cell.o particle.o vec3.o -o nbody

cell.o: cell.cu simple_cell.o cell.h particle.o vec3.o
	$(CC) $(CFLAG) -dc cell.cu simple_cell.o particle.o vec3.o

simple_cell.o: simple_cell.cu simple_cell.h particle.o vec3.o
	$(CC) $(CFLAG) -dc simple_cell.cu particle.o vec3.o

particle.o: particle.cu particle.h vec3.o
	$(CC) $(CFLAG) -dc particle.cu vec3.o

vec3.o: vec3.cu vec3.h
	$(CC) $(CFLAG) -dc vec3.cu

generator:
	$(CC) $(CFLAG) generator.cpp -o generator

clean:
	rm -f nbody generator *.o
