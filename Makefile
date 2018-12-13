CC=c++
CFLAG=-g -Wall -std=c++11

all: nbody generator

nbody: main.cpp cell.cpp simple_cell.o cell.o particle.o vec3.o
	$(CC) $(CFLAG) main.cpp cell.o simple_cell.o particle.o vec3.o -o nbody

cell.o: cell.cpp simple_cell.o cell.h particle.o vec3.o
	$(CC) $(CFLAG) cell.cpp simple_cell.o particle.o vec3.o -c -o cell.o

simple_cell.o: simple_cell.cpp simple_cell.h particle.o vec3.o
	$(CC) $(CFLAG) simple_cell.cpp particle.o vec3.o -c -o simple_cell.o

particle.o: particle.cpp particle.h vec3.o
	$(CC) $(CFLAG) particle.cpp vec3.o -c -o particle.o

vec3.o: vec3.cpp vec3.h
	$(CC) $(CFLAG) vec3.cpp -c -o vec3.o

generator: generator.cpp
	$(CC) $(CFLAG) generator.cpp -o generator

clean:
	rm -f nbody generator *.o