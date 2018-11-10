CC=c++
CFLAG=-g -Wall -std=c++11

all: nbody generator

nbody:
	$(CC) $(CFLAG) main.cpp -o nbody

generator:
	$(CC) $(CFLAG) generator.cpp -o generator

clean:
	rm -f nbody generator *.o
