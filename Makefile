CC=c++
CFLAG=-g -Wall -std=c++11

all:
	$(CC) $(CFLAG) main.cpp -o nbody

clean:
	rm -f nbody