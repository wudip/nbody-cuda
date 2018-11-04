CC=c++
CFLAG=-g -Wall -std=c++11

all:
	$(CC) $(CFLAG) vec3.cpp particle.cpp test.cpp -o app

clean:
	rm -f app