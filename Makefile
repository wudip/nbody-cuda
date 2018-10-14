LIBRARIES = -L/usr/local/lib -lglfw3 -lGL -lGLU -lX11 -lXxf86vm -lXcursor -lXinerama -lXrandr -pthread -lXi -ldl
CONFIG = -I/usr/local/include -std=c++14

all: nbody

nbody: nbody.o
	g++ $(CONFIG) nbody.o $(LIBRARIES) -o nbody

nbody.o: src/main.cpp
	g++ $(CONFIG) src/main.cpp $(LIBRARIES) -c -o nbody.o

clean:
	rm -f nbody *.o
