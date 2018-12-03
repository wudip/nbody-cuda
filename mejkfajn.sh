nvcc -dc -std=c++11 -arch compute_35 *.cu && nvcc -lcuda -lcudart -std=c++11 -arch sm_35 *.o -o nbody
