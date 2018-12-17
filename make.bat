nvcc -dc -std=c++11 -arch compute_35 vec3.cu particle.cu cell.cu simple_cell.cu main.cu && nvcc -lcuda -lcudart -std=c++11 -arch sm_35 vec3.obj particle.obj cell.obj simple_cell.obj main.obj -o nbody
