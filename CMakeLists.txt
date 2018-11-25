cmake_minimum_required(VERSION 3.12)
project(nbody_bakaj_wudi)

set(CMAKE_CXX_STANDARD 14)

add_executable(nbody main.cpp particle.cpp simple_cell.cpp cell.cpp vec3.cpp)
add_executable(generator generator.cpp)
