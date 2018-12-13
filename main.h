#ifndef BAKAJ_WUDI_CUDA_MAIN_H
#define BAKAJ_WUDI_CUDA_MAIN_H

#include <iostream>
#include <vector>
#include <utility>

#include "vec3.h"
#include "particle.h"

std::pair<Particle*, unsigned int> loadParticles(istream &input);

std::vector<Vec3<double>> nbody(const std::vector<Particle> *particles);

void nbodyBarnesHut(Particle *particles, unsigned int nOfParticles, Cell &cell);

__global__ void moveParticles(Particle *particles, const Vec3<double> *forces, unsigned int nOfParticles);

void printParticles(const Particle *particles, unsigned int size, std::ostream &out);

double *computeParticleBoundaries(const Particle *particles, unsigned int size);

#endif //BAKAJ_WUDI_CUDA_MAIN_H
