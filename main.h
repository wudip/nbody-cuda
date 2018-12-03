#ifndef BAKAJ_WUDI_CUDA_MAIN_H
#define BAKAJ_WUDI_CUDA_MAIN_H

#include <iostream>
#include <vector>

#include "vec3.h"
#include "particle.h"

std::vector<Particle> *loadParticles(std::istream &input);

std::vector<Vec3<double>> nbody(const std::vector<Particle> *particles);

Vec3<double> *nbodyBarnesHut(Particle *particles, unsigned int nOfParticles, Cell &cell);

void moveParticles(std::vector<Particle> *particles, const Vec3<double> *forces);

void printParticles(const std::vector<Particle> *particles, std::ostream &out);

double *computeParticleBoundaries(const std::vector<Particle> *particles);

#endif //BAKAJ_WUDI_CUDA_MAIN_H
