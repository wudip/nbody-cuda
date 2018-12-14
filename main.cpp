#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>
#include <utility>

#include "main.h"
#include "cell.h"

// Softening factor squared
#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11
#define WINDOWS_SUCKS 1
#define NUM_REPEATS 5

using namespace std;

int main(int argc, char **argv) {
    // Load input
    pair<Particle*, unsigned int> input;
    if (WINDOWS_SUCKS && argc > 1) {
        std::ifstream ifs;
        ifs.open(argv[1], ifstream::in);
        input = loadParticles(ifs);
        ifs.close();
    } else {
        input = loadParticles(cin);
    }
    Particle* particles = input.first;
    unsigned int size = input.second;
    if (particles == nullptr) return 1;

    // Actual algorithm
    clock_t clk_start = clock();
    for (int i = 0; i < NUM_REPEATS; ++i) {
        // Create octree
        double *particleBoundaries = computeParticleBoundaries(particles, size);
        Cell octree(particleBoundaries, particleBoundaries + 3);
        delete[] particleBoundaries;
        for (unsigned int partIndex = 0; partIndex < size; ++partIndex) {
            octree.add(particles + partIndex);
        }
        octree.updateCenter();
        // Solve N-body problem
        nbodyBarnesHut(particles, size, octree);
    }
    clock_t clk_end = clock();
    //printParticles(particles, size, cout);
    delete[] particles;
    long time = (clk_end - clk_start) / (CLOCKS_PER_SEC / 1000);
    cout << "Time: " << time << " ms" << endl;
    return 0;
}

pair<Particle *, unsigned int> loadParticles(istream &input) {
    // Loading vector
    vector<Particle> particles;
    double x;
    double y;
    double z;
    double mass;
    while (input >> x && input >> y && input >> z && input >> mass) {
        particles.push_back(Particle(x, y, z, mass));
    }

    // Move particles from vector to array
    unsigned int size = (unsigned int) particles.size();
    Particle* particleArr = size == 0 ? nullptr : new Particle[size];
    for (auto pit = particles.begin(); pit < particles.end(); ++pit) {
        particleArr[pit - particles.begin()] = * pit;
    }

    return {particleArr, size};
}

vector<Vec3<double>> nbody(const vector<Particle> *particles) {
    vector<Vec3<double>> forces;
    forces.reserve(particles->size());
    for (auto it1 = particles->begin(); it1 < particles->end(); ++it1) {
        Vec3<double> res(0, 0, 0);
        for (auto it2 = particles->begin(); it2 < particles->end(); ++it2) {
            if (it1 == it2) continue;
            Vec3<double> diff = it1->getPosition() - it2->getPosition();
            double bottom = pow(diff.sqrSize() + SOFTENING_FACTOR_SQR, 1.5);
            double massTotal = GRAVITATION_CONSTANT * it1->mass * it2->mass;
            res += diff * massTotal / bottom;
        }
        forces.push_back(res);
    }
    return forces;
}

void nbodyBarnesHut(Particle *particles, unsigned int nOfParticles, Cell &cell) {
    unsigned int nOfCells;
    pair<SimpleCell *, unsigned int *> serialized = cell.serialize(particles, nOfCells);
    SimpleCell *flatTree = serialized.first;
    unsigned int *partPositions = serialized.second;
    Vec3<double>* forces = new Vec3<double>[nOfParticles];

    #pragma acc parallel loop copy(partPositions[0:nOfParticles], flatTree[0:nOfCells], forces[0:nOfParticles], particles[0:nOfParticles])
    for(unsigned int index = 0; index < nOfParticles; ++index) {
        SimpleCell *particleCell = flatTree + partPositions[index];
        Vec3<double> force = particleCell->getForce(particles) * 10000;
        Vec3<double> acceleration = force / particles[index].mass;
        forces[index] = acceleration;
    }
    moveParticles(particles, forces, nOfParticles);

    delete[] flatTree;
    delete[] partPositions;
}

void moveParticles(Particle *particles, const Vec3<double> *forces, unsigned int n) {
    auto particleIt = particles;
    auto forceIt = forces;
    for(unsigned int i = 0; i < n; i++) {
        Vec3<double> acceleration = *forceIt / particleIt->mass;
        particleIt->accelerate(acceleration);
        particleIt->updatePosition();
        ++particleIt;
        ++forceIt;
    }
}

/**
 * Computes coordinates of smallest possible box containing all the particles
 * @return array of minimum point and maximum point: [min_x, min_y, min_z, max_x, max_y, max_z]
 */
double *computeParticleBoundaries(const Particle* particles, unsigned int size) {
    double *result = new double[6];
    for (int i = 0; i < 3; ++i) {
        result[i] = INFINITY;
        result[3 + i] = -INFINITY;
    }
    for (const Particle* it = particles; it < particles + size; ++it) {
        Vec3<double> pos = it->getPosition();
        for (int dim = 0; dim < 3; ++dim) {
            double coord = pos.getDim(dim);
            if (coord < result[dim]) {
                result[dim] = coord;
            }
            if (coord > result[3 + dim]) {
                result[3 + dim] = coord;
            }
        }
    }
    return result;
}

void printParticles(const Particle *particles, unsigned int size, ostream &out) {
    for (unsigned int i = 0; i < size; ++i) {
        out << particles[i] << endl;
    }
}

