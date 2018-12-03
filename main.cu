#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstdio>

#include "main.h"
#include "cell.h"

// Softening factor squared
#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11
#define WINDOWS_SUCKS 1

using namespace std;

int main(int argc, char **argv) {
    vector<Particle> *particles;
    if (WINDOWS_SUCKS && argc > 1) {
        std::ifstream ifs;
        ifs.open(argv[1], ifstream::in);
        particles = loadParticles(ifs);
        ifs.close();
    } else {
        particles = loadParticles(cin);
    }
    clock_t clk_start = clock();
    for (int i = 0; i < 1; ++i) {
        // Create octree
        double *particleBoundaries = computeParticleBoundaries(particles);

        // move particles from vector to array
        unsigned int size = (unsigned int) particles->size();
        Particle *particleArr = new Particle[size];
        for (auto pit = particles->begin(); pit < particles->end(); ++pit) {
            particleArr[pit - particles->begin()] = *pit;
        }

        Cell octree(particleBoundaries, particleBoundaries + 3);
        delete[] particleBoundaries;
        for (unsigned int partIndex = 0; partIndex < size; ++partIndex) {
            octree.add(particleArr + partIndex);
        }
        octree.updateCenter();


        Vec3<double> *forces = nbodyBarnesHut(particleArr, size, octree);

        // delete new particles
        for (auto pit = particles->begin(); pit < particles->end(); ++pit) {
            *pit = particleArr[pit - particles->begin()];
        }
        delete[] particleArr;

        //vector<Vec3<double>> forces = nbody(particles);
        moveParticles(particles, forces);
        delete forces;
    }
    clock_t clk_end = clock();
    cout << "Time: " << (clk_end - clk_start) << " ms" << endl;
    printParticles(particles, cout);
    delete particles;
    return 0;
}

vector<Particle> *loadParticles(istream &input) {
    vector<Particle> *particles = new vector<Particle>();
    double x;
    double y;
    double z;
    double mass;
    while (input >> x && input >> y && input >> z && input >> mass) {
        particles->push_back(Particle(x, y, z, mass));
    }
    return particles;
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

__global__ void nbodyBarnesHutCuda(
        SimpleCell * cells,
        unsigned int * partPositions,
        Particle * particles,
        unsigned int nOfParticles,
        unsigned int offset)
{
    unsigned int index = offset + threadIdx.x + blockIdx.x * blockDim.x;
    printf("Jezinka\n");
    if (index >= nOfParticles) return;
    printf("Jelen c. %d\n", partPositions[index]);
    SimpleCell *particleCell = cells + partPositions[index];
    printf("Smolicek pacholicek\n");
    Vec3<double> force = particleCell->getForce(particles);
    printf("Force: %lf %lf %lf\n", force.x, force.y, force.z);
    Vec3<double> acceleration = force / particles[index].mass;
    printf("Acc: %lf\n", acceleration);
    particles[index].accelerate(acceleration);
    particles[index].updatePosition();
}

Vec3<double> *nbodyBarnesHut(Particle *particles, unsigned int nOfParticles, Cell &cell) {
    pair<SimpleCell *, unsigned int *> serialized = cell.serialize(particles);
    SimpleCell *flatTree = serialized.first;
    unsigned int *partPositions = serialized.second;

    SimpleCell *flatTreeCuda;
    unsigned int * partPositionsCuda;
    Particle * particlesCuda;

    cudaMalloc((void**)&flatTreeCuda, (nOfParticles) * sizeof(SimpleCell));
    cudaMalloc((void**)&partPositionsCuda, (nOfParticles) * sizeof(unsigned int));
    cudaMalloc((void**)&particlesCuda, (nOfParticles) * sizeof(Particle));

    cudaMemcpy(flatTreeCuda, flatTree, sizeof(SimpleCell)*(nOfParticles), cudaMemcpyHostToDevice);
    cudaMemcpy(partPositionsCuda, partPositions, sizeof(unsigned int)*(nOfParticles), cudaMemcpyHostToDevice);
    cudaMemcpy(particlesCuda, particles, sizeof(Particle)*(nOfParticles), cudaMemcpyHostToDevice);


    // for (unsigned int index = 0; index < nOfParticles; ++index) {
        nbodyBarnesHutCuda <<<nOfParticles / 1024, 1024>>>(flatTreeCuda, partPositions, particles, nOfParticles, 0);
    // }

    cudaDeviceSynchronize();

    cudaMemcpy(particles, particlesCuda, sizeof(Particle)*(nOfParticles), cudaMemcpyDeviceToHost);

    cudaFree(flatTreeCuda);
    cudaFree(partPositionsCuda);
    cudaFree(particlesCuda);

    delete[] flatTree;
    delete[] partPositions;

    Vec3<double> *forces = new Vec3<double>[nOfParticles];
    for (unsigned int index = 0; index < nOfParticles; ++index) {
        forces[index] = Vec3<double>(0, 0, 0);
    }
    return forces;
}

void moveParticles(vector<Particle> *particles, const Vec3<double> *forces) {
    auto particleIt = particles->begin();
    auto forceIt = forces;
    while (particleIt < particles->end()) {
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
double *computeParticleBoundaries(const vector<Particle> *particles) {
    double *result = new double[6];
    for (int i = 0; i < 3; ++i) {
        result[i] = INFINITY;
        result[3 + i] = -INFINITY;
    }
    for (auto it = particles->begin(); it < particles->end(); ++it) {
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

void printParticles(const vector<Particle> *particles, ostream &out) {
    for (auto it = particles->begin(); it < particles->end(); ++it) {
        cout << *it << endl;
    }
}
