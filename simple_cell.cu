#include <cstdio>
#include <cmath>
#include "simple_cell.h"

#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11

#define IS_NOT_EMPTY(index) ((unsigned int) -1) != index

SimpleCell::SimpleCell() {}

SimpleCell::SimpleCell(unsigned int offset, unsigned int particle, const Particle &center, unsigned int parent,
                       unsigned int *sub) :
        offset(offset), particle(particle), center(center), parent(parent) {
    for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
        subtree[i] = sub[i];
    }
}

__device__ const SimpleCell *SimpleCell::getCell(unsigned int position) const {
    int diff = position - offset;
    return this + diff;
}

__device__ void SimpleCell::getForceSiblings(const Particle &refParticle, Vec3<double> &forces) const {
    for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
        unsigned int sibIndex = getCell(parent)->subtree[i];
        if (sibIndex == offset) continue;
        const SimpleCell *sibling = getCell(sibIndex);
        const Particle &sibPart = sibling->center;
        addToForces(forces, refParticle, sibPart);
    }
}

__device__ Vec3<double> SimpleCell::getForce(const Particle *particles) const {
    Vec3<double> force(0, 0, 0);
    const SimpleCell *c = this;
    while (IS_NOT_EMPTY(c->parent)) {
        getForceSiblings(particles[c->particle], force);
        c = getCell(c->parent);
    }
    return force;
}

__device__ void addToForces(Vec3<double> &forces, const Particle &particle, const Particle &sibPart) {
    Vec3<double> diff = particle.getPosition() - sibPart.getPosition();
    double bottom = pow(diff.sqrSize() + SOFTENING_FACTOR_SQR, 1.5);
    double massTotal = GRAVITATION_CONSTANT * particle.mass * sibPart.mass;
    forces += diff * massTotal / bottom;
}
