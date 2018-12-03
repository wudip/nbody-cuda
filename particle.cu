#include "particle.h"

class Cell;

Particle::Particle() :
        position(Vec3<double>()), velocity(Vec3<double>()), mass(0) {};

Particle::Particle(double x, double y, double z, double mass) :
        position(Vec3<double>(x, y, z)), velocity(Vec3<double>(0, 0, 0)), mass(mass) {};

Particle::Particle(Vec3<double> position, double mass) :
        position(position), velocity(Vec3<double>(0, 0, 0)), mass(mass) {};

__host__ __device__ const Vec3<double> &Particle::getPosition() const {
    return position;
}

/**
 * Updates position of the particle according to its velocity
 */
__host__ __device__ void Particle::updatePosition() {
    position += velocity;
}

__host__ __device__ void Particle::accelerate(Vec3<double> acceleration) {
    velocity += acceleration;
}

std::ostream &operator<<(std::ostream &o, const Particle &particle) {
    o << particle.position << " " << particle.mass;
    return o;
}
