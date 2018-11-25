#include "particle.h"

class Cell;

Particle::Particle():
    position(Vec3<double>()), velocity(Vec3<double>()), mass(0), cellIndex((unsigned int) -1)
{};

Particle::Particle(double x, double y, double z, double mass):
    position(Vec3<double>(x, y, z)), velocity(Vec3<double>(0, 0, 0)), mass(mass), cellIndex((unsigned int) -1)
{};

Particle::Particle(Vec3<double> position, double mass):
    position(position), velocity(Vec3<double>(0, 0, 0)), mass(mass), cellIndex((unsigned int) -1)
{};

const Vec3<double>& Particle::getPosition() const {
    return position;
}

/**
 * Updates position of the particle according to its velocity
 */
void Particle::updatePosition() {
    position += velocity;
}

void Particle::accelerate(Vec3<double> acceleration) {
    velocity += acceleration;
}

std::ostream& operator<<( std::ostream &o, const Particle &particle) {
    o << particle.position << " " << particle.mass;
    return o;
}
