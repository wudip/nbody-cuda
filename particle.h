#ifndef BAKAJ_WUDI_CUDA_PARTICLE_H
#define BAKAJ_WUDI_CUDA_PARTICLE_H

#include "vec3.h"

class Cell;

class Particle {
private:
    Vec3<double> position;
    Vec3<double> velocity;
public:
    double mass;
    Cell *cell;

    Particle();

    Particle(double x, double y, double z, double mass);

    Particle(Vec3<double> position, double mass);

    Particle &operator=(const Particle &p) = default;

    #pragma acc routine seq
    const Vec3<double> &getPosition() const;

    /**
     * Updates position of the particle according to its velocity
     */
    void updatePosition();

    void accelerate(Vec3<double> acceleration);

    friend std::ostream &operator<<(std::ostream &o, const Particle &particle);
};


#endif //BAKAJ_WUDI_CUDA_PARTICLE_H
