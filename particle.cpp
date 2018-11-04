#include "vec3.cpp"

class Particle {
    private:
        Vec3<double> direction;
        double mass;
    public:
        Particle() : direction(Vec3<double>()), mass(0) { };
        Particle(Vec3<double> direction, double mass) : direction(direction), mass(mass) { };
};