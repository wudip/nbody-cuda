#include "vec3.cpp"

class Particle {
    private:
    public:
        Vec3<double> direction;
        double mass;
        Particle() : direction(Vec3<double>()), mass(0) { };
        Particle(double x, double y, double z, double mass) : direction(Vec3<double>(x, y, z)), mass(mass) { };
        Particle(Vec3<double> direction, double mass) : direction(direction), mass(mass) { };

        friend std::ostream& operator<<( std::ostream &o, const Particle &particle) {
            o << particle.direction << " " << particle.mass;
            return o;
        }
};