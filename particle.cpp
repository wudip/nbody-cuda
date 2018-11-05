#include "vec3.cpp"

class Particle {
    private:
    public:
        Vec3<double> position;
        Vec3<double> velocity;
        double mass;
        Particle() : position(Vec3<double>()), velocity(Vec3<double>()), mass(0) { };
        Particle(double x, double y, double z, double mass) : position(Vec3<double>(x, y, z)), velocity(Vec3<double>(0, 0, 0)), mass(mass) { };
        Particle(Vec3<double> position, double mass) : position(position), velocity(Vec3<double>(0, 0, 0)), mass(mass) { };

        friend std::ostream& operator<<( std::ostream &o, const Particle &particle) {
            o << particle.position << " " << particle.mass;
            return o;
        }
};
