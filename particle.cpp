#include "vec3.cpp"

class Cell;

class Particle {
    private:
        Vec3<double> position;
        Vec3<double> velocity;
    public:
        double mass;
        Cell * cell;
        Particle() : position(Vec3<double>()), velocity(Vec3<double>()), mass(0) { };
        Particle(double x, double y, double z, double mass) : position(Vec3<double>(x, y, z)), velocity(Vec3<double>(0, 0, 0)), mass(mass) { };
        Particle(Vec3<double> position, double mass) : position(position), velocity(Vec3<double>(0, 0, 0)), mass(mass) { };
        Particle& operator=(const Particle & p) = default;

        const Vec3<double>& getPosition() const {
            return position;
        }

        /**
         * Updates position of the particle according to its velocity
         */
        void updatePosition() {
            position += velocity;
        }

        void accelerate(Vec3<double> acceleration) {
            velocity += acceleration;
        }

        friend std::ostream& operator<<( std::ostream &o, const Particle &particle) {
            o << particle.position << " " << particle.mass;
            return o;
        }
};
