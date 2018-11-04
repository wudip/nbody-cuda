#include <iostream>
#include <cmath>
 
template <class T> class Vec3 {
    public:
        T x, y, z;

        Vec3() : x(0), y(0), z(0) { };
        Vec3(T x, T y, T z) {
            this->x = x;
            this->y = y;
            this->z = z;
        }

        void set(const T &x, const T &y, const T &z) {
            this->x = x;
            this->y = y;
            this->z = z;
        }
 
        void normalise() {
            T magnitude = sqrt((x * x) + (y * y) + (z * z));
            if (magnitude != 0)
            {
                x /= magnitude;
                y /= magnitude;
                z /= magnitude;
            }
        }

        void square() {
            x*=x;
            y*=y;
            z*=z;
        }

        static T dotProduct(const Vec3 &a, const Vec3 &b) {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
 
        T dotProduct(const Vec3 &vec) const {
            return x * vec.x + y * vec.y + z * vec.z;
        }

        static T getDistance(const Vec3 &a, const Vec3 &b) {
            T dx = b.x - a.x;
            T dy = b.y - a.y;
            T dz = b.z - a.z; 
            return sqrt(dx * dx + dy * dy + dz * dz);
        }

 
        Vec3 operator+(const Vec3 &vector) const {
            return Vec3<T>(x + vector.x, y + vector.y, z + vector.z);
        }
 
        void operator+=(const Vec3 &vector) {
            x += vector.x;
            y += vector.y;
            z += vector.z;
        }
 
        Vec3 operator-(const Vec3 &vector) const {
            return Vec3<T>(x - vector.x, y - vector.y, z - vector.z);
        }
 
        void operator-=(const Vec3 &vector) {
            x -= vector.x;
            y -= vector.y;
            z -= vector.z;
        }
 
        Vec3 operator*(const Vec3 &vector) const {
            return Vec3<T>(x * vector.x, y * vector.y, z * vector.z);
        }
 
        Vec3 operator*(const T &value) const {
            return Vec3<T>(x * value, y * value, z * value);
        }
 
        void operator*=(const T &value) {
            x *= value;
            y *= value;
            z *= value;
        }
 
        Vec3 operator/(const T &value) const {
            return Vec3<T>(x / value, y / value, z / value);
        }
 
        void operator/=(const T &value) {
            x /= value;
            y /= value;
            z /= value;
        }

        template <typename U>
        friend std::ostream& operator<<( std::ostream &o, const Vec3<U> &vector) {
            o << "Vec3(" << vector.x << ", " << vector.y << ", " << vector.z << ")";
            return o;
        }
};