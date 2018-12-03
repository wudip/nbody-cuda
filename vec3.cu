
#include <cmath>
#include "vec3.h"

template<typename T>
__host__ __device__ Vec3<T>::Vec3() : x(0), y(0), z(0) {};

template<typename T>
__host__ __device__ Vec3<T>::Vec3(T x, T y, T z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

template<typename T>
void Vec3<T>::set(const T &x, const T &y, const T &z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

template<typename T>
void Vec3<T>::normalise() {
    T magnitude = sqrt((x * x) + (y * y) + (z * z));
    if (magnitude != 0) {
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
    }
}

template<typename T>
void Vec3<T>::square() {
    x *= x;
    y *= y;
    z *= z;
}

template<typename T>
T Vec3<T>::dotProduct(const Vec3<T> &a, const Vec3<T> &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<typename T>
T Vec3<T>::dotProduct(const Vec3<T> &vec) const {
    return x * vec.x + y * vec.y + z * vec.z;
}

template<typename T>
__host__ __device__ T Vec3<T>::getDistance(const Vec3<T> &a, const Vec3<T> &b) {
    T dx = b.x - a.x;
    T dy = b.y - a.y;
    T dz = b.z - a.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}


template<typename T>
__host__ __device__ Vec3<T> Vec3<T>::operator+(const Vec3<T> &vector) const {
    return Vec3<T>(x + vector.x, y + vector.y, z + vector.z);
}

template<typename T>
__host__ __device__ void Vec3<T>::operator+=(const Vec3<T> &vector) {
    x += vector.x;
    y += vector.y;
    z += vector.z;
}

template<typename T>
__host__ __device__ Vec3<T> Vec3<T>::operator-(const Vec3<T> &vector) const {
    return Vec3<T>(x - vector.x, y - vector.y, z - vector.z);
}

template<typename T>
__host__ __device__ void Vec3<T>::operator-=(const Vec3<T> &vector) {
    x -= vector.x;
    y -= vector.y;
    z -= vector.z;
}

template<typename T>
__host__ __device__ Vec3<T> Vec3<T>::operator*(const Vec3<T> &vector) const {
    return Vec3<T>(x * vector.x, y * vector.y, z * vector.z);
}

template<typename T>
__host__ __device__ Vec3<T> Vec3<T>::operator*(const T &value) const {
    return Vec3<T>(x * value, y * value, z * value);
}

template<typename T>
__host__ __device__ void Vec3<T>::operator*=(const T &value) {
    x *= value;
    y *= value;
    z *= value;
}

// template<typename T>
// Vec3<T> Vec3<T>::operator/(const T &value) const {
//     return Vec3<T>(x / value, y / value, z / value);
// }

template<typename T>
__host__ __device__ Vec3<T> Vec3<T>::operator/(const T &value) const {
    return Vec3<T>(x / value, y / value, z / value);
}

template<typename T>
__host__ __device__ void Vec3<T>::operator/=(const T &value) {
    x /= value;
    y /= value;
    z /= value;
}

template<typename T>
__host__ __device__ T Vec3<T>::sqrSize() {
    return x * x + y * y + z * z;
}

template<typename T>
double Vec3<T>::getDim(int dimension) const {
    if (dimension == 0) return x;
    if (dimension == 1) return y;
    if (dimension == 2) return z;
    throw "Dimension out of range";
}

template<typename T>
void Vec3<T>::setDim(int dimension, const double value) {
    if (dimension == 0) x = value;
    else if (dimension == 1) y = value;
    else if (dimension == 2) z = value;
    else throw "Dimension out of range";
}

template
class Vec3<double>;
