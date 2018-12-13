#ifndef BAKAJ_WUDI_CUDA_VEC3_H
#define BAKAJ_WUDI_CUDA_VEC3_H

#include <iostream>

template<class T>
class Vec3 {
public:
    T x, y, z;

    #pragma acc routine seq
    Vec3();

    #pragma acc routine seq
    Vec3(T x, T y, T z);

    void set(const T &x, const T &y, const T &z);

    void normalise();

    void square();

    static T dotProduct(const Vec3 &a, const Vec3 &b);

    T dotProduct(const Vec3 &vec) const;

    #pragma acc routine seq
    static T getDistance(const Vec3 &a, const Vec3 &b);

    #pragma acc routine seq
    Vec3 operator+(const Vec3 &vector) const;

    #pragma acc routine seq
    void operator+=(const Vec3 &vector);

    #pragma acc routine seq
    Vec3 operator-(const Vec3 &vector) const;

    #pragma acc routine seq
    void operator-=(const Vec3 &vector);

    #pragma acc routine seq
    Vec3 operator*(const Vec3 &vector) const;

    #pragma acc routine seq
    Vec3 operator*(const T &value) const;

    #pragma acc routine seq
    void operator*=(const T &value);

    // Vec3 operator/(const T &value) const;

    #pragma acc routine seq
    Vec3 operator/(const T &value) const;

    #pragma acc routine seq
    void operator/=(const T &value);

    #pragma acc routine seq
    T sqrSize();

    double getDim(int dimension) const;

    void setDim(int dimension, const double value);

    template<typename U>
    friend std::ostream &operator<<(std::ostream &o, const Vec3<U> &vector) {
        o << vector.x << " " << vector.y << " " << vector.z;
        return o;
    }
};

#endif //BAKAJ_WUDI_CUDA_VEC3_H
