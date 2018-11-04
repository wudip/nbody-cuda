#include <iostream>

#include "vec3.cpp"

using namespace std;

int main(int argc, char** argv) {
    Vec3<double> A(1.5, 2.8, 6.2);
    Vec3<double> B(8.5, 4.5, 1);

    cout << "A: " << A << endl;
    cout << "B: " << B << endl;
    cout << "A+B: " << A+B << endl;
    cout << "A*B: " << A*B << endl;
    cout << "A*0.5: " << A*0.5 << endl;
    A.normalise();
    cout << "A.normalise: " << A << endl;
    B.square();
    cout << "B.square: " << B << endl;

    return 0;
}