#include <iostream>
#include <random>

#define NUMBER_OF_DIMENSIONS 3

using namespace std;

struct arguments {
    arguments() {
        containsError = false;
        shouldContinue = true;
    }

    double positionVariance;
    double massAvg;
    double massVariance;
    int particleSum;
    bool containsError;
    bool shouldContinue;
};

arguments *readArguments(int argc, char **argv) {
    arguments *a = new arguments();
    a->particleSum = 5; // TODO
    a->positionVariance = 0.3; // TODO
    a->massAvg = 5; // TODO
    a->massVariance = 6; // TODO
    return a;
}

void generateParticle(
        normal_distribution<double> &positionDistribution,
        normal_distribution<double> &massDistribution,
        default_random_engine &randomGenerator
) {
    for (int dimension = 0; dimension < NUMBER_OF_DIMENSIONS; ++dimension) {
        double coord = positionDistribution(randomGenerator);
        cout << coord << " ";
    }
    double mass = massDistribution(randomGenerator);
    cout << mass << endl;
}

int main(int argc, char **argv) {
    arguments *a = readArguments(argc, argv);
    if (a->containsError) {
        return 1;
    }
    if (!a->shouldContinue) {
        return 0;
    }
    default_random_engine randomGenerator;
    normal_distribution<double> positionDistribution(0, a->positionVariance);
    normal_distribution<double> massDistribution(a->massAvg, a->massVariance);
    for (int i = 0; i < a->particleSum; ++i) {
        generateParticle(positionDistribution, massDistribution, randomGenerator);
    }
    return 0;
}
