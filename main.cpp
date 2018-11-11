#include <iostream>
#include <vector>
#include <math.h>
#include <ctime>

#include "cell.cpp"

// Softening factor squared
#define SOFTENING_FACTOR_SQR 0.5
// #define GRAVITATION_CONSTANT 6.67408e-11
#define GRAVITATION_CONSTANT 6.67300E-11

using namespace std;

vector<Particle> * loadParticles(istream& input);
vector<Vec3<double>> nbody(const vector<Particle> * particles);
void moveParticles(vector<Particle> * particles, const vector<Vec3<double>>& forces);
void printParticles(const vector<Particle> * particles, ostream& out);

int main(int argc, char** argv) {
  vector<Particle> * particles = loadParticles(cin);
  clock_t start = clock();
  for(int i = 0; i < 1; ++i) {
    double min[3] = {-5., -5., -5.};
    double max[3] = {20., 20., 20.};
    Cell octree(min, max);
    for(auto it = particles->begin(); it < particles->end(); ++it) {
        octree.add(*it);
    }
    vector<Vec3<double>> forces = nbody(particles);
    moveParticles(particles, forces);
  }
  clock_t end = clock();
  cout << "Time: " << (end - start) << " ms" << endl;
  printParticles(particles, cout);
  delete particles;
  return 0;
}

vector<Particle> * loadParticles(istream& input) {
  vector<Particle> * particles = new vector<Particle>();
  double x;
  double y;
  double z;
  double mass;
  while(input >> x && input >> y && input >> z && input >> mass) {
    particles->push_back(Particle(x, y, z, mass));
  }
  return particles;
}

vector<Vec3<double>> nbody(const vector<Particle> * particles) {
  vector<Vec3<double>> forces;
  for (auto it1 = particles->begin(); it1 < particles->end(); ++it1) {
    Vec3<double> res(0, 0, 0);
    for (auto it2 = particles->begin(); it2 < particles->end(); ++it2) {
      if (it1 == it2) continue;
      Vec3<double> diff = it1->position - it2->position;
      double bottom = pow(diff.sqrSize() + SOFTENING_FACTOR_SQR, 1.5);
      double massTotal = GRAVITATION_CONSTANT * it1->mass * it2->mass;
      res += diff * massTotal / bottom;
    }
    forces.push_back(res);
  }
  return forces;
}

void moveParticles(vector<Particle> * particles, const vector<Vec3<double>>& forces) {
  auto particleIt = particles->begin();
  auto forceIt = forces.begin();
  while(particleIt < particles->end() && forceIt < forces.end()) {
    // particleIt->position += *forceIt / particleIt->mass;
    particleIt->position += *forceIt;
    ++particleIt;
    ++forceIt;
  }
}

void printParticles(const vector<Particle> * particles, ostream& out) {
  for(auto it = particles->begin(); it < particles->end(); ++it) {
    cout << *it << endl;
  }
}
