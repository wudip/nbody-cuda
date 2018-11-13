#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

#include "cell.cpp"

// Softening factor squared
#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11

using namespace std;

vector<Particle> * loadParticles(istream& input);
vector<Vec3<double>> nbody(const vector<Particle> * particles);
void moveParticles(vector<Particle> * particles, const vector<Vec3<double>>& forces);
void printParticles(const vector<Particle> * particles, ostream& out);
vector<Vec3<double>> nbodyBarnesHut(Cell & cell);

int main(int argc, char** argv) {
  vector<Particle> * particles = loadParticles(cin);
  clock_t clk_start = clock();
  for(int i = 0; i < 1; ++i) {
    double min[3] = {-5., -5., -5.};
    double max[3] = {20., 20., 20.};
    Cell octree(min, max);
    for(auto it = particles->begin(); it < particles->end(); ++it) {
        octree.add(&*it);
    }
    vector<Vec3<double>> forcesEstimation = nbodyBarnesHut(octree);
    vector<Vec3<double>> forces = nbody(particles);
    moveParticles(particles, forces);
  }
  clock_t clk_end = clock();
  cout << "Time: " << (clk_end - clk_start) << " ms" << endl;
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
      Vec3<double> diff = it1->getPosition() - it2->getPosition();
      double bottom = pow(diff.sqrSize() + SOFTENING_FACTOR_SQR, 1.5);
      double massTotal = GRAVITATION_CONSTANT * it1->mass * it2->mass;
      res += diff * massTotal / bottom;
    }
    forces.push_back(res);
  }
  return forces;
}

vector<Vec3<double>> nbodyBarnesHut(Cell & cell) {
    return vector<Vec3<double>>();
}

void moveParticles(vector<Particle> * particles, const vector<Vec3<double>>& forces) {
  auto particleIt = particles->begin();
  auto forceIt = forces.begin();
  while(particleIt < particles->end() && forceIt < forces.end()) {
    Vec3<double> acceleration = *forceIt / particleIt->mass;
    particleIt->accelerate(acceleration);
    particleIt->updatePosition();
    ++particleIt;
    ++forceIt;
  }
}

void printParticles(const vector<Particle> * particles, ostream& out) {
  for(auto it = particles->begin(); it < particles->end(); ++it) {
    cout << *it << endl;
  }
}
