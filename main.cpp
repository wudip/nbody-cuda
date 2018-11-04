#include <iostream>
#include <vector>
#include <math.h>

// Softening factor squared
#define SOFTENING_FACTOR_SQR 0.5

using namespace std;

struct Particle {
  double x;
  double y;
  double mass;
  Particle(double _x, double _y, double _mass): x(_x), y(_y), mass(_mass) {};
};

struct MoveVector {
  double x;
  double y;
  MoveVector(): x(0), y(0) {};
  MoveVector(double _x, double _y): x(_x), y(_y) {};
};

vector<Particle> loadParticles(istream& input);
vector<MoveVector> nbody(const vector<Particle>& particles);
void moveParticles(vector<Particle>& particles, const vector<MoveVector>& forces);
void printParticles(const vector<Particle>& particles, ostream& out);

int main(int argc, char** argv) {
  vector<Particle> particles = loadParticles(cin);
  vector<MoveVector> forces = nbody(particles);
  moveParticles(particles, forces);
  printParticles(particles, cout);
  return 0;
}

vector<Particle> loadParticles(istream& input) {
  vector<Particle> particles;
  double x;
  double y;
  double mass;
  while(input >> x && input >> y && input >> mass) {
    particles.push_back(Particle(x, y, mass));
  }
  return particles;
}


vector<MoveVector> nbody(const vector<Particle>& particles) {
  vector<MoveVector> forces;
  for (auto it1 = particles.begin(); it1 < particles.end(); ++it1) {
    double forceX = 0;
    double forceY = 0;
    for (auto it2 = particles.begin(); it2 < particles.end(); ++it2) {
      if (it1 == it2) continue;
      double distX = it1->x - it2->x;
      double distY = it1->y - it2->y;
      double vectorSizeSqr = distX * distX + distY * distY;
      double bottom = pow(vectorSizeSqr + SOFTENING_FACTOR_SQR, 1.5);
      double massTotal = it1->mass * it2->mass;
      forceX -= (massTotal * distX) / bottom;
      forceY -= (massTotal * distY) / bottom;
    }
    forces.push_back(MoveVector(forceX, forceY));
  }
  return forces;
}

void moveParticles(vector<Particle>& particles, const vector<MoveVector>& forces) {
  auto particleIt = particles.begin();
  auto forceIt = forces.begin();
  while(particleIt < particles.end() && forceIt < forces.end()) {
    particleIt->x += forceIt->x;
    particleIt->y += forceIt->y;
    ++particleIt;
    ++forceIt;
  }
}

void printParticles(const vector<Particle>& particles, ostream& out) {
  for(auto it = particles.begin(); it < particles.end(); ++it) {
    double x = it->x;
    double y = it->y;
    double mass = it->mass;
    out << x << " " << y << " " << mass << endl;
  }
}
