#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>

#include "cell.cpp"

// Softening factor squared
#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11
#define WINDOWS_SUCKS 1

using namespace std;

vector<Particle> * loadParticles(istream& input);
vector<Vec3<double>> nbody(const vector<Particle> * particles);
void moveParticles(vector<Particle> * particles, const vector<Vec3<double>>& forces);
void printParticles(const vector<Particle> * particles, ostream& out);
vector<Vec3<double>> nbodyBarnesHut(const vector<Particle> * particles, Cell & cell);
double * computeParticleBoundaries(const vector<Particle> * particles);

int main(int argc, char** argv) {
  vector<Particle> * particles;
  if (WINDOWS_SUCKS && argc > 1) {
      std::ifstream ifs;
      ifs.open(argv[1], ifstream::in);
      particles = loadParticles(ifs);
      ifs.close();
  }
  else {
      particles = loadParticles(cin);
  }
  clock_t clk_start = clock();
  for(int i = 0; i < 1000; ++i) {
    double * particleBoundaries = computeParticleBoundaries(particles);
    Cell octree(particleBoundaries, particleBoundaries + 3);
    delete[] particleBoundaries;
    for(auto it = particles->begin(); it < particles->end(); ++it) {
        octree.add(&*it);
    }
    octree.updateCenter();
    vector<Vec3<double>> forces = nbodyBarnesHut(particles, octree);
    //vector<Vec3<double>> forces = nbody(particles);
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
  forces.reserve(particles->size());
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

vector<Vec3<double>> nbodyBarnesHut(const vector<Particle> * particles, Cell & cell) {
    vector<Vec3<double>> forces;
    forces.reserve(particles->size());
    #pragma acc parallel loop
    for (auto partit = particles->begin(); partit < particles->end(); ++partit) {
        Vec3<double> f = partit->cell->getForce();
        int index = (int)(partit - particles->begin());
        forces.at(index) = f;
    }
    return forces;
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

/**
 * Computes coordinates of smallest possible box containing all the particles
 * @return array of minimum point and maximum point: [min_x, min_y, min_z, max_x, max_y, max_z]
 */
double * computeParticleBoundaries(const vector<Particle> * particles) {
    double * result = new double[6];
    for (int i = 0; i < 3; ++i) {
        result[i] = INFINITY;
        result[3 + i] = -INFINITY;
    }
    for (auto it = particles->begin(); it < particles->end(); ++it) {
        Vec3<double> pos = it->getPosition();
        for (int dim = 0; dim < 3; ++dim) {
            int coord = pos.getDim(dim);
            if (coord < result[dim]) {
                result[dim] = coord;
            }
            if (coord > result[3 + dim]) {
                result[3 + dim] = coord;
            }
        }
    }
    return result;
}

void printParticles(const vector<Particle> * particles, ostream& out) {
  for(auto it = particles->begin(); it < particles->end(); ++it) {
    cout << *it << endl;
  }
}
