#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>

#include <iostream>
#include <cmath>

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

    T sqrSize() {
        return x * x + y * y + z * z;
    }

    template <typename U>
    friend std::ostream& operator<<( std::ostream &o, const Vec3<U> &vector) {
        o << vector.x << " " << vector.y << " " << vector.z;
        return o;
    }

    double getDim(int dimension) const {
        if (dimension == 0) return x;
        if (dimension == 1) return y;
        if (dimension == 2) return z;
        throw "Dimension out of range";
    }

    void setDim(int dimension, const double value) {
        if (dimension == 0) x = value;
        else if (dimension == 1) y = value;
        else if (dimension == 2) z = value;
        else throw "Dimension out of range";
    }
};

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


#define NUM_OF_SUBCELLS 8
#define NUM_OF_DIMENSIONS 3

#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11

using namespace std;

void addToForces(Vec3<double>& forces, Particle* particle, Particle* sibPart) {
    Vec3<double> diff = particle->getPosition() - sibPart->getPosition();
    double bottom = pow(diff.sqrSize() + SOFTENING_FACTOR_SQR, 1.5);
    double massTotal = GRAVITATION_CONSTANT * particle->mass * sibPart->mass;
    forces += diff * massTotal / bottom;
}

/**
 * Octree cell
 */
class Cell {
    Particle * part;
    Cell * parent;
    Cell * subtree[NUM_OF_SUBCELLS];
    Vec3<double> minPoint;
    Vec3<double> maxPoint;
    /**
     * Splits the cell into 2^{@link #NUM_OF_DIMENSIONS} cells. They are stored in {@link #subtree} variable.
     */
    void split();
    void printCell(ostream & ost, int & id) const;
    void getCenter(Vec3<double> & center, double & mass) const;
    void getForceSiblings(const Cell * c, Vec3<double>& forces) const {
        if (parent) {
            for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
                Cell * sibling = parent->subtree[i];
                if (sibling == this) continue;
                Particle* sibPart = sibling->getCenter();
                Particle* particle = c->part;
                addToForces(forces, particle, sibPart);
            }
        }
    }

public:
    Cell() = delete;
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS]);
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS], Cell * daddy);
    ~Cell();
    void add(Particle *);
    Cell * getSubcell(Particle * particle);
    void printGraph(ostream & ost) const;
    Particle * getCenter() const;
    Vec3<double> getForce() const {
        Vec3<double> potato(0, 0, 0);
        const Cell * c = this;
        while (c->parent) {
            getForceSiblings(c, potato);
            c = c->parent;
        }
        return potato;
    }
};

Cell::Cell(const double *boundMin, const double *boundMax, Cell * daddy):
        part(nullptr), parent(daddy), subtree{nullptr, nullptr, nullptr}
{
    for (int dim = 0; dim < NUM_OF_DIMENSIONS; ++dim) {
        minPoint.setDim(dim, boundMin[dim]);
        maxPoint.setDim(dim, boundMax[dim]);
    }
}

Cell::Cell(const double *boundMin, const double *boundMax):
        Cell(boundMin, boundMax, nullptr)
{
}

Cell::~Cell() {
    if (subtree[0]) {
        for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
            delete subtree[i];
        }
    }
}

void Cell::add(Particle * particle) {
    if (subtree[0]) {
        Cell * cell = getSubcell(particle);
        cell->add(particle);
    }
    else if (part) {
        split();
        add(part);
        add(particle);
    }
    else {
        part = particle;
        particle->cell = this;
    }
}

void Cell::split() {
    double tempBoundMin[NUM_OF_DIMENSIONS];
    double tempBoundMax[NUM_OF_DIMENSIONS];
    for (int subcell = 0; subcell < NUM_OF_SUBCELLS; ++subcell) {
        for (int dimension = 0; dimension < NUM_OF_DIMENSIONS; ++dimension) {
            int dimThBit = (1 << dimension) & subcell; // if true, will get center + max
            double min = minPoint.getDim(dimension);
            double max = maxPoint.getDim(dimension);
            double center = (min + max) / 2;
            if (dimThBit) {
                tempBoundMin[dimension] = center;
                tempBoundMax[dimension] = max;
            }
            else {
                tempBoundMin[dimension] = min;
                tempBoundMax[dimension] = center;
            }
        }
        subtree[subcell] = new Cell(tempBoundMin, tempBoundMax, this);
    }
}

Cell *Cell::getSubcell(Particle * particle) {
    int index = 0;
    for (int dim = 0; dim < NUM_OF_DIMENSIONS; ++dim) {
        double min = minPoint.getDim(dim);
        double max = maxPoint.getDim(dim);
        double center = (min + max) / 2;
        double posDim = particle->getPosition().getDim(dim);
        if (posDim > center) {
            index += (int) pow(2, dim);
        }
    }
    return subtree[index];
}

// DEBUG

void Cell::printCell(ostream & ost, int & id) const {
    int myId = id++;
    ost << "a" << myId << "[label=\"" << minPoint << "\\n" << maxPoint << "\\n";
    if (part) ost << part->getPosition();
    else ost << "---";
    ost << "\"]" << endl;
    if (!subtree[0]) return;
    for (int subcell = 0; subcell < NUM_OF_SUBCELLS; ++subcell) {
        ost << "a" << myId << " -> a" << id << endl;
        subtree[subcell]->printCell(ost, id);
    }
}

void Cell::printGraph(ostream & ost) const {
    ost << "digraph koprovkaaknedliky {" << endl;
    int id = 0;
    printCell(ost, id);
    ost << "}" << endl;
}

Particle * Cell::getCenter() const {
    Vec3<double> center(0, 0, 0);
    double mass = 0;
    if (part || subtree[0]) {
        getCenter(center, mass);
        center /= mass;
    }
    Particle * particle = new Particle(center.x, center.y, center.z, mass);
    return particle;
}

void Cell::getCenter(Vec3<double> & center, double & mass) const {
    if (part) {
        Vec3<double> position = part->getPosition();
        double m = part->mass;
        center += position * m;
        mass += m;
    }
    if (subtree[0]) {
        for (int i = 0; i < NUM_OF_DIMENSIONS; ++i) {
            subtree[i]->getCenter(center, mass);
        }
    }
}

// Softening factor squared
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
        forces.push_back(f);
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
