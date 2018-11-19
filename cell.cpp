#include <iostream>
#include <cmath>

#include "particle.cpp"

#define NUM_OF_SUBCELLS 8
#define NUM_OF_DIMENSIONS 3

#define SOFTENING_FACTOR_SQR 0.5
#define GRAVITATION_CONSTANT 6.67300E-11

using namespace std;

void addToForces(Vec3<double>& forces, Particle* particle, Particle& sibPart) {
    Vec3<double> diff = particle->getPosition() - sibPart.getPosition();
    double bottom = pow(diff.sqrSize() + SOFTENING_FACTOR_SQR, 1.5);
    double massTotal = GRAVITATION_CONSTANT * particle->mass * sibPart.mass;
    forces += diff * massTotal / bottom;
}

/**
 * Octree cell
 */
class Cell {
    Particle * part;
    Particle center;
    Cell * parent;
    Cell * subtree[NUM_OF_SUBCELLS];
    Vec3<double> minPoint;
    Vec3<double> maxPoint;
    /**
     * Splits the cell into 2^{@link #NUM_OF_DIMENSIONS} cells. They are stored in {@link #subtree} variable.
     */
    void split();
    void printCell(ostream & ost, int & id) const;
    void getForceSiblings(const Cell * c, Vec3<double>& forces) const {
        if (parent) {
            for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
                Cell * sibling = parent->subtree[i];
                if (sibling == this) continue;
                Particle& sibPart = sibling->center;
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
    void updateCenter();
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

void Cell::updateCenter() {
    Vec3<double> cpos();
    double mass = 0;
    if (part) {
        Vec3<double> position = part->getPosition();
        double m = part->mass;
        cpos += position * m;
        mass += m;
    }
    if (subtree[0]) {
        for (int i = 0; i < NUM_OF_DIMENSIONS; ++i) {
            subtree[i]->updateCenter();
            Particle& subCenter = subtree[i]->center;
            cpos += subCenter.getPosition() * subCenter.mass;
            mass += subCenter.mass;
        }
    }
    if (mass > 0) {
        cpos /= mass;
    }
    center = Particle(cpos.x, cpos.y, cpos.z, mass);
}
