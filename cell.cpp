#include <iostream>
#include <math.h>

#include "particle.cpp"

#define NUM_OF_SUBCELLS 8
#define NUM_OF_DIMENSIONS 3

using namespace std;

/**
 * Octree cell
 */
class Cell {
    Particle * part;
    Cell * subtree[NUM_OF_SUBCELLS];
    Vec3<double> minPoint;
    Vec3<double> maxPoint;
    /**
     * Splits the cell into 2^{@link #NUM_OF_DIMENSIONS} cells. They are stored in {@link #subtree} variable.
     */
    void split();
    void splitDimension(int dimension, int index, double boundMin[NUM_OF_DIMENSIONS], double boundMax[NUM_OF_DIMENSIONS]);

public:
    Cell() = delete;
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS]);
    ~Cell();
    void add(Particle *);
    Cell * getSubcell(Particle * particle);

};

Cell::Cell(const double *boundMin, const double *boundMax):
    part(nullptr), subtree{nullptr, nullptr, nullptr}
{
    for (int dim = 0; dim < NUM_OF_DIMENSIONS; ++dim) {
        minPoint.setDim(dim, boundMin[dim]);
        maxPoint.setDim(dim, boundMax[dim]);
    }
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
    }
}

void Cell::split() {
    double * tempBoundMin = new double[NUM_OF_DIMENSIONS];
    double * tempBoundMax = new double[NUM_OF_DIMENSIONS];
    splitDimension(0, 0, tempBoundMin, tempBoundMax);
    delete [] tempBoundMin;
    delete [] tempBoundMax;
}

void Cell::splitDimension(int dimension, int index, double *boundMin, double *boundMax) {
    if (dimension == NUM_OF_DIMENSIONS) {
        subtree[index] = new Cell(boundMin, boundMax);
        return;
    }
    double min = minPoint.getDim(dimension);
    double max = maxPoint.getDim(dimension);
    double center = (min + max) / 2;
    boundMin[dimension] = min;
    boundMax[dimension] = center;
    splitDimension(dimension + 1, index, boundMin, boundMax);
    boundMin[dimension] = center;
    boundMax[dimension] = max;
    splitDimension(dimension + 1, index + (int) pow(2, dimension), boundMin, boundMax);
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
