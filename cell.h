#ifndef BAKAJ_WUDI_CUDA_CELL_H
#define BAKAJ_WUDI_CUDA_CELL_H

#include <iostream>

#include "vec3.h"
#include "particle.h"

#define NUM_OF_SUBCELLS 8
#define NUM_OF_DIMENSIONS 3

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
    void printCell(std::ostream & ost, int & id) const;
    void getForceSiblings(const Cell * c, Vec3<double>& forces) const;

public:
    Cell() = delete;
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS]);
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS], Cell * daddy);
    ~Cell();
    void add(Particle *);
    Cell * getSubcell(Particle * particle);
    void printGraph(std::ostream & ost) const;
    void updateCenter();
    Vec3<double> getForce() const;
};

void addToForces(Vec3<double>& forces, Particle* particle, Particle& sibPart);

#endif //BAKAJ_WUDI_CUDA_CELL_H
