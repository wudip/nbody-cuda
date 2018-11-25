#ifndef BAKAJ_WUDI_CUDA_CELL_H
#define BAKAJ_WUDI_CUDA_CELL_H

#include <iostream>

#include "vec3.h"
#include "particle.h"
#include "simple_cell.h"

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
    void serialize(SimpleCell* cellList, const Particle * partBeginning, unsigned int& index, unsigned int daddy) const;

public:
    Cell();
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS]);
    Cell(const double boundMin[NUM_OF_DIMENSIONS], const double boundMax[NUM_OF_DIMENSIONS], Cell * daddy);
    Cell(const Cell& buddy);
    ~Cell();
    Cell& operator=(const Cell& buddy);
    void add(Particle *);
    Cell * getSubcell(Particle * particle);
    void printGraph(std::ostream & ost) const;
    void updateCenter();
    /**
     * @return number of cells in subtree (including this one)
     */
    unsigned int getNumOfNodes() const;
    /**
     * Creates new array of all cells in subtree
     * @param partBeginning beginning of array of particles (necessary to get position of the particle)
     */
    SimpleCell* serialize(const Particle * partBeginning) const;
};

#endif //BAKAJ_WUDI_CUDA_CELL_H
