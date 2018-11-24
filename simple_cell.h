#ifndef BAKAJ_WUDI_CUDA_SIMPLE_CELL_H
#define BAKAJ_WUDI_CUDA_SIMPLE_CELL_H

/**
 * Simplified cell used in particle position computation
 * It contains no pointer whatsoever so it could be used easily on GPU
 * Position of parent/children is computed by position in the array (so whole tree must be saved in an array according
 * to the indexes.
 */

// TODO: code duplication
#define NUM_OF_SUBCELLS 8
#define NUM_OF_DIMENSIONS 3

class SimpleCell {
    /**
     * Position of particle in particle array (or {@code (unsigned int) -1} if it contains no particle)
     */
    unsigned int particle;

    /**
     * Center of mass of the cell
     */
    Particle center; // TODO: maybe unnecessary velocity and cellIndex

    /**
     * Position of parent in the cell array (or {@code (unsigned int) -1} if no parent is present)
     */
    unsigned int parent;

    /**
     * Positions of all children in the cell array (first item contains {@code (unsigned int) -1} if it has no children)
     */
    unsigned int subtree[NUM_OF_SUBCELLS];

};

#endif //BAKAJ_WUDI_CUDA_SIMPLE_CELL_H
