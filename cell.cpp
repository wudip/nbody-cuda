#include <iostream>
#include <cmath>

#include "cell.h"

using namespace std;

Cell::Cell() {}

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

Cell::Cell(const Cell & buddy):
    part(buddy.part), center(buddy.center), parent(buddy.parent), minPoint(buddy.minPoint), maxPoint(buddy.maxPoint)
{
    for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
        subtree[i] = buddy.subtree[i];
    }
}

Cell::~Cell() {
    if (subtree[0]) {
        for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
            delete subtree[i];
        }
    }
}

Cell& Cell::operator=(const Cell& buddy)       {
    part = buddy.part;
    center = buddy.center;
    parent = buddy.parent;
    for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
        subtree[i] = buddy.subtree[i];
    }
    minPoint = buddy.minPoint;
    maxPoint = buddy.maxPoint;
    return *this;

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
    Vec3<double> cpos;
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

unsigned int Cell::getNumOfNodes() const {
    unsigned int numOfCells = 1;
    if (subtree[0]) {
        for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
            numOfCells += subtree[i]->getNumOfNodes();
        }
    }
    return numOfCells;
}

Cell* Cell::serialize() const {
    Cell * cells = new Cell[getNumOfNodes()];
    unsigned int index = 0;
    serialize(cells, index);
    //cout << getNumOfNodes() << " nodes, " << index << " filled" << endl;
    return cells;
}

void SimpleCell::serialize(SimpleCell* cellList, unsigned int &index) const {
    cellList[index] = SimpleCell();
    if (part) part->cellIndex = index;
    index++;
    if (subtree[0]) {
        for (int i = 0; i < NUM_OF_SUBCELLS; ++i) {
            subtree[i]->serialize(cellList, index);
        }
    }
}
