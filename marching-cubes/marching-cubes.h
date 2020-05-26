#pragma once

#include "as/as-math-ops.hpp"

#include <vector>

namespace mc
{

struct Point
{
    float val_;
    as::vec3_t position_;
};

struct Cell
{
    Point points_[8];
};

struct Triangle
{
    as::vec3_t verts_[3];
};

Point*** createPointVolume(int dimension);
Cell*** createCellVolume(int dimension);

void generatePointData(
    Point*** points, int dimension, const as::vec3_t& offset, float scale);
void generateCellData(Cell*** cells, Point*** points, int dimension);

void destroyPointVolume(Point*** points, int dimension);
void destroyCellVolume(Cell*** cells, int dimension);

std::vector<Triangle> march(Cell*** cells, int dimension, float threshold);

} // namespace mc
