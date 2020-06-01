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

struct CellValues
{
    float values_[8];
};

struct CellPositions
{
    as::vec3_t points_[8];
};

struct Triangle
{
    Triangle() = default;
    Triangle(const as::vec3_t& a, const as::vec3_t& b, const as::vec3_t& c)
        : verts_{a, b, c}
    {
    }

    as::vec3_t verts_[3];
};

Point*** createPointVolume(int dimension);
CellValues*** createCellValues(int dimension);
CellPositions*** createCellPositions(int dimension);

void generatePointData(
    Point*** points, int dimension, const as::vec3_t& offset, float scale,
    const as::vec3_t& cam);
void generateCellData(
    CellPositions*** cellPositions, CellValues*** cellValues, Point*** points,
    int dimension);

void destroyPointVolume(Point*** points, int dimension);

void destroyCellValues(CellValues*** cells, int dimension);
void destroyCellPositions(CellPositions*** cells, int dimension);

std::vector<Triangle> march(
    CellPositions*** cellPositions, CellValues*** cellValues, int dimension,
    float threshold);

} // namespace mc
