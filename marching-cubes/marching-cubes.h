#pragma once

#include "as/as-math-ops.hpp"

#include <vector>

struct point_t
{
    float val_;
    as::vec3_t position_;
};

struct cell_t
{
    point_t points_[8];
};

struct triangle_t
{
    as::vec3_t verts_[3];
};

point_t*** create_point_volume(int dimension);
cell_t*** create_cell_volume(int dimension);

void generate_point_data(
    point_t*** points, int dimension, const as::vec3_t& offset, float scale);
void generate_cell_data(cell_t*** cells, point_t*** points, int dimension);

void destroy_point_volume(point_t*** points, int dimension);
void destroy_cell_volume(cell_t*** cells, int dimension);

std::vector<triangle_t> march(cell_t*** cells, int dimension, float threshold);
