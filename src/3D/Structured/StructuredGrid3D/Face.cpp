#include "Face.h"
#include "StructuredGrid3D.h"

const std::array<Face::Direction, 6> Face::DIRECTIONS = {
    Face::I_POS,
    Face::I_NEG,
    Face::J_POS,
    Face::J_NEG,
    Face::K_POS,
    Face::K_NEG
};

Face::Face(const StructuredGrid3D &grid, Face::Direction f, Label i, Label j, Label k)
    :
      _grid(grid),
      _i(f == I_POS ? i + 1 : i),
      _j(f == J_POS ? j + 1 : j),
      _k(f == K_POS ? k + 1 : k)
{
    std::array<Point3D, 4> nodes;

    switch(f)
    {
    case I_POS: case I_NEG:
        nodes = {
            _grid.node(_i, j, k),
            _grid.node(_i, j + 1, k),
            _grid.node(_i, j + 1, k + 1),
            _grid.node(_i, j, k + 1)
        };
        break;

    case J_POS: case J_NEG:
        nodes = {
            _grid.node(i, _j, k),
            _grid.node(i, _j, k + 1),
            _grid.node(i + 1, _j, k + 1),
            _grid.node(i + 1, _j, k)
        };
        break;

    case K_POS: case K_NEG:
        nodes = {
            _grid.node(i, j, _k),
            _grid.node(i + 1, j, _k),
            _grid.node(i + 1, j + 1, _k),
            _grid.node(i, j + 1, _k)
        };
        break;
    }

    _shape = Rectangle(std::accumulate(nodes.begin(), nodes.end(), Point3D(0., 0., 0.)) / 4., cross(nodes[1] - nodes[0], nodes[3] - nodes[0]));
}
