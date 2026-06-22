#pragma once
#include <vector>
#include <Eigen/Dense>

struct Node {
    int id;
    Eigen::Vector3d position;
};

enum class ElementType {
    TETRA4,
    TRI3,
    BEAM2
};

struct MeshElement {
    int id;
    ElementType type;
    std::vector<int> nodeIds;
    int    materialId = -1;
    double thickness  = 0.0;
    int    sectionId  = -1;
};

struct Mesh {
    std::vector<Node>        nodes;
    std::vector<MeshElement> elements;

    int totalDOF() const { return (int)nodes.size() * 3; }
};