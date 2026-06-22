#pragma once
#include <Eigen/Dense>
#include <vector>

class Element {
public:
    virtual ~Element() = default;

    virtual Eigen::MatrixXd stiffnessMatrix(
        const std::vector<Eigen::Vector3d>& nodes,
        double E,
        double nu,
        double thickness = 0.0
    ) const = 0;

    virtual int numNodes()   const = 0;
    virtual int dofPerNode() const = 0;
    int totalDOF() const { return numNodes() * dofPerNode(); }
};