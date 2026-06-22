#pragma once
#include "elements/Element.hpp"

class ShellElement : public Element {
public:
    Eigen::MatrixXd stiffnessMatrix(
        const std::vector<Eigen::Vector3d>& nodes,
        double E,
        double nu,
        double thickness = 0.005
    ) const override;

    int numNodes()   const override { return 3; }
    int dofPerNode() const override { return 6; }

private:
    // Membrane contribution (in-plane)
    Eigen::MatrixXd membraneStiffness(
        const std::vector<Eigen::Vector3d>& p,
        double E, double nu, double t) const;

    // Bending contribution (out-of-plane)
    Eigen::MatrixXd bendingStiffness(
        const std::vector<Eigen::Vector3d>& p,
        double E, double nu, double t) const;
};