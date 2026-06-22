#pragma once
#include "elements/Element.hpp"

class BeamElement : public Element {
public:
    // Beam cross-section
    struct Section {
        double A;    // area              [m^2]
        double Iy;   // inertia about y   [m^4]
        double Iz;   // inertia about z   [m^4]
        double J;    // torsion constant  [m^4]
    };

    explicit BeamElement(const Section& s) : m_section(s) {}

    Eigen::MatrixXd stiffnessMatrix(
        const std::vector<Eigen::Vector3d>& nodes,
        double E,
        double nu,
        double thickness = 0.0
    ) const override;

    int numNodes()   const override { return 2; }
    int dofPerNode() const override { return 6; }

private:
    Section m_section;

    // Local Ke in the beam reference frame
    Eigen::MatrixXd localStiffness(double L, double E, double G) const;

    // 12x12 rotation matrix from local to global
    Eigen::MatrixXd rotationMatrix(const Eigen::Vector3d& p0,
                                   const Eigen::Vector3d& p1) const;
};