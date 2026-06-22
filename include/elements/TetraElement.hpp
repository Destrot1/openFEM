#pragma once
#include "elements/Element.hpp"

class TetraElement : public Element {
public:
    Eigen::MatrixXd stiffnessMatrix(
        const std::vector<Eigen::Vector3d>& nodes,
        double E,
        double nu,
        double thickness = 0.0
    ) const override;

    int numNodes()   const override { return 4; }
    int dofPerNode() const override { return 3; }

    // Constant stress over the element (Tet4 has linear shape
    // functions, so strain/stress don't vary within it): sigma = D*B*ue.
    // ue is the element's 12 nodal displacements (ux,uy,uz per node, in
    // the same order as `nodes`). Returns [sxx,syy,szz,sxy,syz,szx], in
    // the same units as E (MPa if E is in MPa).
    Eigen::VectorXd stress(
        const std::vector<Eigen::Vector3d>& nodes,
        double E,
        double nu,
        const Eigen::VectorXd& ue
    ) const;
};