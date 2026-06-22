#include "elements/BeamElement.hpp"
#include <stdexcept>
#include <cmath>

Eigen::MatrixXd BeamElement::localStiffness(
    double L, double E, double G) const
{
    double A  = m_section.A;
    double Iy = m_section.Iy;
    double Iz = m_section.Iz;
    double J  = m_section.J;

    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(12, 12);

    // -- Axial tension/compression (DOF 0 and 6) ------------------
    double ea = E * A / L;
    Ke(0,0)= ea; Ke(0,6)=-ea;
    Ke(6,0)=-ea; Ke(6,6)= ea;

    // -- Torsion (DOF 3 and 9) ------------------------------------
    double gj = G * J / L;
    Ke(3,3)= gj; Ke(3,9)=-gj;
    Ke(9,3)=-gj; Ke(9,9)= gj;

    // -- Bending in the XZ plane -> about Y axis -> uses Iy --------
    // DOF involved: w (2,8) and rot_y (4,10)
    double eiy = E * Iy / (L*L*L);
    Ke(2,2) = 12*eiy;    Ke(2,4)  = 6*L*eiy;
    Ke(2,8) =-12*eiy;    Ke(2,10) = 6*L*eiy;
    Ke(4,2) = 6*L*eiy;   Ke(4,4)  = 4*L*L*eiy;
    Ke(4,8) =-6*L*eiy;   Ke(4,10) = 2*L*L*eiy;
    Ke(8,2) =-12*eiy;    Ke(8,4)  =-6*L*eiy;
    Ke(8,8) = 12*eiy;    Ke(8,10) =-6*L*eiy;
    Ke(10,2)= 6*L*eiy;   Ke(10,4) = 2*L*L*eiy;
    Ke(10,8)=-6*L*eiy;   Ke(10,10)= 4*L*L*eiy;

    // -- Bending in the XY plane -> about Z axis -> uses Iz --------
    // DOF involved: v (1,7) and rot_z (5,11)
    double eiz = E * Iz / (L*L*L);
    Ke(1,1) = 12*eiz;    Ke(1,5)  =-6*L*eiz;
    Ke(1,7) =-12*eiz;    Ke(1,11) =-6*L*eiz;
    Ke(5,1) =-6*L*eiz;   Ke(5,5)  = 4*L*L*eiz;
    Ke(5,7) = 6*L*eiz;   Ke(5,11) = 2*L*L*eiz;
    Ke(7,1) =-12*eiz;    Ke(7,5)  = 6*L*eiz;
    Ke(7,7) = 12*eiz;    Ke(7,11) = 6*L*eiz;
    Ke(11,1)=-6*L*eiz;   Ke(11,5) = 2*L*L*eiz;
    Ke(11,7)= 6*L*eiz;   Ke(11,11)= 4*L*L*eiz;

    return Ke;
}

Eigen::MatrixXd BeamElement::rotationMatrix(
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& p1) const
{
    // Local x axis = beam direction
    Eigen::Vector3d ex = (p1 - p0).normalized();

    // Local y axis: pick a reference vector not parallel to ex
    Eigen::Vector3d ref = (std::abs(ex.z()) < 0.9)
                          ? Eigen::Vector3d(0,0,1)
                          : Eigen::Vector3d(0,1,0);

    Eigen::Vector3d ez = ex.cross(ref).normalized();
    Eigen::Vector3d ey = ez.cross(ex).normalized();

    // 3x3 rotation matrix local -> global
    Eigen::Matrix3d R3;
    R3.row(0) = ex;
    R3.row(1) = ey;
    R3.row(2) = ez;

    // Expand to 12x12 (four 3x3 blocks on the diagonal)
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(12, 12);
    for (int i = 0; i < 4; ++i)
        R.block<3,3>(3*i, 3*i) = R3;

    return R;
}

Eigen::MatrixXd BeamElement::stiffnessMatrix(
    const std::vector<Eigen::Vector3d>& p,
    double E, double nu, double) const
{
    if ((int)p.size() != 2)
        throw std::invalid_argument("Beam2 needs exactly 2 nodes");

    double L = (p[1] - p[0]).norm();
    if (L < 1e-14)
        throw std::runtime_error("Degenerate beam (zero length)");

    // Shear modulus
    double G = E / (2.0 * (1.0 + nu));

    // Ke in the local frame
    Eigen::MatrixXd Ke_local = localStiffness(L, E, G);

    // Rotate to the global frame: Ke_global = R^T * Ke_local * R
    Eigen::MatrixXd R = rotationMatrix(p[0], p[1]);
    return R.transpose() * Ke_local * R;
}