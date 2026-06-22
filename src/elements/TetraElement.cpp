#include "elements/TetraElement.hpp"
#include <stdexcept>

namespace {

// Strain-displacement matrix B (6x12) and the tetrahedron volume --
// shared by stiffnessMatrix() and stress(), both need the same
// shape-function derivatives. B relates nodal displacements to the 6
// strains: [exx, eyy, ezz, gxy, gyz, gxz].
Eigen::MatrixXd tetraBMatrix(const std::vector<Eigen::Vector3d>& p, double& volumeOut) {
    if ((int)p.size() != 4)
        throw std::invalid_argument("Tet4 needs exactly 4 nodes");

    // Build the matrix C of homogeneous coordinates: C * [a,b,c,d]^T = u
    // for each node. Its inverse gives the shape function derivatives.
    Eigen::Matrix4d C;
    for (int i = 0; i < 4; ++i)
        C.row(i) << 1.0, p[i].x(), p[i].y(), p[i].z();

    double detC = C.determinant();
    if (std::abs(detC) < 1e-14)
        throw std::runtime_error("Degenerate tetrahedron (zero volume)");

    volumeOut = std::abs(detC) / 6.0;

    Eigen::Matrix4d Ci = C.inverse();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 12);
    for (int i = 0; i < 4; ++i) {
        double bi = Ci(1, i);  // dNi/dx
        double ci = Ci(2, i);  // dNi/dy
        double di = Ci(3, i);  // dNi/dz
        int    col = 3 * i;

        B(0, col)   = bi;
        B(1, col+1) = ci;
        B(2, col+2) = di;
        B(3, col)   = ci;  B(3, col+1) = bi;
        B(4, col+1) = di;  B(4, col+2) = ci;
        B(5, col)   = di;  B(5, col+2) = bi;
    }
    return B;
}

// Constitutive matrix D (6x6) -- linear isotropic 3D.
Eigen::MatrixXd tetraDMatrix(double E, double nu) {
    double k = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double a = k * (1.0 - nu);
    double b = k * nu;
    double c = k * (1.0 - 2.0 * nu) / 2.0;

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
    D(0,0)=a; D(0,1)=b; D(0,2)=b;
    D(1,0)=b; D(1,1)=a; D(1,2)=b;
    D(2,0)=b; D(2,1)=b; D(2,2)=a;
    D(3,3)=c; D(4,4)=c; D(5,5)=c;
    return D;
}

} // namespace

Eigen::MatrixXd TetraElement::stiffnessMatrix(
    const std::vector<Eigen::Vector3d>& p,
    double E, double nu, double) const
{
    double volume;
    Eigen::MatrixXd B = tetraBMatrix(p, volume);
    Eigen::MatrixXd D = tetraDMatrix(E, nu);

    // Ke = B^T * D * B * V -- for Tet4, B is constant so the integral
    // over the element volume is trivial.
    return B.transpose() * D * B * volume;
}

Eigen::VectorXd TetraElement::stress(
    const std::vector<Eigen::Vector3d>& p,
    double E, double nu,
    const Eigen::VectorXd& ue) const
{
    double volume;
    Eigen::MatrixXd B = tetraBMatrix(p, volume);
    Eigen::MatrixXd D = tetraDMatrix(E, nu);
    return D * B * ue;
}
