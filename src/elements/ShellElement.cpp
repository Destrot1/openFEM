#include "elements/ShellElement.hpp"
#include <stdexcept>
#include <cmath>

// -- Utility: triangle area ----------------------------------------
static double triArea(const Eigen::Vector3d& p0,
                      const Eigen::Vector3d& p1,
                      const Eigen::Vector3d& p2) {
    return 0.5 * ((p1 - p0).cross(p2 - p0)).norm();
}

// -- Utility: triangle local reference frame -----------------------
// Returns [e1, e2, e3] where e3 is the normal
static void localFrame(const Eigen::Vector3d& p0,
                       const Eigen::Vector3d& p1,
                       const Eigen::Vector3d& p2,
                       Eigen::Vector3d& e1,
                       Eigen::Vector3d& e2,
                       Eigen::Vector3d& e3) {
    e1 = (p1 - p0).normalized();
    e3 = ((p1 - p0).cross(p2 - p0)).normalized();
    e2 = e3.cross(e1);
}

Eigen::MatrixXd ShellElement::membraneStiffness(
    const std::vector<Eigen::Vector3d>& p,
    double E, double nu, double t) const
{
    // Project the nodes onto the local plane
    Eigen::Vector3d e1, e2, e3;
    localFrame(p[0], p[1], p[2], e1, e2, e3);

    // Local 2D coordinates
    auto toLocal = [&](const Eigen::Vector3d& pt) -> Eigen::Vector2d {
        Eigen::Vector3d d = pt - p[0];
        return { d.dot(e1), d.dot(e2) };
    };

    Eigen::Vector2d r0 = toLocal(p[0]);
    Eigen::Vector2d r1 = toLocal(p[1]);
    Eigen::Vector2d r2 = toLocal(p[2]);

    double x1=r0.x(), y1=r0.y();
    double x2=r1.x(), y2=r1.y();
    double x3=r2.x(), y3=r2.y();

    double A = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    if (A < 1e-14)
        throw std::runtime_error("Degenerate shell triangle");

    // CST (Constant Strain Triangle) shape function derivatives
    double b1 = y2 - y3,  b2 = y3 - y1,  b3 = y1 - y2;
    double c1 = x3 - x2,  c2 = x1 - x3,  c3 = x2 - x1;

    // Membrane B matrix (3x6)
    Eigen::MatrixXd Bm = Eigen::MatrixXd::Zero(3, 6);
    Bm(0,0)=b1; Bm(0,2)=b2; Bm(0,4)=b3;
    Bm(1,1)=c1; Bm(1,3)=c2; Bm(1,5)=c3;
    Bm(2,0)=c1; Bm(2,1)=b1;
    Bm(2,2)=c2; Bm(2,3)=b2;
    Bm(2,4)=c3; Bm(2,5)=b3;
    Bm /= (2.0 * A);

    // Plane-stress D matrix (3x3)
    double k = E / (1.0 - nu * nu);
    Eigen::MatrixXd Dm = Eigen::MatrixXd::Zero(3, 3);
    Dm(0,0)=k;     Dm(0,1)=k*nu;
    Dm(1,0)=k*nu;  Dm(1,1)=k;
    Dm(2,2)=k*(1.0-nu)/2.0;

    return Bm.transpose() * Dm * Bm * A * t;
}

Eigen::MatrixXd ShellElement::bendingStiffness(
    const std::vector<Eigen::Vector3d>& p,
    double E, double nu, double t) const
{
    Eigen::Vector3d e1, e2, e3;
    localFrame(p[0], p[1], p[2], e1, e2, e3);

    auto toLocal = [&](const Eigen::Vector3d& pt) -> Eigen::Vector2d {
        Eigen::Vector3d d = pt - p[0];
        return { d.dot(e1), d.dot(e2) };
    };

    Eigen::Vector2d r0 = toLocal(p[0]);
    Eigen::Vector2d r1 = toLocal(p[1]);
    Eigen::Vector2d r2 = toLocal(p[2]);

    double x1=r0.x(), y1=r0.y();
    double x2=r1.x(), y2=r1.y();
    double x3=r2.x(), y3=r2.y();

    double A = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

    double b1 = y2-y3, b2 = y3-y1, b3 = y1-y2;
    double c1 = x3-x2, c2 = x1-x3, c3 = x2-x1;

    // Bending B matrix (3x9) -- curvatures from rotational DOFs
    // DOF per node: [w, rotx, roty] -> 9 DOF total
    Eigen::MatrixXd Bf = Eigen::MatrixXd::Zero(3, 9);
    Bf(0,2)=b1; Bf(0,5)=b2; Bf(0,8)=b3;
    Bf(1,1)=c1; Bf(1,4)=c2; Bf(1,7)=c3;
    Bf(2,1)=b1; Bf(2,2)=c1;
    Bf(2,4)=b2; Bf(2,5)=c2;
    Bf(2,7)=b3; Bf(2,8)=c3;
    Bf /= (2.0 * A);

    // Bending D matrix
    double d = E * t*t*t / (12.0 * (1.0 - nu*nu));
    Eigen::MatrixXd Df = Eigen::MatrixXd::Zero(3, 3);
    Df(0,0)=d;     Df(0,1)=d*nu;
    Df(1,0)=d*nu;  Df(1,1)=d;
    Df(2,2)=d*(1.0-nu)/2.0;

    return Bf.transpose() * Df * Bf * A;
}

Eigen::MatrixXd ShellElement::stiffnessMatrix(
    const std::vector<Eigen::Vector3d>& p,
    double E, double nu, double t) const
{
    if ((int)p.size() != 3)
        throw std::invalid_argument("Shell Tri3 needs exactly 3 nodes");

    // Final 18x18 Ke -- assembles membrane (6x6) and bending (9x9)
    // into the correct DOFs per node: [ux,uy,uz, rotx,roty,rotz]
    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(18, 18);

    Eigen::MatrixXd Km = membraneStiffness(p, E, nu, t);  // 6x6
    Eigen::MatrixXd Kb = bendingStiffness(p, E, nu, t);   // 9x9

    // Membrane -> in-plane translational DOFs (ux, uy per node)
    // node i: global DOFs 6i, 6i+1
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    Ke(6*i+a, 6*j+b) += Km(2*i+a, 2*j+b);

    // Bending -> [uz, rotx, roty] per node
    // node i: global DOFs 6i+2, 6i+3, 6i+4
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    Ke(6*i+2+a, 6*j+2+b) += Kb(3*i+a, 3*j+b);

    // TODO: the drilling DOF (rotz, local index 5 per node) has no
    // stiffness -> zero diagonal terms. Add a small drilling stiffness
    // to avoid singularities when shells are assembled in the global system.

    return Ke;
}