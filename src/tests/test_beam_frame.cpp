#include "elements/BeamElement.hpp"
#include "fem/Assembler.hpp"
#include "fem/Solver.hpp"
#include "mesh/Mesh.hpp"
#include <fmt/core.h>
#include <Eigen/Dense>

// -----------------------------------------------------------------
// Test: 3-beam structure (planar frame)
//
//   Node 0 (left fixed)  -- elem 0 --  Node 1 (free node)
//                                          |
//                                        elem 2
//                                          |
//   Node 2 (right fixed) -- elem 1 --  Node 1 (free node)
//                                      Node 3 (bottom fixed)
//
//   Loads on Node 1: F downward, M about z
// -----------------------------------------------------------------

void test_beam_frame() {
    fmt::print("\n=== TEST: Beam Frame ===\n");

    // -- Data ------------------------------------------------------
    double E   = 210e9;          // Pa
    double nu  = 0.3;            // Poisson
    double A   = 1e-3;           // m^2  (1000 mm^2)
    double J   = 1e-5;           // m^4
    double Iy  = J;              // equal to J for simplicity
    double Iz  = J;
    double L   = 0.5;            // m   (500 mm)
    double F   = 1e5;            // N
    double M   = 5e4;            // Nm

    // -- Nodes -----------------------------------------------------
    // Node 0: left fixed     (-L, 0, 0)
    // Node 1: central free     (0, 0, 0)
    // Node 2: right fixed    (+L, 0, 0)
    // Node 3: bottom fixed     (0,-L, 0)
    Mesh mesh;
    mesh.nodes = {
        {0, Eigen::Vector3d(-L,  0, 0)},
        {1, Eigen::Vector3d( 0,  0, 0)},
        {2, Eigen::Vector3d( L,  0, 0)},
        {3, Eigen::Vector3d( 0, -L, 0)},
    };

    // -- Section ---------------------------------------------------
    BeamElement::Section sec { A, Iy, Iz, J };

    // -- Elements --------------------------------------------------
    // elem 0: node 0 -> node 1  (left horizontal beam)
    // elem 1: node 1 -> node 2  (right horizontal beam)
    // elem 2: node 1 -> node 3  (vertical beam)
    mesh.elements = {
        {0, ElementType::BEAM2, {0, 1}, 0, 0.0, 0},
        {1, ElementType::BEAM2, {1, 2}, 0, 0.0, 0},
        {2, ElementType::BEAM2, {1, 3}, 0, 0.0, 0},
    };

    // -- Material --------------------------------------------------
    std::vector<Material> materials = {
        {0, "steel", E, nu, 7850}
    };

    // -- Sections --------------------------------------------------
    std::vector<BeamSection> sections = {
        {0, "S1", A, Iy, Iz, J}
    };

    // -- Assemble K ------------------------------------------------
    Assembler assembler;
    auto K = assembler.assemble(mesh, materials, sections);

    int totalDOF = assembler.totalDOF();
    fmt::print("Total DOF: {}\n", totalDOF);

    // -- Force vector f --------------------------------------------
    // Node 1 has offset = 6 (every node has 6 DOF)
    // Node 1 DOFs: [ux=6, uy=7, uz=8, rotx=9, roty=10, rotz=11]
    Eigen::VectorXd f = Eigen::VectorXd::Zero(totalDOF);
    int node1offset = 1 * 6;
    f(node1offset + 1) = -F;   // Fy = -F (downward)
    f(node1offset + 5) =  M;   // Mz = M  (moment about z)

    // -- Boundary conditions ---------------------------------------
    // Fixed supports: node 0, node 2, node 3 -> all 6 DOF locked
    std::vector<int> fixedDOFs;
    for (int node : {0, 2, 3}) {
        int offset = node * 6;
        for (int d = 0; d < 6; ++d)
            fixedDOFs.push_back(offset + d);
    }

    // -- Solve -----------------------------------------------------
    Solver solver;
    auto result = solver.solve(K, f, fixedDOFs);

    if (!result.converged) {
        fmt::print("FAILED: {}\n", result.message);
        return;
    }

    // -- Results at the free node (node 1) -------------------------
    fmt::print("\n=== Results at free node (node 1) ===\n");
    fmt::print("  ux     = {:+.6e} m\n",   result.displacements[node1offset + 0]);
    fmt::print("  uy     = {:+.6e} m\n",   result.displacements[node1offset + 1]);
    fmt::print("  uz     = {:+.6e} m\n",   result.displacements[node1offset + 2]);
    fmt::print("  rotx   = {:+.6e} rad\n", result.displacements[node1offset + 3]);
    fmt::print("  roty   = {:+.6e} rad\n", result.displacements[node1offset + 4]);
    fmt::print("  rotz   = {:+.6e} rad\n", result.displacements[node1offset + 5]);

    // TODO: add an analytical reference value and assert the result is
    // within tolerance, so this becomes a real pass/fail regression test
    // instead of just printing numbers.
}

int main() {
    test_beam_frame();
    return 0;
}