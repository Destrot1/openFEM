#pragma once
#include "fem/Model.hpp"
#include "fem/Assembler.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Mesher.hpp"
#include <Eigen/Dense>
#include <vector>

// Converts entity-scoped BoundaryCondition/Load (dim/tag only) into the
// DOF-level data (global indices, force vector) Assembler/Solver understand.
class LoadResolver {
public:
    // Returns the global DOF indices to lock, one per locked DOF across all bcs.
    // DOFs that don't exist on a node (e.g. rotation on a 3-DOF node) are skipped.
    static std::vector<int> resolveFixedDOFs(
        const std::vector<BoundaryCondition>& bcs,
        const Assembler&                       assembler,
        const Mesher&                          mesher
    );

    // Builds the global force vector, splitting each load across its entity's
    // nodes by tributary area/length, not evenly (mesh density would skew that).
    static Eigen::VectorXd resolveForces(
        const Mesh&                mesh,
        const std::vector<Load>&  loads,
        const Assembler&           assembler,
        const Mesher&              mesher
    );
};
