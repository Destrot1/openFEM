#pragma once
#include "fem/Model.hpp"
#include "fem/Assembler.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Mesher.hpp"
#include <Eigen/Dense>
#include <vector>

// Translates entity-scoped BoundaryCondition/Load (which only know
// "this applies to entity (dim, tag)") into the DOF-level language the
// Assembler/Solver actually speak (global row indices, a force
// vector). Assembler/Solver never need to know what a "face" or a
// "constraint" is -- this class is the only place that bridges the two
// worlds, using Mesher's live entity queries (entityNodeIndices/
// surfaceTriangles/surfaceLines -- these work for ANY entity, not just
// the ones configured in 'assign') plus the Assembler's node -> DOF
// layout. mesher must still be alive (not shutdown()) for both methods.
class LoadResolver {
public:
    // One global DOF index per locked degree of freedom, across all
    // bcs. A DOF that doesn't exist for a given node (e.g. asking to
    // lock a rotation on a node that only has 3 DOF) is silently
    // skipped -- there's nothing there to lock.
    static std::vector<int> resolveFixedDOFs(
        const std::vector<BoundaryCondition>& bcs,
        const Assembler&                       assembler,
        const Mesher&                          mesher
    );

    // A global force vector (size = assembler.totalDOF()), built by
    // distributing each load's total force across the nodes of its
    // entity, weighted by tributary area (faces) or length (edges) --
    // not split evenly per node, which would depend on how fine the
    // mesh happens to be there. mesher must still be alive (not
    // shutdown()) since this queries it for the real triangles/segments
    // of each load's entity.
    static Eigen::VectorXd resolveForces(
        const Mesh&                mesh,
        const std::vector<Load>&  loads,
        const Assembler&           assembler,
        const Mesher&              mesher
    );
};
