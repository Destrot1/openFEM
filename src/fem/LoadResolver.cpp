#include "fem/LoadResolver.hpp"
#include <fmt/core.h>
#include <map>

std::vector<int> LoadResolver::resolveFixedDOFs(
    const std::vector<BoundaryCondition>& bcs,
    const Assembler&                       assembler,
    const Mesher&                          mesher)
{
    std::vector<int> fixed;
    const auto& nodeOffset = assembler.nodeOffset();
    const auto& nodeDOF    = assembler.nodeDOF();

    for (const auto& bc : bcs) {
        auto nodeIndices = mesher.entityNodeIndices(bc.dim, bc.entityId);
        if (nodeIndices.empty()) {
            fmt::print("[LoadResolver] WARNING: no mesh nodes found for "
                       "BC on (dim={}, tag={}) -- skipped.\n", bc.dim, bc.entityId);
            continue;
        }

        for (int nodeIdx : nodeIndices) {
            int offset = nodeOffset[nodeIdx];
            int ndof   = nodeDOF[nodeIdx];
            // Only consider DOFs this node actually has -- a node with
            // 3 DOF (pure tetra corner) has no rotations to lock.
            for (int d = 0; d < ndof && d < 6; ++d)
                if (bc.locked[d])
                    fixed.push_back(offset + d);
        }
    }
    return fixed;
}

Eigen::VectorXd LoadResolver::resolveForces(
    const Mesh&               mesh,
    const std::vector<Load>& loads,
    const Assembler&          assembler,
    const Mesher&             mesher)
{
    Eigen::VectorXd f = Eigen::VectorXd::Zero(assembler.totalDOF());
    const auto& nodeOffset = assembler.nodeOffset();

    for (const auto& ld : loads) {
        auto nodeIndices = mesher.entityNodeIndices(ld.dim, ld.entityId);
        if (nodeIndices.empty()) {
            fmt::print("[LoadResolver] WARNING: no mesh nodes found for "
                       "load on (dim={}, tag={}) -- skipped.\n", ld.dim, ld.entityId);
            continue;
        }

        // Tributary weight per node: area share (faces) or length
        // share (edges) -- NOT an equal split, so the result doesn't
        // depend on how fine the mesh happens to be in that spot (see
        // project discussion: a node touched by more/bigger elements
        // represents more of the surface, and should carry more of
        // the load).
        std::map<int, double> weight;
        for (int n : nodeIndices) weight[n] = 0.0;

        Eigen::Vector3d normalSum = Eigen::Vector3d::Zero();

        if (ld.dim == 2) {
            for (const auto& tri : mesher.surfaceTriangles(ld.dim, ld.entityId)) {
                const Eigen::Vector3d& p0 = mesh.nodes[tri[0]].position;
                const Eigen::Vector3d& p1 = mesh.nodes[tri[1]].position;
                const Eigen::Vector3d& p2 = mesh.nodes[tri[2]].position;
                // Cross product magnitude = 2x triangle area -- using
                // it directly (not normalized) makes bigger triangles
                // contribute proportionally more to the average normal
                // below, which is what we want.
                Eigen::Vector3d cross = (p1 - p0).cross(p2 - p0);
                double area = 0.5 * cross.norm();
                weight[tri[0]] += area / 3.0;
                weight[tri[1]] += area / 3.0;
                weight[tri[2]] += area / 3.0;
                normalSum += cross;
            }
        } else if (ld.dim == 1) {
            for (const auto& seg : mesher.surfaceLines(ld.dim, ld.entityId)) {
                const Eigen::Vector3d& p0 = mesh.nodes[seg[0]].position;
                const Eigen::Vector3d& p1 = mesh.nodes[seg[1]].position;
                double len = (p1 - p0).norm();
                weight[seg[0]] += len / 2.0;
                weight[seg[1]] += len / 2.0;
            }
        }
        // dim == 3 (a whole solid) or no elements found above: no
        // tributary measure available -- fall back to an equal split
        // below (every weight stays 0, handled right after this).

        double totalWeight = 0.0;
        for (const auto& [n, w] : weight) totalWeight += w;
        if (totalWeight <= 0.0) {
            for (auto& [n, w] : weight) w = 1.0;
            totalWeight = (double)weight.size();
        }

        Eigen::Vector3d totalForce;
        if (ld.type == LoadType::PRESSURE) {
            // TODO: the average normal's SIGN depends on gmsh's
            // triangle winding order, which isn't guaranteed to point
            // outward from the solid -- it may need to be flipped.
            // Verify the direction in the viewer before trusting the
            // sign of a pressure load.
            Eigen::Vector3d normal = normalSum.normalized();
            double totalArea = totalWeight; // sum of the per-node area shares = total area
            totalForce = normal * ld.values[0] * totalArea;
        } else { // FORCE_NODAL
            totalForce = Eigen::Vector3d(ld.values[0], ld.values[1], ld.values[2]);
        }

        for (const auto& [n, w] : weight) {
            Eigen::Vector3d nodalForce = totalForce * (w / totalWeight);
            int offset = nodeOffset[n];
            f(offset + 0) += nodalForce.x();
            f(offset + 1) += nodalForce.y();
            f(offset + 2) += nodalForce.z();
        }
    }
    return f;
}
