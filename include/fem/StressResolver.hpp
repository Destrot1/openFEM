#pragma once
#include "mesh/Mesh.hpp"
#include "fem/Model.hpp"
#include <Eigen/Dense>
#include <vector>
#include <array>

// Constant stress per element, TETRA4 only -- shell/beam stress
// recovery isn't needed anywhere in this project yet, so it isn't
// implemented.
struct ElementStress {
    int               elementId;
    std::vector<int>  nodeIds;
    int               materialId;
    std::array<double,6> sigma; // sxx, syy, szz, sxy, syz, szx (MPa)
    double            vonMises; // MPa
};

struct NodalStressResult {
    // One entry per mesh node, averaged over every element touching
    // it -- gives a smooth contour for display. Averaging hides the
    // real peak at a sharp stress concentration (e.g. a hole's edge),
    // which is why maxVonMisesRaw/safetyFactor below come from the
    // unaveraged per-element values instead.
    std::vector<std::array<double,3>> nodalSigmaNormal; // sxx, syy, szz
    std::vector<double>               nodalVonMises;

    double maxVonMisesRaw    = 0; // true peak across elements, MPa
    double yieldStrengthUsed = 0; // yield strength of the material at the peak element, MPa
    double safetyFactor      = 0; // yieldStrengthUsed / maxVonMisesRaw, 0 if not computable
};

class StressResolver {
public:
    // sigma = D*B*ue per element, using the already-solved global
    // displacement vector and the Assembler's node->DOF offset map to
    // pull out each element's own nodal displacements.
    static std::vector<ElementStress> computeElementStresses(
        const Mesh&                  mesh,
        const std::vector<Material>& materials,
        const Eigen::VectorXd&       displacements,
        const std::vector<int>&      nodeOffset
    );

    static NodalStressResult averageToNodes(
        const std::vector<Material>&      materials,
        const std::vector<ElementStress>& elementStresses,
        int                                nodeCount
    );
};
