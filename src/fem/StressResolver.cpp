#include "fem/StressResolver.hpp"
#include "elements/TetraElement.hpp"
#include <cmath>

namespace {
double vonMisesFrom(const std::array<double,6>& s) {
    double sxx=s[0], syy=s[1], szz=s[2], sxy=s[3], syz=s[4], szx=s[5];
    return std::sqrt(
        0.5 * ((sxx-syy)*(sxx-syy) + (syy-szz)*(syy-szz) + (szz-sxx)*(szz-sxx))
        + 3.0 * (sxy*sxy + syz*syz + szx*szx)
    );
}
}

std::vector<ElementStress> StressResolver::computeElementStresses(
    const Mesh&                  mesh,
    const std::vector<Material>& materials,
    const Eigen::VectorXd&       displacements,
    const std::vector<int>&      nodeOffset)
{
    std::vector<ElementStress> result;
    TetraElement tetra;

    for (const auto& elem : mesh.elements) {
        if (elem.type != ElementType::TETRA4) continue;
        if (elem.materialId < 0 || elem.materialId >= (int)materials.size()) continue;
        const Material& mat = materials[elem.materialId];

        std::vector<Eigen::Vector3d> nodePos;
        nodePos.reserve(elem.nodeIds.size());
        for (int nid : elem.nodeIds) nodePos.push_back(mesh.nodes[nid].position);

        Eigen::VectorXd ue(12);
        for (int i = 0; i < 4; ++i) {
            int off = nodeOffset[elem.nodeIds[i]];
            ue(i*3 + 0) = displacements(off + 0);
            ue(i*3 + 1) = displacements(off + 1);
            ue(i*3 + 2) = displacements(off + 2);
        }

        Eigen::VectorXd sigma6 = tetra.stress(nodePos, mat.E, mat.nu, ue);

        ElementStress es;
        es.elementId  = elem.id;
        es.nodeIds    = elem.nodeIds;
        es.materialId = elem.materialId;
        for (int i = 0; i < 6; ++i) es.sigma[i] = sigma6(i);
        es.vonMises = vonMisesFrom(es.sigma);
        result.push_back(es);
    }
    return result;
}

NodalStressResult StressResolver::averageToNodes(
    const std::vector<Material>&      materials,
    const std::vector<ElementStress>& elementStresses,
    int                                nodeCount)
{
    NodalStressResult result;
    result.nodalSigmaNormal.assign(nodeCount, {0.0, 0.0, 0.0});
    result.nodalVonMises.assign(nodeCount, 0.0);
    std::vector<int> touchCount(nodeCount, 0);

    const ElementStress* peak = nullptr;
    for (const auto& es : elementStresses) {
        for (int nid : es.nodeIds) {
            result.nodalSigmaNormal[nid][0] += es.sigma[0];
            result.nodalSigmaNormal[nid][1] += es.sigma[1];
            result.nodalSigmaNormal[nid][2] += es.sigma[2];
            result.nodalVonMises[nid]       += es.vonMises;
            touchCount[nid]++;
        }
        if (!peak || es.vonMises > peak->vonMises) peak = &es;
    }

    for (int i = 0; i < nodeCount; ++i) {
        if (touchCount[i] == 0) continue;
        result.nodalSigmaNormal[i][0] /= touchCount[i];
        result.nodalSigmaNormal[i][1] /= touchCount[i];
        result.nodalSigmaNormal[i][2] /= touchCount[i];
        result.nodalVonMises[i]       /= touchCount[i];
    }

    // Safety factor uses the material of whichever element actually
    // produced the peak stress -- with one material in the model this
    // is moot, but it stays correct for a mixed-material assembly too.
    if (peak) {
        result.maxVonMisesRaw = peak->vonMises;
        if (peak->materialId >= 0 && peak->materialId < (int)materials.size()) {
            result.yieldStrengthUsed = materials[peak->materialId].yieldStrength;
            if (result.yieldStrengthUsed > 0)
                result.safetyFactor = result.yieldStrengthUsed / result.maxVonMisesRaw;
        }
    }

    return result;
}
