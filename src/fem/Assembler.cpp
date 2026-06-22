#include "fem/Assembler.hpp"
#include <fmt/core.h>
#include <stdexcept>

// -- Compute how many DOFs each node has ---------------------------
// A node shared between a tetra (3 DOF) and a beam (6 DOF)
// gets 6 DOF -- we always take the maximum
void Assembler::computeDOFMap(const Mesh& mesh) {
    int n = (int)mesh.nodes.size();
    m_nodeDOF.assign(n, 0);

    for (const auto& elem : mesh.elements) {
        int dof = (elem.type == ElementType::TETRA4) ? 3 : 6;
        for (int nid : elem.nodeIds) {
            if (nid < 0 || nid >= n)
                throw std::runtime_error("Invalid node ID in element");
            m_nodeDOF[nid] = std::max(m_nodeDOF[nid], dof);
        }
    }

    // Nodes not connected to any element -> 3 DOF by default
    for (int i = 0; i < n; ++i)
        if (m_nodeDOF[i] == 0) m_nodeDOF[i] = 3;

    // Cumulative offset
    m_nodeOffset.resize(n + 1);
    m_nodeOffset[0] = 0;
    for (int i = 0; i < n; ++i)
        m_nodeOffset[i+1] = m_nodeOffset[i] + m_nodeDOF[i];

    m_totalDOF = m_nodeOffset[n];
    fmt::print("[Assembler] Nodes: {}  Total DOF: {}\n", n, m_totalDOF);
}

// -- Insert Ke into the triplet list -------------------------------
// The triplet list is the standard way to build sparse matrices
// in Eigen: accumulate (row, column, value) entries and then call
// setFromTriplets once at the end
void Assembler::assembleElement(
    const MeshElement&                   elem,
    const Eigen::MatrixXd&               Ke,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    int dofPerNode = (elem.type == ElementType::TETRA4) ? 3 : 6;
    int nNodes     = (int)elem.nodeIds.size();

    for (int i = 0; i < nNodes; ++i) {
        int nodeI  = elem.nodeIds[i];
        int offsetI = m_nodeOffset[nodeI];

        for (int j = 0; j < nNodes; ++j) {
            int nodeJ   = elem.nodeIds[j];
            int offsetJ = m_nodeOffset[nodeJ];

            // Insert the dofPerNode x dofPerNode block
            for (int di = 0; di < dofPerNode; ++di)
                for (int dj = 0; dj < dofPerNode; ++dj)
                    triplets.emplace_back(
                        offsetI + di,
                        offsetJ + dj,
                        Ke(i * dofPerNode + di,
                           j * dofPerNode + dj)
                    );
        }
    }
}

// -- Main assembly -------------------------------------------------
Eigen::SparseMatrix<double> Assembler::assemble(
    const Mesh&                     mesh,
    const std::vector<Material>&    materials,
    const std::vector<BeamSection>& sections)
{
    computeDOFMap(mesh);

    // FEM elements
    TetraElement tetra;
    ShellElement shell;

    // Triplet list -- collects all contributions
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.elements.size() * 144);  // conservative estimate

    int processed = 0;
    int skipped   = 0;

    for (const auto& elem : mesh.elements) {

        // Get the node positions of the element
        std::vector<Eigen::Vector3d> nodePos;
        nodePos.reserve(elem.nodeIds.size());
        for (int nid : elem.nodeIds) {
            if (nid < 0 || nid >= (int)mesh.nodes.size())
                throw std::runtime_error("Node ID out of range");
            nodePos.push_back(mesh.nodes[nid].position);
        }

        // Get the material
        if (elem.materialId < 0 ||
            elem.materialId >= (int)materials.size()) {
            ++skipped;
            continue;
        }
        const Material& mat = materials[elem.materialId];

        // Compute Ke based on element type
        Eigen::MatrixXd Ke;
        try {
            switch (elem.type) {

                case ElementType::TETRA4:
                    Ke = tetra.stiffnessMatrix(nodePos, mat.E, mat.nu);
                    break;

                case ElementType::TRI3:
                    Ke = shell.stiffnessMatrix(
                             nodePos, mat.E, mat.nu, elem.thickness);
                    break;

                case ElementType::BEAM2: {
                    if (elem.sectionId < 0 ||
                        elem.sectionId >= (int)sections.size()) {
                        ++skipped;
                        continue;
                    }
                    const BeamSection& sec = sections[elem.sectionId];
                    BeamElement::Section bs {
                        sec.area, sec.Iy, sec.Iz, sec.J
                    };
                    BeamElement beam(bs);
                    Ke = beam.stiffnessMatrix(nodePos, mat.E, mat.nu);
                    break;
                }

                default:
                    ++skipped;
                    continue;
            }
        }
        catch (const std::exception& e) {
            fmt::print("[Assembler] Warning -- element {}: {}\n",
                       elem.id, e.what());
            ++skipped;
            continue;
        }

        assembleElement(elem, Ke, triplets);
        ++processed;
    }

    fmt::print("[Assembler] Processed: {}  Skipped: {}\n",
               processed, skipped);

    // Build the sparse matrix from the triplet list
    Eigen::SparseMatrix<double> K(m_totalDOF, m_totalDOF);
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();

    fmt::print("[Assembler] K size: {}x{}  Non-zeros: {}\n",
               K.rows(), K.cols(), K.nonZeros());

    return K;
}