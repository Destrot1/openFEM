#pragma once
#include "mesh/Mesh.hpp"
#include "fem/Model.hpp"
#include "elements/TetraElement.hpp"
#include "elements/ShellElement.hpp"
#include "elements/BeamElement.hpp"
#include <Eigen/Sparse>
#include <vector>

class Assembler {
public:
    Eigen::SparseMatrix<double> assemble(
        const Mesh&                            mesh,
        const std::vector<Material>&           materials,
        const std::vector<BeamSection>&        sections
    );

    int totalDOF() const { return m_totalDOF; }

    // Node -> global DOF offset map (built during assemble()).
    // Needed to translate node-level BCs/loads into global DOF indices.
    const std::vector<int>& nodeOffset() const { return m_nodeOffset; }
    const std::vector<int>& nodeDOF()    const { return m_nodeDOF; }

private:
    int              m_totalDOF = 0;
    std::vector<int> m_nodeDOF;
    std::vector<int> m_nodeOffset;

    void computeDOFMap(const Mesh& mesh);

    void assembleElement(
        const MeshElement&                   elem,
        const Eigen::MatrixXd&               Ke,
        std::vector<Eigen::Triplet<double>>& triplets
    );
};