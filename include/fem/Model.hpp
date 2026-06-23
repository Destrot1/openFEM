#pragma once
#include <string>
#include <vector>
#include <array>
#include "mesh/Mesh.hpp"

// Unit system used everywhere: mm, N, MPa, tonne -- matches STEP files
// (authored in mm), since geometry is never converted to meters.
struct Material {
    int         id;
    std::string name;
    double      E             = 210e3;   // MPa (210 GPa)
    double      nu            = 0.3;     // Poisson coeff. 
    double      rho           = 7.85e-9; // tonne/mm^3
    double      yieldStrength = 0;       // MPa, 0 = not set (safety factor not computable)
};

struct BeamSection {
    int         id;
    std::string name;
    double      area = 0; // mm^2
    double      Iy   = 0; // mm^4
    double      Iz   = 0; // mm^4
    double      J    = 0; // mm^4
};

// entityId is a gmsh tag, unique only within its own dimension --
// dim is needed too, to identify an entity without ambiguity.
struct EntityAssignment {
    int         dim        = -1;
    int         entityId   = -1;
    ElementType femType    = ElementType::TETRA4;
    int         materialId = -1;
    double      thickness  = 0.0;  // mm
    int         sectionId  = -1;
    double      meshSize   = 10.0; // mm
};

enum class BCType { FIXED, PINNED, CUSTOM };

// locked[ux,uy,uz,rotx,roty,rotz] = true means that DOF is fixed to zero.
// FIXED locks all 6, PINNED locks translations only, CUSTOM picks per-DOF.
struct BoundaryCondition {
    int    dim      = -1;
    int    entityId = -1;
    BCType type     = BCType::FIXED;
    bool   locked[6] = {false,false,false,false,false,false};
};

enum class LoadType { FORCE_NODAL, PRESSURE, GRAVITY };

// FORCE_NODAL values: Fx, Fy, Fz in N.
// PRESSURE values[0]: MPa (= N/mm^2, consistent with the mm/MPa system).
struct Load {
    int      dim      = -1;
    int      entityId = -1;
    LoadType type     = LoadType::FORCE_NODAL;
    double   values[3] = {0,0,0};
};

struct FEMModel {
    std::string name;
    std::string stepFile;

    Mesh mesh;

    std::vector<Material>          materials;
    std::vector<BeamSection>       sections;
    std::vector<EntityAssignment>  assignments;
    std::vector<BoundaryCondition> bcs;
    std::vector<Load>              loads;

    bool geometryLoaded  = false;
    bool meshGenerated   = false;
    bool physicsAssigned = false;
    bool solved          = false;

    std::vector<double> displacements;

    // Per-node [ux,uy,uz] in mm, same order as mesh.nodes -- filled by
    // cmdSolve() after solving; empty until then.
    std::vector<std::array<double,3>> nodalDisplacements;

    // Per-node stress (MPa), averaged over touching elements for a smooth
    // contour. Filled by cmdSolve(); empty until then.
    std::vector<std::array<double,3>> nodalStressNormal; // sxx, syy, szz
    std::vector<double>               nodalVonMises;

    // True unaveraged peak von Mises stress and the resulting safety factor --
    // the nodal average above would hide a sharp local stress peak.
    double maxVonMisesRaw    = 0; // MPa
    double yieldStrengthUsed = 0; // yield strength of the material at the peak element, MPa
    double safetyFactor      = 0; // yieldStrengthUsed / maxVonMisesRaw, 0 if not computable
};