#pragma once
#include <string>
#include <vector>
#include <array>
#include "mesh/Mesh.hpp"

// Unit system used EVERYWHERE in this project: mm, N, MPa, tonne.
// Why: STEP files from mechanical CAD are almost always authored in
// millimeters, and we never convert the geometry (StepImporter/Mesher
// read node coordinates exactly as the file has them) -- so every other
// quantity (material stiffness, section properties, mesh size) has to
// be expressed in a unit system that's consistent with millimeters,
// not meters. MPa = N/mm^2 makes this work out exactly: 1 MPa = 1e6 Pa,
// matching the mm/m^3 scale difference. This is the same convention
// many CAD-driven FEM tools (Abaqus, Ansys, etc.) default to for this
// exact reason.
// NOTE: a STEP file authored in a different unit (inches, meters...)
// is NOT auto-detected or converted -- this project assumes mm.
struct Material {
    int         id;
    std::string name;
    double      E             = 210e3;   // MPa (210 GPa)
    double      nu            = 0.3;
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

// entityId is a gmsh tag, unique only within its own dimension -- so
// dim is required alongside it to identify an entity unambiguously
// (e.g. (dim=2, entityId=5) and (dim=3, entityId=5) are different entities).
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

// locked[i] = true means that DOF is constrained to zero; false means
// it's free. Order: ux, uy, uz, rotx, roty, rotz. FIXED/PINNED are just
// convenient presets for this same mask (see cmdBoundary): FIXED locks
// all 6, PINNED locks only the 3 translations. CUSTOM lets the user
// pick each one individually -- no magnitude, only locked-or-free
// (a nonzero prescribed displacement isn't supported, just lock/free).
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

    // One [ux,uy,uz] per node (mm), same order as mesh.nodes -- a
    // convenience view of `displacements` (which is flat, indexed by
    // global DOF) for anything that wants "this node moved by this
    // much", like the results viewer. Filled by cmdSolve() right after
    // solving; empty until then.
    std::vector<std::array<double,3>> nodalDisplacements;

    // Per-node stress (MPa), averaged over every element touching that
    // node -- smooths the contour for display, same idea as
    // nodalDisplacements above. Filled by cmdSolve(); empty until then.
    std::vector<std::array<double,3>> nodalStressNormal; // sxx, syy, szz
    std::vector<double>               nodalVonMises;

    // The TRUE (unaveraged) peak von Mises stress across all elements,
    // and the resulting safety factor against the yield strength of
    // the material at that element -- the nodal average above is only
    // for the contour plot; it would understate a sharp stress
    // concentration (e.g. at a hole's edge).
    double maxVonMisesRaw    = 0; // MPa
    double yieldStrengthUsed = 0; // yield strength of the material at the peak element, MPa
    double safetyFactor      = 0; // yieldStrengthUsed / maxVonMisesRaw, 0 if not computable
};