#pragma once
#include "mesh/Mesh.hpp"
#include "fem/Model.hpp"
#include <string>
#include <vector>
#include <map>
#include <array>

// A geometry entity as gmsh sees it -- (dim, tag) is the same ID gmsh
// itself uses, and the same one EntityAssignment/BC/Load store.
struct GeometryEntity {
    int dim = -1;   // 3 = solid, 2 = surface, 1 = curve
    int tag = -1;   // unique only within its own dimension

    std::string name;

    enum class FEMHint { TETRA, SHELL, BEAM, NONE } femHint = FEMHint::NONE;

    double bbox[6] = {0,0,0,0,0,0};
};

// Bundles everything the web viewer's results pages need, grouped since
// exportWebView()'s parameter list kept growing. Empty before 'solve'.
struct ResultsExport {
    std::vector<std::array<double,3>> displacements;  // ux, uy, uz per node (mm)
    std::vector<std::array<double,3>> stressNormal;    // sxx, syy, szz per node (MPa, nodally averaged)
    std::vector<double>               vonMises;        // per node (MPa, nodally averaged)
    double yieldStrength  = 0; // MPa, yield strength of the material at the peak element
    double maxVonMisesRaw = 0; // true unaveraged peak across elements, MPa
    double safetyFactor   = 0; // yieldStrength / maxVonMisesRaw, 0 if not computable
};

class Mesher {
public:
    // Phase A: imports the BREP into gmsh and discovers its entities
    // (dim, tag, bbox, femHint). Leaves gmsh ready for meshAndExtract().
    const std::vector<GeometryEntity>& importGeometry(const std::string& brepFile);

    // Phase B: re-imports the original BREP first (a previous mesh attempt
    // replaces gmsh's model), then meshes at the assignments' size and reads it back.
    Mesh meshAndExtract(const std::vector<EntityAssignment>& assignments);

    // Call once when done with gmsh (e.g. on 'quit'). Safe even if gmsh
    // was never initialized.
    void shutdown();

    // Exports gmsh's current geometry/mesh as scene.json for the web viewer,
    // plus bcs/loads/results for display -- false if there's nothing to export yet.
    bool exportWebView(const std::string& dir,
                       const std::vector<BoundaryCondition>& bcs,
                       const std::vector<Load>& loads,
                       const ResultsExport& results);

    const std::vector<GeometryEntity>& entities() const { return m_entities; }

    // Every entity in the model, unfiltered (unlike entities()) -- needed
    // for bc/loads, since those target one face/edge, not a whole solid.
    std::vector<GeometryEntity> allBoundaryEntities() const;

    // Live triangles/line segments of one entity, in our own node indices --
    // used by LoadResolver for tributary-area/length weighting of loads.
    std::vector<std::array<int,3>> surfaceTriangles(int dim, int tag) const;
    std::vector<std::array<int,2>> surfaceLines(int dim, int tag) const;

    // Just the node indices of one entity, no connectivity -- works for ANY
    // entity, used by LoadResolver for boundary conditions.
    std::vector<int> entityNodeIndices(int dim, int tag) const;

private:
    std::vector<GeometryEntity> m_entities;
    bool m_gmshReady = false;

    // Remembered from the last importGeometry() call, so meshAndExtract()
    // can always re-import the original CAD geometry fresh.
    std::string m_brepPath;

    // Convert a gmsh element type into our ElementType
    ElementType gmshTypeToFEM(int gmshType);

    // Tries meshing in a child process with a time limit (gmsh has no safe
    // way to interrupt it) -- false on timeout/failure, true + a .msh file on success.
    bool meshWithTimeout(int timeoutSeconds, const std::string& tmpMshPath);
};
