#pragma once
#include "mesh/Mesh.hpp"
#include "fem/Model.hpp"
#include <string>
#include <vector>
#include <map>
#include <array>

// A geometry entity as gmsh sees it, discovered by Mesher::importGeometry().
// (dim, tag) is the same identifier gmsh itself uses, and the same one
// EntityAssignment/BoundaryCondition/Load store -- no translation needed
// anywhere in the pipeline.
struct GeometryEntity {
    int dim = -1;   // 3 = solid, 2 = surface, 1 = curve
    int tag = -1;   // unique only within its own dimension

    std::string name;

    enum class FEMHint { TETRA, SHELL, BEAM, NONE } femHint = FEMHint::NONE;

    double bbox[6] = {0,0,0,0,0,0};
};

// Bundles everything the web viewer's "Results"/"Stress" pages need --
// grouped together since exportWebView()'s parameter list kept growing
// every time a new result field was added after a solve. All empty by
// default (e.g. before 'solve' has run); exportWebView() omits a field
// from scene.json entirely when its size doesn't match the node count.
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
    // Phase A: import the BREP file into gmsh and discover its entities
    // (dim, tag, bounding box, suggested femHint). Leaves gmsh
    // initialized, ready for meshAndExtract().
    const std::vector<GeometryEntity>& importGeometry(const std::string& brepFile);

    // Phase B: must be called after importGeometry(). Re-imports the
    // ORIGINAL BREP first (see m_brepPath) -- a successful sized mesh
    // attempt replaces gmsh's in-memory model with one reconstructed
    // from the mesh alone (gmsh::open() of the .msh result, see the
    // .cpp), which has no CAD topology left to re-mesh at a different
    // size. Without this re-import, calling meshAndExtract() a second
    // time (e.g. re-running 'assign' with a different mesh size) would
    // silently keep the previous mesh instead of regenerating it.
    // Applies the mesh size from each assignment, generates the mesh,
    // and reads it back as our own Mesh. Leaves gmsh initialized (does
    // NOT finalize it), so the result can still be viewed with
    // showGUI() afterwards.
    Mesh meshAndExtract(const std::vector<EntityAssignment>& assignments);

    // Opens gmsh's own graphical viewer and blocks until its window is
    // closed. Shows whatever gmsh currently has loaded -- just the
    // geometry if called after importGeometry(), the actual mesh if
    // called after meshAndExtract(). Useful to see exactly which
    // (dim, tag) is which entity, instead of guessing from bbox numbers.
    void showGUI();

    // Call once, when the whole program is done with gmsh (e.g. on CLI
    // "quit"). Safe to call even if gmsh was never initialized.
    void shutdown();

    // Exports whatever gmsh currently has (geometry or mesh) as a small
    // web-based view: dir/scene.json -- our own format (nodes,
    // triangles, lines, entity labels with their center of mass), NOT
    // one of gmsh's file formats (this build's writer doesn't support
    // OBJ, and STL would lose node sharing/element types/entity info).
    // If no mesh exists yet, generates a quick 2D-only preview mesh
    // first, just so there is something to show -- this is NOT the
    // real mesh that meshAndExtract will later produce from the user's
    // chosen sizes.
    // bcs/loads (can be empty) are included in the export so the
    // viewer can mark which entities have a constraint/load and show
    // a legend -- purely for display, this does not feed into solving.
    // results (can be default-constructed, e.g. before 'solve' has
    // run) carries displacements/stress/von Mises/safety factor, one
    // entry per node, same order as the exported nodes.
    // Returns false if there is nothing to export yet.
    bool exportWebView(const std::string& dir,
                       const std::vector<BoundaryCondition>& bcs,
                       const std::vector<Load>& loads,
                       const ResultsExport& results);

    const std::vector<GeometryEntity>& entities() const { return m_entities; }

    // Every surface/curve/solid currently in the model, INCLUDING ones
    // that belong to a higher-dimension entity (e.g. a solid's own
    // faces) -- unlike entities(), nothing is filtered out here. Use
    // this for boundary conditions/loads: a BC or a load is applied to
    // one specific face/edge of a solid, not to the solid as a whole,
    // so the "free entities only" list used for assigning element
    // types isn't enough there.
    std::vector<GeometryEntity> allBoundaryEntities() const;

    // The real triangles/line segments of one entity, queried live from
    // gmsh (must still be initialized -- i.e. called after
    // meshAndExtract(), before shutdown()), translated into OUR OWN
    // node indices (matching the Mesh returned by meshAndExtract()).
    // Used by LoadResolver for tributary-area/length weighting of
    // loads -- a face's nodes alone aren't enough for that, we also
    // need to know which nodes each triangle/segment actually connects.
    std::vector<std::array<int,3>> surfaceTriangles(int dim, int tag) const;
    std::vector<std::array<int,2>> surfaceLines(int dim, int tag) const;

    // Just the node indices (no connectivity) belonging to one entity,
    // queried live the same way as above. This works for ANY entity --
    // not just the ones configured in 'assign' -- which is what
    // LoadResolver needs for boundary conditions (it only needs to know
    // which nodes to lock, not how they connect).
    std::vector<int> entityNodeIndices(int dim, int tag) const;

private:
    std::vector<GeometryEntity> m_entities;
    bool m_gmshReady = false;

    // Remembered from the last importGeometry() call, so
    // meshAndExtract() can re-import the original CAD geometry fresh
    // every time it runs, instead of meshing whatever gmsh's model
    // happens to currently be (which may have been replaced by a
    // previous mesh attempt -- see meshAndExtract()).
    std::string m_brepPath;

    // Convert a gmsh element type into our ElementType
    ElementType gmshTypeToFEM(int gmshType);

    // Attempts gmsh::model::mesh::generate(3) in a CHILD PROCESS, with
    // a time limit. gmsh has no safe way to interrupt a mesh generation
    // call mid-computation from within the same process (it could leave
    // its internal state corrupted) -- a separate process is the only
    // reliable way to abort a stuck attempt without risking the whole
    // program. On success, writes the result to tmpMshPath and returns
    // true; the caller still has to gmsh::open(tmpMshPath) to load it
    // back into the parent's own model. Returns false (parent's model
    // untouched, still just the geometry, no file written) if the
    // child timed out or failed, so the caller can fall back to an
    // unconstrained mesh instead.
    bool meshWithTimeout(int timeoutSeconds, const std::string& tmpMshPath);
};
