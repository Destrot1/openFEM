/*
 * ================================================================
 *  MESHING PIPELINE -- how a STEP file becomes an FEM-ready Mesh
 * ================================================================
 *
 *  1. StepImporter::load(stepFile)               [OpenCASCADE]
 *     Reads the STEP file and keeps the shape in memory (g_shape).
 *     Does NOT discover entities -- that now happens after the
 *     geometry reaches gmsh (step 3), so entities are identified
 *     directly by gmsh's own (dim, tag), with no separate numbering
 *     to keep in sync with anything else.
 *
 *  2. StepImporter::exportToBRep(path)            [OpenCASCADE]
 *     Writes the SAME in-memory shape (g_shape) to a .brep file.
 *     Why BREP and not just reusing the original STEP file?
 *     STEP is a generic exchange format: any reader has to
 *     RECONSTRUCT the topology (faces/edges/vertices and how they
 *     connect) from its text entities, and that reconstruction is
 *     not guaranteed identical between independent reads -- tolerances,
 *     sewing, etc. can differ. BREP is a direct, lossless dump of
 *     OpenCASCADE's own internal topological structure: writing it
 *     and reading it back reproduces the exact same shape, with the
 *     exact same face boundaries, no reconstruction, no ambiguity.
 *     IMPORTANT: at this point the .brep file is still PURE geometry.
 *     No tetrahedra, triangles or segments exist yet.
 *
 *  3. Mesher::importGeometry(brepFile)
 *     Imports the BREP into gmsh (gmsh::model::occ::importShapes).
 *     From here on, gmsh owns the geometry and assigns each entity a
 *     (dim, tag): dim=3 solids, dim=2 surfaces, dim=1 curves. We read
 *     that list directly (getEntities) and use gmsh's own (dim, tag)
 *     as THE identifier for everything downstream -- EntityAssignment
 *     / BoundaryCondition / Load all store dim + entityId=tag. No
 *     translation table needed, because nothing identifies entities
 *     independently before this point anymore.
 *     Still no mesh elements at this point -- only entities discovered.
 *
 *  4. CLI "assign"
 *     User picks femType + material + meshSize per entity (shown from
 *     importGeometry()'s result) -> builds the EntityAssignment list,
 *     reusing the same dim + tag already shown.
 *
 *  5. Mesher::meshAndExtract(assignments)
 *
 *     a) gmsh::model::mesh::setSize({{dim,tag}}, meshSize) per
 *        assignment -- sets the target element size (in mm -- see the
 *        unit system note in fem/Model.hpp) directly on that one
 *        entity, no ambiguity.
 *
 *     b) gmsh::model::mesh::generate(3), attempted in a child process
 *        with a time limit (meshWithTimeout), as a generic safety net
 *        against a mesh request that's slow/heavy for whatever reason
 *        (e.g. an accidentally tiny size relative to the part, or a
 *        genuinely large/complex model) -- not because of any known
 *        gmsh issue. THIS is the only place where tetrahedra/triangles/
 *        segments are actually created. gmsh decides the element type
 *        from each entity's GEOMETRIC DIMENSION, not from our femType:
 *            solid   (3D) -> tetrahedra
 *            surface (2D) -> triangles
 *            curve   (1D) -> line segments
 *        shell support is currently disabled in the CLI (see
 *        importGeometry()) -- a solid always becomes tetrahedra here)
 *
 *     c) Read back nodes (gmsh::model::mesh::getNodes, no dim/tag ->
 *        all nodes at once) -> fill our Node list.
 *
 *     d) Read back elements PER ASSIGNMENT (gmsh::model::mesh::
 *        getElements(..., dim, tag)) -> each element is tagged with
 *        THAT assignment's materialId/thickness/sectionId, never a
 *        different entity's by mistake.
 *
 *  6. Result: a Mesh (nodes + elements), ready for Assembler::assemble().
 *     BoundaryCondition/Load are resolved into actual DOFs/forces later,
 *     by LoadResolver (fem/LoadResolver.hpp) -- it queries Mesher
 *     directly (entityNodeIndices/surfaceTriangles/surfaceLines) for
 *     whichever entity each one targets, live, rather than going
 *     through a node map built ahead of time -- this works for ANY
 *     entity, not just the ones configured in 'assign'.
 * ================================================================
 */

#include "mesh/Mesher.hpp"
#include <gmsh.h>
#include <fmt/core.h>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <algorithm>
#include <map>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <csignal>
#include <cmath>

// Works like a dictionary: looks up a gmsh element type number and
// returns our own matching ElementType (e.g. gmsh's "4" -> our TETRA4).
ElementType Mesher::gmshTypeToFEM(int gmshType) {
    switch (gmshType) {
        case 4: return ElementType::TETRA4; // 4-node tetrahedron
        case 2: return ElementType::TRI3;   // 3-node triangle
        case 1: return ElementType::BEAM2;  // 2-node line
        default:
            throw std::runtime_error(
                fmt::format("Unknown gmsh element type: {}", gmshType)
            );
    }
}

std::vector<GeometryEntity> Mesher::allBoundaryEntities() const {
    std::vector<GeometryEntity> result;
    if (!m_gmshReady)
        return result;

    for (int dim = 3; dim >= 1; --dim) {
        std::vector<std::pair<int,int>> ents;
        gmsh::model::getEntities(ents, dim);

        for (const auto& [d, tag] : ents) {
            GeometryEntity ent;
            ent.dim = d;
            ent.tag = tag;
            ent.name = fmt::format(
                "{}_{}",
                dim == 3 ? "Solid" : dim == 2 ? "Face" : "Edge",
                tag
            );
            // Shell support is disabled in the CLI for now (see
            // importGeometry() for why) -- never hint at it here either.
            ent.femHint = dim == 3 ? GeometryEntity::FEMHint::TETRA :
                          dim == 1 ? GeometryEntity::FEMHint::BEAM :
                                     GeometryEntity::FEMHint::NONE;
            // Can throw ("Empty bounding box") for an entity with no
            // elements/points to measure -- this happens for example
            // after a sized mesh attempt succeeds and we gmsh::open()
            // its .msh result, which can leave a stray entity behind
            // with nothing in it. Skip that one entity rather than
            // letting the whole program crash on an unhandled exception.
            try {
                gmsh::model::getBoundingBox(
                    dim, tag,
                    ent.bbox[0], ent.bbox[2], ent.bbox[4],
                    ent.bbox[1], ent.bbox[3], ent.bbox[5]
                );
            } catch (const std::exception&) {
                continue;
            }
            result.push_back(ent);
        }
    }
    return result;
}

std::vector<std::array<int,3>> Mesher::surfaceTriangles(int dim, int tag) const {
    std::vector<std::array<int,3>> result;
    if (!m_gmshReady)
        return result;

    // Rebuild the gmsh-tag -> our-index map the same way meshAndExtract()
    // did. This works because gmsh returns nodes in the same order every
    // time you ask, as long as nothing got re-meshed in between -- true
    // here since we're called after meshAndExtract() and before shutdown().
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, params;
    gmsh::model::mesh::getNodes(nodeTags, coords, params);
    std::map<std::size_t, int> nodeIndexMap;
    for (std::size_t i = 0; i < nodeTags.size(); ++i)
        nodeIndexMap[nodeTags[i]] = (int)i;

    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

    for (std::size_t t = 0; t < elemTypes.size(); ++t) {
        if (elemTypes[t] != 2)   // triangles only
            continue;
        for (std::size_t e = 0; e < elemTags[t].size(); ++e) {
            result.push_back({
                nodeIndexMap[elemNodeTags[t][e*3 + 0]],
                nodeIndexMap[elemNodeTags[t][e*3 + 1]],
                nodeIndexMap[elemNodeTags[t][e*3 + 2]]
            });
        }
    }
    return result;
}

std::vector<std::array<int,2>> Mesher::surfaceLines(int dim, int tag) const {
    std::vector<std::array<int,2>> result;
    if (!m_gmshReady)
        return result;

    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, params;
    gmsh::model::mesh::getNodes(nodeTags, coords, params);
    std::map<std::size_t, int> nodeIndexMap;
    for (std::size_t i = 0; i < nodeTags.size(); ++i)
        nodeIndexMap[nodeTags[i]] = (int)i;

    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

    for (std::size_t t = 0; t < elemTypes.size(); ++t) {
        if (elemTypes[t] != 1)   // line segments only
            continue;
        for (std::size_t e = 0; e < elemTags[t].size(); ++e) {
            result.push_back({
                nodeIndexMap[elemNodeTags[t][e*2 + 0]],
                nodeIndexMap[elemNodeTags[t][e*2 + 1]]
            });
        }
    }
    return result;
}

bool Mesher::exportWebView(const std::string& dir,
                           const std::vector<BoundaryCondition>& bcs,
                           const std::vector<Load>& loads,
                           const ResultsExport& results) {
    if (!m_gmshReady) {
        fmt::print("[Mesher] Nothing to show yet -- load a geometry first.\n");
        return false;
    }

    // Always export whatever gmsh currently has. If no mesh exists yet
    // (the user hasn't run 'mesh'), make a quick 2D-only preview
    // (surfaces only, not the full 3D volume) just so there is
    // something to render -- this is NOT the real mesh, that one is
    // still made later by meshAndExtract() from the user's chosen
    // sizes. If a real mesh already exists, this block is skipped and
    // we export that instead.
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags);
    if (elemTypes.empty()) {
        fmt::print("[Mesher] No mesh yet -- generating a quick preview...\n");

        // Size it off the model's own bounding box -- this runs right
        // after 'load', before 'assign'/'mesh' have given us anything
        // to go on, so there's no per-entity size yet to honor. Clamped
        // on both ends so it stays fast on a huge model and doesn't
        // degenerate into thousands of tiny triangles on a small one.
        double xmin, ymin, zmin, xmax, ymax, zmax;
        gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);
        double diag = std::sqrt((xmax-xmin)*(xmax-xmin) +
                                 (ymax-ymin)*(ymax-ymin) +
                                 (zmax-zmin)*(zmax-zmin));
        double previewSize = std::clamp(diag / 40.0, 0.5, 20.0);
        gmsh::option::setNumber("Mesh.MeshSizeMax", previewSize);

        gmsh::model::mesh::generate(2);

        // Don't let this leak into a later meshAndExtract() call --
        // it sets its own ceiling from the user's assignments, but
        // only when at least one assignment requests one.
        gmsh::option::setNumber("Mesh.MeshSizeMax", 1e22);

        elemTypes.clear(); elemTags.clear(); elemNodeTags.clear();
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags);
    }

    // We write our own simple JSON instead of using gmsh's file writer:
    // this build's writer doesn't support OBJ at all, and STL would
    // lose node sharing and only handle triangles -- no tets, no
    // entity/label info, nothing we could reuse later for results.

    // -- All nodes, in gmsh's own order ------------------------------
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

    std::map<std::size_t, int> nodeIndexMap;
    for (std::size_t i = 0; i < nodeTags.size(); ++i)
        nodeIndexMap[nodeTags[i]] = (int)i;

    nlohmann::json nodesJson = nlohmann::json::array();
    for (std::size_t i = 0; i < nodeTags.size(); ++i)
        nodesJson.push_back({nodeCoords[3*i], nodeCoords[3*i+1], nodeCoords[3*i+2]});

    // -- Triangles/lines PER ENTITY, queried directly from gmsh -------
    // When gmsh meshes a solid (generate(3)), it first meshes the
    // surfaces bounding it (generate(2) internally) and those surface
    // elements stay in the mesh, each still tagged to the real (dim=2)
    // surface it belongs to -- so we don't need to guess the solid's
    // outer "skin" by counting shared tet faces; we just ask gmsh for
    // each surface's own elements directly, already correctly tagged.
    // This also means each face/edge in the viewer can be highlighted
    // on its own when hovered (see web/index.html) -- the user can
    // point at the actual surface instead of hunting for a tiny dot.
    nlohmann::json entitiesJson = nlohmann::json::array();
    for (const auto& ent : allBoundaryEntities()) {
        if (ent.dim == 3)
            continue;   // a solid has no surface mesh of its own to draw

        std::vector<int> et;
        std::vector<std::vector<std::size_t>> etags, enodes;
        gmsh::model::mesh::getElements(et, etags, enodes, ent.dim, ent.tag);

        nlohmann::json tris = nlohmann::json::array();
        nlohmann::json lns  = nlohmann::json::array();

        for (std::size_t t = 0; t < et.size(); ++t) {
            int nodesPerElem = (et[t] == 2) ? 3 : (et[t] == 1) ? 2 : 0;
            if (nodesPerElem == 0)
                continue;

            for (std::size_t e = 0; e < etags[t].size(); ++e) {
                nlohmann::json idx = nlohmann::json::array();
                for (int n = 0; n < nodesPerElem; ++n) {
                    std::size_t gmshNodeTag = enodes[t][e * nodesPerElem + n];
                    idx.push_back(nodeIndexMap[gmshNodeTag]);
                }
                if (nodesPerElem == 3) tris.push_back(idx);
                else                   lns.push_back(idx);
            }
        }

        if (tris.empty() && lns.empty())
            continue;   // nothing meshed for this entity yet

        entitiesJson.push_back({
            {"dim",  ent.dim},
            {"tag",  ent.tag},
            {"name", ent.name},
            {"x", (ent.bbox[0] + ent.bbox[1]) / 2.0},
            {"y", (ent.bbox[2] + ent.bbox[3]) / 2.0},
            {"z", (ent.bbox[4] + ent.bbox[5]) / 2.0},
            {"triangles", tris},
            {"lines", lns}
        });
    }

    // -- bcs/loads, purely for the viewer to mark them and build a
    // legend -- this does NOT feed into solving, just display. -------
    nlohmann::json bcsJson = nlohmann::json::array();
    for (const auto& bc : bcs) {
        double xmin=0,ymin=0,zmin=0,xmax=0,ymax=0,zmax=0;
        try {
            gmsh::model::getBoundingBox(bc.dim, bc.entityId, xmin,ymin,zmin,xmax,ymax,zmax);
        } catch (const std::exception&) {
            continue; // entity gone/empty -- see allBoundaryEntities() comment
        }
        const char* typeName = bc.type == BCType::FIXED  ? "fixed" :
                                bc.type == BCType::PINNED ? "pinned" :
                                                             "custom";
        bcsJson.push_back({
            {"dim", bc.dim}, {"tag", bc.entityId}, {"type", typeName},
            {"locked", std::vector<bool>(bc.locked, bc.locked + 6)},
            {"x", (xmin+xmax)/2.0}, {"y", (ymin+ymax)/2.0}, {"z", (zmin+zmax)/2.0}
        });
    }

    nlohmann::json loadsJson = nlohmann::json::array();
    for (const auto& ld : loads) {
        double xmin=0,ymin=0,zmin=0,xmax=0,ymax=0,zmax=0;
        try {
            gmsh::model::getBoundingBox(ld.dim, ld.entityId, xmin,ymin,zmin,xmax,ymax,zmax);
        } catch (const std::exception&) {
            continue; // entity gone/empty -- see allBoundaryEntities() comment
        }
        const char* typeName = ld.type == LoadType::FORCE_NODAL ? "force" :
                                ld.type == LoadType::PRESSURE   ? "pressure" :
                                                                   "gravity";
        loadsJson.push_back({
            {"dim", ld.dim}, {"tag", ld.entityId}, {"type", typeName},
            {"values", std::vector<double>(ld.values, ld.values + 3)},
            {"x", (xmin+xmax)/2.0}, {"y", (ymin+ymax)/2.0}, {"z", (zmin+zmax)/2.0}
        });
    }

    // Per-node displacements (ux,uy,uz in mm), for the results viewer
    // to compute the deformed shape and color by field. Only included
    // if the count matches the nodes we just exported -- e.g. empty
    // before 'solve' has run, or if the mesh changed since then.
    nlohmann::json dispJson = nlohmann::json::array();
    if (results.displacements.size() == nodeTags.size()) {
        for (const auto& d : results.displacements)
            dispJson.push_back({d[0], d[1], d[2]});
    }

    // Same idea for the per-node stress fields -- nodally averaged
    // (see ResultsExport/StressResolver), so they're smooth contours
    // here; maxVonMisesRaw/safetyFactor below are the unaveraged true
    // peak instead, for the actual safety check.
    nlohmann::json stressJson;
    if (results.stressNormal.size() == nodeTags.size() &&
        results.vonMises.size()     == nodeTags.size()) {
        nlohmann::json sigmaJson = nlohmann::json::array();
        for (const auto& s : results.stressNormal)
            sigmaJson.push_back({s[0], s[1], s[2]});
        stressJson["sigma"]          = sigmaJson;
        stressJson["vonMises"]       = results.vonMises;
        stressJson["yieldStrength"]  = results.yieldStrength;
        stressJson["maxVonMisesRaw"] = results.maxVonMisesRaw;
        stressJson["safetyFactor"]   = results.safetyFactor;
    }

    nlohmann::json out;
    out["nodes"]         = nodesJson;
    out["entities"]      = entitiesJson;
    out["bcs"]           = bcsJson;
    out["loads"]         = loadsJson;
    out["displacements"] = dispJson;
    out["stress"]        = stressJson;

    std::ofstream f(dir + "/scene.json");
    f << out.dump();

    fmt::print("[Mesher] Web view exported to: {}\n", dir);
    return true;
}

const std::vector<GeometryEntity>&
Mesher::importGeometry(const std::string& brepFile)
{
    fmt::print("[Mesher] Loading geometry: {}\n", brepFile);
    m_brepPath = brepFile;

    if (!m_gmshReady) {
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", 2);
        m_gmshReady = true;
    } else {
        // A geometry was already loaded in this same gmsh session --
        // clear it first so entities from the previous load don't pile
        // up together with the new ones.
        gmsh::clear();
    }

    gmsh::vectorpair outDimTags;
    gmsh::model::occ::importShapes(brepFile, outDimTags);
    gmsh::model::occ::synchronize();

    m_entities.clear();

    // Solids first, then surfaces, then curves -- same display order
    // as the old StepImporter-based discovery.
    for (int dim = 3; dim >= 1; --dim) {
        std::vector<std::pair<int,int>> ents;
        gmsh::model::getEntities(ents, dim);

        for (const auto& [d, tag] : ents) {
            // Keep only entities NOT already part of a higher-dimension
            // one -- the gmsh equivalent of the old "free face" / "free
            // edge" check (which used OpenCASCADE's TopExp_Explorer). A
            // solid's own boundary faces/edges get meshed automatically
            // as part of meshing the solid; they are not independent
            // entities the user should assign separately.
            if (dim < 3) {
                std::vector<int> upward, downward;
                gmsh::model::getAdjacencies(dim, tag, upward, downward);
                if (!upward.empty())
                    continue;
            }

            GeometryEntity ent;
            ent.dim = dim;
            ent.tag = tag;
            ent.name = fmt::format(
                "{}_{}",
                dim == 3 ? "Solid" : dim == 2 ? "Face" : "Edge",
                tag
            );
            // Shell support is disabled in the CLI for now (see
            // importGeometry() for why) -- never hint at it here either.
            ent.femHint = dim == 3 ? GeometryEntity::FEMHint::TETRA :
                          dim == 1 ? GeometryEntity::FEMHint::BEAM :
                                     GeometryEntity::FEMHint::NONE;

            // gmsh's own mesher does NOT need this bounding box for
            // anything -- it meshes straight from the topology/geometry
            // and the mesh size we set, regardless of bbox. We ask gmsh
            // for it ourselves, purely to show approximate dimensions to
            // the user -- gmsh's own mesher doesn't need it.
            //
            // (A "thin solid -> suggest shell" heuristic used to live
            // here. Shell support needs real mid-surface extraction to
            // mean anything geometrically -- a hint without that
            // capability behind it was more confusing than useful, so
            // it's disabled for now rather than half-implemented. See
            // cmdAssign() in InteractiveCLI.cpp.)
            gmsh::model::getBoundingBox(
                dim, tag,
                ent.bbox[0], ent.bbox[2], ent.bbox[4],
                ent.bbox[1], ent.bbox[3], ent.bbox[5]
            );

            m_entities.push_back(ent);
        }
    }

    fmt::print("[Mesher] {} entities found.\n", m_entities.size());
    return m_entities;
}

std::vector<int> Mesher::entityNodeIndices(int dim, int tag) const {
    std::vector<int> result;
    if (!m_gmshReady)
        return result;

    // Same tag -> our-index rebuild as surfaceTriangles()/surfaceLines().
    std::vector<std::size_t> allNodeTags;
    std::vector<double> allCoords, allParams;
    gmsh::model::mesh::getNodes(allNodeTags, allCoords, allParams);
    std::map<std::size_t, int> nodeIndexMap;
    for (std::size_t i = 0; i < allNodeTags.size(); ++i)
        nodeIndexMap[allNodeTags[i]] = (int)i;

    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, params;
    // includeBoundary=true: by default gmsh excludes nodes that also
    // belong to this entity's boundary (e.g. nodes shared with the
    // edges around a face). A small/oddly-shaped face can end up with
    // ALL its nodes on that boundary and none "exclusively its own" --
    // without this flag, getNodes() would then return an empty list
    // even though the face clearly has real triangles (this is exactly
    // what happened in testing: a face with real mesh elements still
    // came back with zero nodes here, silently skipping its BC/load).
    gmsh::model::mesh::getNodes(nodeTags, coords, params, dim, tag, true);

    result.reserve(nodeTags.size());
    for (auto t : nodeTags)
        result.push_back(nodeIndexMap.at(t));
    return result;
}

bool Mesher::meshWithTimeout(int timeoutSeconds, const std::string& tmpMshPath) {
    pid_t pid = fork();

    if (pid < 0) {
        // Couldn't fork at all -- just mesh directly, no timeout
        // protection possible, better than refusing to work.
        fmt::print("[Mesher] fork() failed -- meshing without a timeout.\n");
        gmsh::model::mesh::generate(3);
        gmsh::write(tmpMshPath);
        return true;
    }

    if (pid == 0) {
        // Child: this is a separate copy of the whole process (and of
        // gmsh's in-memory state) from the moment of fork() -- nothing
        // it does here can corrupt the parent. Whatever happens
        // (success, crash, hang killed by the parent), the parent's
        // own gmsh model is untouched until it explicitly opens our
        // result file.
        try {
            gmsh::model::mesh::generate(3);
            gmsh::write(tmpMshPath);
            _exit(0);
        } catch (...) {
            _exit(1);
        }
    }

    // Parent: poll instead of a single blocking wait, so we can give
    // up and kill the child after timeoutSeconds.
    for (int i = 0; i < timeoutSeconds * 10; ++i) {
        int status = 0;
        pid_t r = waitpid(pid, &status, WNOHANG);
        if (r == pid)
            return WIFEXITED(status) && WEXITSTATUS(status) == 0;
        usleep(100000); // 100ms
    }

    kill(pid, SIGKILL);
    waitpid(pid, nullptr, 0);  // reap the zombie
    return false;
}

Mesh Mesher::meshAndExtract(const std::vector<EntityAssignment>& assignments)
{
    // Always start from the ORIGINAL CAD geometry, not whatever gmsh's
    // model currently has -- a previous successful sized mesh replaced
    // it with a mesh-only reconstruction (gmsh::open() below), which
    // has no CAD topology left to re-mesh at a different size. Without
    // this, calling meshAndExtract() again (e.g. re-running 'assign'
    // with a different mesh size) would silently keep the old mesh.
    if (!m_brepPath.empty())
        importGeometry(m_brepPath);

    // -- Mesh size: a global ceiling, gmsh refines freely below it --
    // The smallest size requested across all assignments becomes a
    // global maximum (Mesh.MeshSizeMax, in mm), and gmsh's own default
    // heuristic (curvature-based refinement near small features,
    // coarser elsewhere) stays free to refine further wherever it
    // judges necessary. This is NOT yet true per-entity control --
    // every entity shares the same ceiling -- that needs a Distance/
    // Threshold field per entity, which is future work (see also the
    // idea of refining near bc/loads entities, same mechanism).
    double maxSize = 0.0;
    for (const auto& a : assignments)
        if (maxSize == 0.0 || a.meshSize < maxSize)
            maxSize = a.meshSize;
    if (maxSize > 0.0)
        gmsh::option::setNumber("Mesh.MeshSizeMax", maxSize);

    // -- Generate the mesh, with a safety net -----------------------
    // generate(3) = full 3D mesh
    // automatically generates: points -> curves -> surfaces -> volumes
    //
    // Try it with the user's requested size (maxSize set above) first,
    // but bounded by a timeout, as a generic guard against a request
    // that's much heavier than intended (e.g. a size far too small for
    // the part) taking unreasonably long. If it doesn't finish in time,
    // fall back to gmsh's own unconstrained default sizing -- coarser/
    // finer than asked, but it returns a usable mesh instead of nothing.
    bool sized = meshWithTimeout(20, "/tmp/openFEM_sized_mesh.msh");
    if (sized) {
        gmsh::open("/tmp/openFEM_sized_mesh.msh");
    } else {
        fmt::print("[Mesher] Sized meshing didn't finish in time -- "
                    "falling back to gmsh's own default mesh "
                    "(requested size not honored for this run).\n");
        gmsh::option::setNumber("Mesh.MeshSizeMax", 1e22);
        gmsh::model::mesh::generate(3);
    }

    // Optimize element quality
    gmsh::model::mesh::optimize("Laplace2D");

    // -- Read the nodes --------------------------------------------
    Mesh result;

    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;

    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

    // Map gmsh tag -> our index (gmsh does not guarantee 0-based indices)
    std::map<std::size_t, int> nodeIndexMap;
    result.nodes.reserve(nodeTags.size());

    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        Node n;
        n.id = (int)i;
        n.position = Eigen::Vector3d(
            nodeCoords[3*i],
            nodeCoords[3*i + 1],
            nodeCoords[3*i + 2]
        );
        nodeIndexMap[nodeTags[i]] = (int)i;
        result.nodes.push_back(n);
    }

    // -- Read the elements, one entity (one assignment) at a time --
    // Querying getElements() per (dim, tag) -- instead of once for the
    // whole model -- means each element is tagged with the material/
    // thickness/section of the entity it actually came from, instead
    // of always assignments[0].
    int elemId = 0;
    for (const auto& a : assignments) {
        std::vector<int>                      elemTypes;
        std::vector<std::vector<std::size_t>> elemTags;
        std::vector<std::vector<std::size_t>> elemNodeTags;

        gmsh::model::mesh::getElements(
            elemTypes, elemTags, elemNodeTags, a.dim, a.entityId
        );

        for (std::size_t t = 0; t < elemTypes.size(); ++t) {
            int gmshType = elemTypes[t];

            // Accept only the types we handle
            if (gmshType != 1 && gmshType != 2 && gmshType != 4)
                continue;

            ElementType femType;
            int nodesPerElem;
            try {
                femType = gmshTypeToFEM(gmshType);
                nodesPerElem = (gmshType == 4) ? 4 :
                               (gmshType == 2) ? 3 : 2;
            } catch (...) { continue; }

            for (std::size_t e = 0; e < elemTags[t].size(); ++e) {
                MeshElement elem;
                elem.id   = elemId++;
                elem.type = femType;

                for (int n = 0; n < nodesPerElem; ++n) {
                    std::size_t gmshNodeTag =
                        elemNodeTags[t][e * nodesPerElem + n];
                    elem.nodeIds.push_back(nodeIndexMap[gmshNodeTag]);
                }

                elem.materialId = a.materialId;
                elem.thickness  = a.thickness;
                elem.sectionId  = a.sectionId;

                result.elements.push_back(elem);
            }
        }
    }

    fmt::print("[Mesher] Nodes: {}  Elements: {}\n",
               result.nodes.size(), result.elements.size());

    // gmsh is intentionally left initialized here (not finalized) so
    // the generated mesh can still be viewed with showGUI() afterwards.
    return result;
}

void Mesher::showGUI() {
    if (!m_gmshReady) {
        fmt::print("[Mesher] Nothing to show yet -- load a geometry first.\n");
        return;
    }
    fmt::print("[Mesher] Opening gmsh viewer (close its window to return here)...\n");

    // By default gmsh only draws the wireframe (edges) of the geometry,
    // not the solid surfaces -- these two options turn surface shading
    // on, so the model looks like a solid instead of just lines.
    gmsh::option::setNumber("Geometry.Surfaces", 1);
    gmsh::option::setNumber("Geometry.SurfaceType", 2);

    gmsh::fltk::run();
}

void Mesher::shutdown() {
    if (m_gmshReady) {
        gmsh::finalize();
        m_gmshReady = false;
    }
}
