#include "cli/InteractiveCLI.hpp"
#include "mesh/Mesher.hpp"
#include "fem/Assembler.hpp"
#include "fem/Solver.hpp"
#include "fem/LoadResolver.hpp"
#include "fem/StressResolver.hpp"
#include "io/ModelSerializer.hpp"
#include <fmt/core.h>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>

void InteractiveCLI::printBanner() {
    fmt::print("\n+==============================+\n");
    fmt::print(  "|        openFEM  v0.1         |\n");
    fmt::print(  "+==============================+\n\n");
}

void InteractiveCLI::printStatus() {
    fmt::print("\n--- Pipeline status ----------------------\n");
    fmt::print("  Geometry : {}{}\n", m_model.geometryLoaded  ? "[x]" : "[ ]",
               m_model.geometryLoaded  ? "" : "  -> type 'load'");
    fmt::print("  Physics  : {}{}\n", m_model.physicsAssigned ? "[x]" : "[ ]",
               m_model.physicsAssigned ? "" : "  -> type 'assign'");
    fmt::print("  Mesh     : {}{}\n", m_model.meshGenerated   ? "[x]" : "[ ]",
               m_model.meshGenerated   ? "" : "  -> runs automatically at the end of 'assign'");
    fmt::print("  Solved   : {}{}\n", m_model.solved          ? "[x]" : "[ ]",
               m_model.solved          ? "" : "  -> type 'bc', then 'loads', then 'solve'");
    fmt::print("------------------------------------------\n\n");
}

void InteractiveCLI::printMenu() {
    fmt::print("Available commands:\n");
    fmt::print("  load      -- load STEP file\n");
    fmt::print("  assign    -- assign elements/materials, then meshes automatically\n");
    fmt::print("  mesh      -- re-mesh with the current assignments (only needed to redo it)\n");
    fmt::print("  bc        -- define boundary conditions\n");
    fmt::print("  loads     -- define loads\n");
    fmt::print("  solve     -- assemble K and solve\n");
    fmt::print("  results   -- show results\n");
    fmt::print("  save      -- save model setup to JSON (to reopen/redo later)\n");
    fmt::print("  open      -- load model setup from JSON\n");
    fmt::print("  snapshot  -- save the current viewer state (mesh+results) under a name\n");
    fmt::print("  restore   -- bring back a saved snapshot in the viewer, no re-solving\n");
    fmt::print("  refresh   -- clear the viewer and reset the pipeline, start over\n");
    fmt::print("  status    -- show pipeline status\n");
    fmt::print("  quit      -- exit\n\n");
}

double InteractiveCLI::promptDouble(const std::string& msg, double def) {
    fmt::print("{} [{}] > ", msg, def);
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) return def;
    try { return std::stod(line); }
    catch (...) { return def; }
}

std::string InteractiveCLI::promptString(const std::string& msg,
                                          const std::string& def) {
    if (def.empty()) fmt::print("{} > ", msg);
    else             fmt::print("{} [{}] > ", msg, def);
    std::string line;
    std::getline(std::cin, line);
    return line.empty() ? def : line;
}

bool InteractiveCLI::promptConfirm(const std::string& msg) {
    std::string r = promptString(msg + " [y/N]", "n");
    return (r == "y" || r == "Y");
}

void InteractiveCLI::cmdLoad() {
    std::string path = promptString("STEP file path");
    if (path.empty()) return;

    if (!m_importer.load(path)) {
        fmt::print("  [ERROR] Cannot load: {}\n", path);
        return;
    }
    m_model.stepFile       = path;
    m_model.geometryLoaded = true;

    // Export to BREP and import into gmsh right away, so entities are
    // identified by gmsh's own (dim, tag) from this point on -- see the
    // pipeline notes at the top of Mesher.cpp.
    std::string brepPath = "/tmp/openFEM_model.brep";
    if (!m_importer.exportToBRep(brepPath)) {
        fmt::print("  [ERROR] Cannot export BREP.\n");
        return;
    }
    const auto& entities = m_mesher.importGeometry(brepPath);

    fmt::print("\n=== Geometry entities found ===\n");
    for (const auto& ent : entities) {
        const char* hint =
            ent.femHint == GeometryEntity::FEMHint::TETRA ? "TETRA" :
            ent.femHint == GeometryEntity::FEMHint::BEAM  ? "BEAM"  : "?";
        fmt::print("  [dim {} tag {:2d}] {:20s}  dims: {:6.1f} x {:6.1f} x {:6.1f} mm  -> suggested: {}\n",
                    ent.dim, ent.tag, ent.name,
                    (ent.bbox[1]-ent.bbox[0]),
                    (ent.bbox[3]-ent.bbox[2]),
                    (ent.bbox[5]-ent.bbox[4]),
                    hint);
    }

    autoExportWeb();
}

void InteractiveCLI::cmdAssign() {
    if (!m_model.geometryLoaded) {
        fmt::print("  Load a STEP file first with 'load'.\n");
        return;
    }

    // Re-running 'assign' redoes the assignment from scratch -- without
    // this, a second run APPENDS to the existing lists instead of
    // replacing them: meshAndExtract() would then process every
    // entity's elements once per duplicate assignment (silently
    // doubling them, and later double-counting their stiffness in
    // Assembler::assemble()), and the global mesh-size ceiling would
    // take the min across old AND new sizes, so a coarser re-request
    // would have no visible effect.
    m_model.assignments.clear();
    m_model.materials.clear();
    m_model.sections.clear();

    // E/yieldStrength in MPa, rho in tonne/mm^3 -- see the unit system
    // note in fem/Model.hpp. Aluminium's yield strength is just a
    // typical value for a generic structural alloy (6061-T6), not tied
    // to any specific spec -- use "custom" for a real material.
    std::map<std::string, Material> builtinMaterials = {
        {"C40",       {0, "C40",       206e3, 0.30, 7.85e-9, 420.0}},
        {"aluminium", {1, "aluminium",  69e3, 0.33, 2.70e-9, 270.0}},
        {"custom",    {2, "custom",        0,    0,       0,    0}},
    };

    fmt::print("\n=== Physics assignment ===\n");
    fmt::print("Available materials: C40, aluminium, custom\n\n");

    // Shell is intentionally not offered here -- choosing it wouldn't
    // do anything geometrically useful without real mid-surface
    // extraction behind it (see Mesher::importGeometry()), so it's
    // disabled rather than half-working. The underlying element
    // (ShellElement, ElementType::TRI3, Assembler's TRI3 case) is all
    // still there for when that's revisited -- this only hides it from
    // the CLI.
    for (const auto& ent : m_mesher.entities()) {
        const char* defHint =
            ent.femHint == GeometryEntity::FEMHint::BEAM ? "beam" : "tetra";

        fmt::print("--- Entity [dim {} tag {:2d}] {} ---\n", ent.dim, ent.tag, ent.name);

        std::string typeStr = promptString("  Element type (tetra/beam)", defHint);

        EntityAssignment a;
        a.dim      = ent.dim;
        a.entityId = ent.tag;
        a.meshSize = promptDouble("  Mesh size [mm]", 10.0);

        if (typeStr == "tetra") {
            a.femType = ElementType::TETRA4;
        } else if (typeStr == "beam") {
            a.femType = ElementType::BEAM2;

            // TODO: these 4 values are entered by hand below. A more
            // advanced version could compute them automatically from the
            // real cross-section geometry: area/Iy/Iz are a direct
            // integral over the section profile (OpenCASCADE's
            // BRepGProp can do this once we have that profile). J
            // (torsion constant) is harder -- except for simple shapes
            // (e.g. circular) it requires solving a Saint-Venant torsion
            // sub-problem, not a plain integral. Either way, we would
            // first need to extract the actual 2D cross-section from the
            // geometry (the STEP edge here is only the 1D centerline) --
            // a geometric extraction step similar in spirit (and
            // complexity) to the shell mid-surface extraction this
            // project doesn't do yet either (see Mesher::importGeometry()).
            BeamSection sec;
            sec.id   = (int)m_model.sections.size();
            sec.name = fmt::format("Section_{}", sec.id);
            sec.area = promptDouble("  Area A [mm^2]",            100.0);
            sec.Iy   = promptDouble("  Inertia Iy [mm^4]",      10000.0);
            sec.Iz   = promptDouble("  Inertia Iz [mm^4]",      10000.0);
            sec.J    = promptDouble("  Torsion constant J [mm^4]", 10000.0);
            m_model.sections.push_back(sec);
            a.sectionId = sec.id;
        }

        std::string matName = promptString("  Material", "C40");
        auto it = builtinMaterials.find(matName);
        if (it != builtinMaterials.end()) {
            Material mat = it->second;
            if (matName == "custom") {
                mat.E             = promptDouble("    E [MPa]",          210e3);
                mat.nu            = promptDouble("    nu",                 0.3);
                mat.rho           = promptDouble("    rho [tonne/mm^3]", 7.85e-9);
                mat.yieldStrength = promptDouble("    Yield strength [MPa]", 0.0);
            }
            mat.id = (int)m_model.materials.size();
            m_model.materials.push_back(mat);
            a.materialId = mat.id;
        }

        m_model.assignments.push_back(a);
        fmt::print("  OK assigned\n\n");
    }

    m_model.physicsAssigned = true;
    fmt::print("=== Assignment complete ===\n");

    // Meshing isn't adaptive to bc/loads in this project (size only
    // comes from the assignments above -- see meshAndExtract()), so
    // there's no technical reason to make the user run 'mesh' as a
    // separate step right after this one. 'mesh' is still available
    // on its own if they want to re-mesh with the same assignments.
    cmdMesh();

    // Only useful if you intend to reopen this setup later (via
    // 'open') to redo meshing/solving on it -- it saves the model
    // definition (materials, assignments, bc/loads), not the mesh or
    // any solved results.
    if (promptConfirm("Save configuration now (only needed if you want to reuse it later)?"))
        cmdSave();
}

void InteractiveCLI::cmdMesh() {
    if (!m_model.physicsAssigned) {
        fmt::print("  Assign physics first with 'assign'.\n");
        return;
    }

    // Geometry was already imported into gmsh back in cmdLoad() (that's
    // where importGeometry() ran) -- here we only set sizes, mesh, and
    // read the result back.
    m_model.mesh = m_mesher.meshAndExtract(m_model.assignments);
    m_model.meshGenerated = true;

    fmt::print("\n=== Mesh generated ===\n");
    fmt::print("  Nodes   : {}\n", m_model.mesh.nodes.size());
    fmt::print("  Elements: {}\n", m_model.mesh.elements.size());

    autoExportWeb();
}

void InteractiveCLI::autoExportWeb() {
    if (!m_model.geometryLoaded) return;

    ResultsExport results;
    results.displacements  = m_model.nodalDisplacements;
    results.stressNormal   = m_model.nodalStressNormal;
    results.vonMises       = m_model.nodalVonMises;
    results.yieldStrength  = m_model.yieldStrengthUsed;
    results.maxVonMisesRaw = m_model.maxVonMisesRaw;
    results.safetyFactor   = m_model.safetyFactor;

    m_mesher.exportWebView("/workspace/web", m_model.bcs, m_model.loads, results);
}

void InteractiveCLI::cmdBoundary() {
    if (!m_model.geometryLoaded) {
        fmt::print("  Load a STEP file first with 'load'.\n");
        return;
    }

    fmt::print("\n=== Boundary conditions ===\n");

    do {
        if (!m_model.bcs.empty()) {
            fmt::print("\nCurrently set:\n");
            for (std::size_t i = 0; i < m_model.bcs.size(); ++i) {
                const char* typeName = m_model.bcs[i].type == BCType::FIXED  ? "fixed"  :
                                       m_model.bcs[i].type == BCType::PINNED ? "pinned" : "custom";
                fmt::print("  [{}] dim {} tag {} -- {}\n",
                           i, m_model.bcs[i].dim, m_model.bcs[i].entityId, typeName);
            }
        }

        std::string action = promptString("  Action (add/remove)", "add");
        if (action == "remove") {
            if (m_model.bcs.empty()) {
                fmt::print("  Nothing to remove.\n");
                continue;
            }
            int idx = (int)promptDouble("  Index to remove", -1);
            if (idx < 0 || idx >= (int)m_model.bcs.size()) {
                fmt::print("  Unknown index.\n");
                continue;
            }
            m_model.bcs.erase(m_model.bcs.begin() + idx);
            fmt::print("  Removed.\n");
            continue;
        }

        // Unlike 'assign', boundary conditions target one specific
        // face/edge of a solid, not the solid as a whole -- so we use
        // the unfiltered entity list here (it includes a solid's own
        // faces/edges, not just the "free" top-level ones).
        auto boundaryEntities = m_mesher.allBoundaryEntities();

        fmt::print("\nAvailable entities:\n");
        for (const auto& ent : boundaryEntities)
            fmt::print("  [dim {} tag {:2d}] {:<10} center ({:.1f}, {:.1f}, {:.1f}) mm\n",
                       ent.dim, ent.tag, ent.name,
                       (ent.bbox[0] + ent.bbox[1]) / 2.0,
                       (ent.bbox[2] + ent.bbox[3]) / 2.0,
                       (ent.bbox[4] + ent.bbox[5]) / 2.0);

        int dim = (int)promptDouble("  Entity dim (3=solid/2=face/1=edge)", -1);
        int tag = (int)promptDouble("  Entity tag", -1);
        bool found = false;
        for (const auto& ent : boundaryEntities)
            if (ent.dim == dim && ent.tag == tag) { found = true; break; }
        if (!found) {
            fmt::print("  Unknown (dim, tag).\n");
            continue;
        }

        std::string typeStr = promptString(
            "  BC type (fixed/pinned/custom)", "fixed");

        BoundaryCondition bc;
        bc.dim      = dim;
        bc.entityId = tag;

        if (typeStr == "pinned") {
            bc.type = BCType::PINNED;
            // Translations locked, rotations free.
            bc.locked[0] = bc.locked[1] = bc.locked[2] = true;
            bc.locked[3] = bc.locked[4] = bc.locked[5] = false;
        } else if (typeStr == "custom") {
            bc.type = BCType::CUSTOM;
            bc.locked[0] = promptConfirm("    Lock ux?");
            bc.locked[1] = promptConfirm("    Lock uy?");
            bc.locked[2] = promptConfirm("    Lock uz?");
            bc.locked[3] = promptConfirm("    Lock rotx?");
            bc.locked[4] = promptConfirm("    Lock roty?");
            bc.locked[5] = promptConfirm("    Lock rotz?");
        } else {
            bc.type = BCType::FIXED;
            for (bool& l : bc.locked) l = true;
        }

        m_model.bcs.push_back(bc);
        fmt::print("  OK added\n");

    } while (promptConfirm("Make another change to boundary conditions?"));

    autoExportWeb();
}

void InteractiveCLI::cmdLoads() {
    if (!m_model.geometryLoaded) {
        fmt::print("  Load a STEP file first with 'load'.\n");
        return;
    }

    fmt::print("\n=== Loads ===\n");

    do {
        if (!m_model.loads.empty()) {
            fmt::print("\nCurrently set:\n");
            for (std::size_t i = 0; i < m_model.loads.size(); ++i) {
                const auto& ld = m_model.loads[i];
                const char* typeName = ld.type == LoadType::FORCE_NODAL ? "force" :
                                       ld.type == LoadType::PRESSURE   ? "pressure" : "gravity";
                fmt::print("  [{}] dim {} tag {} -- {}\n", i, ld.dim, ld.entityId, typeName);
            }
        }

        std::string action = promptString("  Action (add/remove)", "add");
        if (action == "remove") {
            if (m_model.loads.empty()) {
                fmt::print("  Nothing to remove.\n");
                continue;
            }
            int idx = (int)promptDouble("  Index to remove", -1);
            if (idx < 0 || idx >= (int)m_model.loads.size()) {
                fmt::print("  Unknown index.\n");
                continue;
            }
            m_model.loads.erase(m_model.loads.begin() + idx);
            fmt::print("  Removed.\n");
            continue;
        }

        // Same reasoning as cmdBoundary(): a load targets one specific
        // face/edge of a solid, not the solid as a whole, so we need
        // the unfiltered entity list here too.
        auto boundaryEntities = m_mesher.allBoundaryEntities();

        fmt::print("\nAvailable entities:\n");
        for (const auto& ent : boundaryEntities)
            fmt::print("  [dim {} tag {:2d}] {:<10} center ({:.1f}, {:.1f}, {:.1f}) mm\n",
                       ent.dim, ent.tag, ent.name,
                       (ent.bbox[0] + ent.bbox[1]) / 2.0,
                       (ent.bbox[2] + ent.bbox[3]) / 2.0,
                       (ent.bbox[4] + ent.bbox[5]) / 2.0);

        int dim = (int)promptDouble("  Entity dim (3=solid/2=face/1=edge)", -1);
        int tag = (int)promptDouble("  Entity tag", -1);
        bool found = false;
        for (const auto& ent : boundaryEntities)
            if (ent.dim == dim && ent.tag == tag) { found = true; break; }
        if (!found) {
            fmt::print("  Unknown (dim, tag).\n");
            continue;
        }

        // Note: gravity (LoadType::GRAVITY in fem/Model.hpp) is
        // deliberately not offered here. It's a body force acting on
        // the whole mass of the model, not something you attach to one
        // chosen entity like a force or a pressure -- it doesn't fit
        // this per-entity flow. Not computed automatically either: if
        // the user wants it, they apply it explicitly.
        std::string typeStr = promptString(
            "  Load type (force/pressure)", "force");

        Load ld;
        ld.dim      = dim;
        ld.entityId = tag;

        if (typeStr == "pressure") {
            ld.type = LoadType::PRESSURE;
            // Single magnitude, stored in values[0] -- values[1]/[2]
            // unused for this type.
            ld.values[0] = promptDouble("    Pressure [MPa]", 0.0);
        } else {
            ld.type = LoadType::FORCE_NODAL;
            ld.values[0] = promptDouble("    Fx [N]", 0.0);
            ld.values[1] = promptDouble("    Fy [N]", 0.0);
            ld.values[2] = promptDouble("    Fz [N]", 0.0);
        }

        m_model.loads.push_back(ld);
        fmt::print("  OK added\n");

    } while (promptConfirm("Make another change to loads?"));

    autoExportWeb();
}

void InteractiveCLI::cmdSolve() {
    if (!m_model.meshGenerated) {
        fmt::print("  Generate the mesh first with 'mesh'.\n");
        return;
    }
    if (m_model.bcs.empty()) {
        fmt::print("  No boundary conditions defined ('bc') -- the system "
                    "would be free to move as a rigid body and the solver "
                    "would fail. Add at least one before solving.\n");
        return;
    }

    Assembler assembler;
    auto K = assembler.assemble(m_model.mesh, m_model.materials, m_model.sections);

    auto fixedDOFs = LoadResolver::resolveFixedDOFs(m_model.bcs, assembler, m_mesher);
    auto f         = LoadResolver::resolveForces(m_model.mesh, m_model.loads, assembler, m_mesher);

    Solver solver;
    auto result = solver.solve(K, f, fixedDOFs);

    if (!result.converged) {
        fmt::print("  [ERROR] Solve failed: {}\n", result.message);
        return;
    }

    m_model.displacements = result.displacements;
    m_model.solved         = true;

    // Per-node [ux,uy,uz] view of the same data, for the results
    // viewer. The first 3 DOFs of every node are always translations
    // (ux,uy,uz), in that order, regardless of whether the node has
    // 3 or 6 total DOF -- every element type (Tetra/Beam/Shell) puts
    // them there in its local stiffness matrix, so the Assembler's
    // global layout inherits the same order.
    const auto& nodeOffset = assembler.nodeOffset();
    m_model.nodalDisplacements.resize(m_model.mesh.nodes.size());
    for (std::size_t i = 0; i < m_model.mesh.nodes.size(); ++i) {
        int off = nodeOffset[i];
        m_model.nodalDisplacements[i] = {
            result.displacements[off + 0],
            result.displacements[off + 1],
            result.displacements[off + 2]
        };
    }

    // Stress recovery (TETRA4 only) from the same displacement vector,
    // then averaged to nodes for a smooth contour -- see
    // StressResolver/fem/Model.hpp for why the true unaveraged peak
    // (maxVonMisesRaw/safetyFactor) is tracked separately.
    Eigen::VectorXd uVec = Eigen::Map<const Eigen::VectorXd>(
        result.displacements.data(), (Eigen::Index)result.displacements.size());
    auto elementStresses = StressResolver::computeElementStresses(
        m_model.mesh, m_model.materials, uVec, nodeOffset);
    auto nodalStress = StressResolver::averageToNodes(
        m_model.materials, elementStresses, (int)m_model.mesh.nodes.size());

    m_model.nodalStressNormal = nodalStress.nodalSigmaNormal;
    m_model.nodalVonMises     = nodalStress.nodalVonMises;
    m_model.maxVonMisesRaw    = nodalStress.maxVonMisesRaw;
    m_model.yieldStrengthUsed = nodalStress.yieldStrengthUsed;
    m_model.safetyFactor      = nodalStress.safetyFactor;

    fmt::print("  Solved -- max displacement: {:.6e}\n",
               *std::max_element(
                   result.displacements.begin(), result.displacements.end(),
                   [](double a, double b) { return std::abs(a) < std::abs(b); }
               ));
    fmt::print("  Peak von Mises stress (unaveraged): {:.3f} MPa\n", m_model.maxVonMisesRaw);
    if (m_model.safetyFactor > 0)
        fmt::print("  Safety factor (yield / peak von Mises): {:.2f}\n", m_model.safetyFactor);
    else
        fmt::print("  Safety factor: n/a (no yield strength set on the material)\n");

    autoExportWeb();
}

void InteractiveCLI::cmdResults() {
    // TODO: implement a text results summary (max/min displacement and
    // stress, safety factor) here -- cmdSolve() already fills
    // m_model.displacements/nodalVonMises/safetyFactor, and the 'web'
    // command's "Results"/"Stress" pages cover the visual side.
    fmt::print("  [TODO] Results not yet implemented.\n");
}

void InteractiveCLI::cmdSave() {
    std::string path = promptString("JSON file path", "model.json");
    ModelSerializer::save(m_model, path);
}

void InteractiveCLI::cmdLoadJson() {
    std::string path = promptString("JSON file path", "model.json");
    if (ModelSerializer::load(m_model, path))
        fmt::print("  Model loaded successfully.\n");
}

namespace { const char* SNAPSHOT_DIR = "/workspace/web/snapshots"; }

void InteractiveCLI::cmdSnapshot() {
    if (!m_model.geometryLoaded) {
        fmt::print("  Load a STEP file first with 'load'.\n");
        return;
    }
    // scene.json already has everything the viewer needs (mesh,
    // displacements, stress, safety factor) -- a snapshot is just a
    // copy of it under a name you choose, kept around so 'restore'
    // can bring that exact view back later without re-solving
    // anything. Re-export first so the snapshot reflects what's
    // actually on screen right now, not a stale prior state.
    autoExportWeb();

    std::string name = promptString("Snapshot name", "snapshot1");
    if (name.empty()) return;

    std::filesystem::create_directories(SNAPSHOT_DIR);
    std::ifstream in("/workspace/web/scene.json", std::ios::binary);
    if (!in) {
        fmt::print("  [ERROR] No scene.json to snapshot yet.\n");
        return;
    }
    std::string dst = std::string(SNAPSHOT_DIR) + "/" + name + ".json";
    std::ofstream out(dst, std::ios::binary);
    out << in.rdbuf();
    fmt::print("  Saved snapshot: {}\n", dst);
}

void InteractiveCLI::cmdRestore() {
    namespace fs = std::filesystem;
    fs::create_directories(SNAPSHOT_DIR);

    fmt::print("\nAvailable snapshots:\n");
    bool any = false;
    for (const auto& entry : fs::directory_iterator(SNAPSHOT_DIR)) {
        if (entry.path().extension() == ".json") {
            fmt::print("  {}\n", entry.path().stem().string());
            any = true;
        }
    }
    if (!any) {
        fmt::print("  (none yet -- use 'snapshot' to create one)\n");
        return;
    }

    std::string name = promptString("Snapshot name to restore", "");
    if (name.empty()) return;

    std::string src = std::string(SNAPSHOT_DIR) + "/" + name + ".json";
    std::ifstream in(src, std::ios::binary);
    if (!in) {
        fmt::print("  [ERROR] No such snapshot: {}\n", name);
        return;
    }
    std::ofstream out("/workspace/web/scene.json", std::ios::binary);
    out << in.rdbuf();
    // The browser polls scene.json on its own (see web/index.html) --
    // this just replaces the file, no separate "push" needed.
    fmt::print("  Restored '{}' -- the browser will pick it up within a "
               "couple seconds.\n", name);
}

void InteractiveCLI::cmdRefresh() {
    // Full reset: a fresh FEMModel (clears every pipeline flag and all
    // assignments/bcs/loads/results), gmsh shut down (so the next
    // 'load' re-initializes cleanly instead of reusing old state), and
    // scene.json deleted so the browser viewer (which polls for it)
    // goes blank instead of keeping the previous mesh on screen.
    m_model = FEMModel{};
    m_mesher.shutdown();
    std::filesystem::remove("/workspace/web/scene.json");
    fmt::print("  Cleared -- viewer and pipeline state reset. Use 'load' to start over.\n");
}

void InteractiveCLI::run() {
    printBanner();
    printMenu();

    while (true) {
        printStatus();
        std::string cmd = promptString("openFEM");

        if      (cmd == "load")   cmdLoad();
        else if (cmd == "assign") cmdAssign();
        else if (cmd == "mesh")   cmdMesh();
        else if (cmd == "bc")     cmdBoundary();
        else if (cmd == "loads")  cmdLoads();
        else if (cmd == "solve")  cmdSolve();
        else if (cmd == "results")cmdResults();
        else if (cmd == "save")   cmdSave();
        else if (cmd == "open")   cmdLoadJson();
        else if (cmd == "snapshot") cmdSnapshot();
        else if (cmd == "restore")  cmdRestore();
        else if (cmd == "refresh")  cmdRefresh();
        else if (cmd == "help")   printMenu();
        else if (cmd == "quit" || cmd == "exit") {
            m_mesher.shutdown();
            fmt::print("Goodbye.\n");
            break;
        }
        else fmt::print("  Unknown command. Type 'help'.\n");
    }
}