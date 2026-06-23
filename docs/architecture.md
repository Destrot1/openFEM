# Architecture

## Directory layout

```
openFEM/
├── CMakeLists.txt         openFEM + test_beam_frame targets
├── docker-compose.yml     fem_dev service: bind-mounts the repo + home dir,
│                          ports 7681 (terminal), 8081 (web viewer)
├── docker/
│   ├── Dockerfile         Ubuntu 24.04 + gmsh/OpenCASCADE/Eigen/ttyd/etc.
│   └── start-gui.sh       Container entrypoint -- builds the project,
│                          serves the viewer, starts the embedded terminal
├── web/
│   ├── index.html         The viewer (Three.js)
│   └── scene.json         Generated at runtime, not committed
├── include/               Public headers, mirrors src/ by subsystem
│   ├── cli/
│   ├── elements/
│   ├── fem/
│   ├── geometry/
│   ├── io/
│   └── mesh/
├── src/                   Implementation, same subsystem split
│   ├── main.cpp
│   ├── cli/
│   ├── elements/
│   ├── fem/
│   ├── geometry/
│   ├── io/
│   ├── mesh/
│   └── tests/             Standalone manual check (own executable)
└── docs/                  This documentation
```

Two CMake targets:
- **`openFEM`** — the full CLI application.
- **`test_beam_frame`** — a standalone executable (elements + fem sources
  only, no CLI/mesh/geometry) that builds a small frame by hand and prints
  its displacements, to check the element/assembly/solver math in isolation.

## Pipeline

```
STEP file → geometry entities → assign element type/material/mesh size
  → mesh (gmsh) → boundary conditions + loads
  → solve (assemble K, factorize, back-substitute) → results
```

Every step that changes the model re-exports its state to the web viewer
on its own (`autoExportWeb()` in the CLI) — the browser polls for changes,
so there's no manual "refresh" step.

## Subsystem-by-subsystem

### `geometry/StepImporter`
Reads a STEP file via OpenCASCADE and converts it to BREP — a direct,
lossless dump of OpenCASCADE's own topology, with no reconstruction
ambiguity (unlike re-reading the STEP text format, which different readers
can reconstruct slightly differently).

### `mesh/Mesher`
Owns the whole gmsh-facing side of the pipeline, in two phases:
- **`importGeometry()`** imports the BREP and discovers its entities, using
  gmsh's own `(dim, tag)` as the identifier everywhere downstream — no
  separate numbering to keep in sync.
- **`meshAndExtract()`** applies the requested mesh size, generates the
  mesh (with a timeout-guarded fallback to gmsh's own default sizing for
  requests that turn out too heavy), and reads the result back as a `Mesh`.

It also exports the live view for the browser (`exportWebView()`) and
answers live entity queries (which nodes/triangles belong to a given face,
used to resolve boundary conditions and loads).

### `fem/Model.hpp`
The central data model: `Material`, `BeamSection`, `EntityAssignment`,
`BoundaryCondition`, `Load`, and `FEMModel` (the whole state of one
analysis: mesh, assignments, results). Unit system used everywhere: **mm,
N, MPa, tonne** — matching STEP files, which are authored in mm and never
converted.

### `elements/`
One class per element type, each implementing `stiffnessMatrix()`:
- **`TetraElement`** (4 nodes, 3 DOF/node) — constant-strain tetrahedron.
- **`ShellElement`** (3 nodes, 6 DOF/node) — CST membrane + bending.
- **`BeamElement`** (2 nodes, 6 DOF/node) — Euler-Bernoulli beam.

### `fem/Assembler`
Builds the global sparse stiffness matrix: maps each node to its global
DOFs, then scatters every element's local stiffness matrix into the right
place.

### `fem/LoadResolver`
Bridges entity-scoped boundary conditions/loads (which only know "this
applies to face X") into the DOF-level language the solver speaks (global
indices, a force vector). Loads are distributed across an entity's nodes
by tributary area/length, not evenly — so the result reflects the real
load distribution rather than an artifact of mesh density.

### `fem/Solver`
Direct sparse solve (`Eigen::SparseLU`) with the penalty method for
boundary conditions. Factorization runs on a background thread so the CLI
can print progress for large systems.

### `fem/StressResolver`
Recovers per-element stress from the solved displacements, then averages
it onto nodes for a smooth contour — while keeping the true (unaveraged)
peak stress and the resulting safety factor against the material's yield
strength, since averaging would otherwise hide a sharp stress
concentration.

### `io/ModelSerializer`
Saves/loads the model **setup** (materials, assignments, bc, loads) to/from
JSON — used by the CLI's `save`/`open`. Mesh and results aren't part of
this round-trip; see `snapshot`/`restore` in [commands.md](commands.md) for
saving a *solved* result instead.

### `cli/InteractiveCLI`
The REPL that drives everything — see [commands.md](commands.md) for the
full command list.

## Web viewer (`web/index.html`)

A single page, Three.js, no build step. Three modes:

- **Setup** — hover a face/edge to see its name; constraints/loads marked
  on the geometry; a small axis gizmo (mirrors the camera's rotation, not
  anchored to the model) in the corner.
- **Deformation** — exaggerated deformed shape, Ux/Uy/Uz/Total in a 4-way
  split, one shared camera.
- **Stress** — same deformed shape, σx/σy/σz/Von Mises, plus the true
  peak stress and safety factor next to the Von Mises legend.

A real terminal (via [ttyd](https://github.com/tsl0922/ttyd)) is docked at
the bottom of the page, always running the CLI.
