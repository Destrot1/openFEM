# Commands

Two separate sets of commands: **Docker** (run from your own terminal, to
start/stop the container) and the **openFEM CLI** (run inside the embedded
terminal at the bottom of `http://localhost:8081/`).

## Docker

Run these from the repository root.

| Command | What it does |
|---|---|
| `docker compose up -d --build` | First start on a new machine, or after changing `Dockerfile`/`start-gui.sh`/`docker-compose.yml`. Builds the image, then starts the container. |
| `docker compose up -d` | Normal start. If the image already exists, this is enough — no need for `--build` every time. |
| `docker compose down` | Stops and removes the container (the build cache volume is kept, so the compiled binary survives). |
| `docker compose down -v` | Full reset — also deletes the named volume holding the compiled C++ project. The next start recompiles from scratch. |
| `docker compose logs -f` | Follow the container's logs (useful if the embedded terminal/viewer don't come up). |

The container's entrypoint (`docker/start-gui.sh`) automatically:
- compiles the C++ project on every start (incremental — fast if nothing
  changed, only slow the very first time),
- serves the web viewer on `:8081`,
- serves the embedded terminal (already running `./openFEM`) on `:7681`.

## openFEM CLI

Typed inside the embedded terminal (or `./openFEM` from `/workspace/build`
in any shell inside the container).

| Command | What it does |
|---|---|
| `load` | Prompts for a STEP file path (e.g. `/Desktop/part.step`), imports the geometry, lists the detected entities. |
| `assign` | Per entity: element type (`tetra`/`beam`), mesh size [mm], material (`C40`/`aluminium`/`custom`). Then **meshes automatically** — no separate step needed. |
| `mesh` | Re-meshes with the current assignments. Only needed to redo it (e.g. after a sized attempt timed out and fell back to gmsh's default sizing). |
| `bc` | Add or remove boundary conditions. `add`: pick an entity (dim/tag) and a type (`fixed`/`pinned`/`custom`, the last lets you lock each of the 6 DOFs individually). `remove`: pick an index from the list shown. |
| `loads` | Add or remove loads, same `add`/`remove` flow as `bc`. `add`: `force` (Fx/Fy/Fz in N) or `pressure` (MPa). |
| `solve` | Assembles the stiffness matrix and solves. Prints the peak von Mises stress and the resulting safety factor (if a yield strength was set on the material). |
| `results` | Not implemented yet — use the web viewer's Deformation/Stress pages instead. |
| `save` | Saves the model **setup** (materials, assignments, bc, loads) to a JSON file — for reopening later with `open`. Does *not* save the mesh or any solved results. |
| `open` | Loads a model setup previously saved with `save`. |
| `snapshot` | Saves the *current viewer state* (mesh + results) under a name, in `web/snapshots/`. |
| `restore` | Lists saved snapshots and brings one back into the viewer — no re-solving needed. |
| `refresh` | Clears the pipeline and the viewer, so you can start a new analysis from a blank state without restarting the container. |
| `status` | Prints the pipeline status (which steps are done, and what to type next for any that aren't). |
| `help` | Prints the command list. |
| `quit` | Exits the CLI (drops the embedded terminal to a plain shell). |

## A full session, start to finish

```
load
/Desktop/my_part.step
assign
tetra
5
C40
n
bc
add
2
3
fixed
n
loads
add
2
7
force
0
0
-5000
n
solve
```

Then open the **Deformation**/**Stress** tabs in the viewer at
`http://localhost:8081/` to see the results.
