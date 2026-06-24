# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is Elle

Elle is a C/C++ microstructural modelling framework for simulating physical processes (grain growth, diffusion, recrystallisation, melting, etc.) in 2D polycrystalline materials. It uses an imake/X11 build system.

## Building

Elle uses `imake` to generate `Makefile`s from `Imakefile`s. The top-level build script:

```sh
./install.sh [batch | wx | x | all]
```

- `batch` — headless executables, no display (→ `binb/`)
- `wx` — wxWidgets/OpenGL display executables (→ `binwx/`)
- `x` — X11/Motif display executables (→ `binx/`)
- `all` — builds all three variants

`bin/` is symlinked to the last selected variant. The script runs `xmkmf`, `make Makefiles`, `make install_base`, then the appropriate install targets (`install_b`, `install_wx`, `install_x`).

### Rebuilding after code changes (without a full reinstall)

From the relevant subdirectory (e.g., `processes/diffusion/`):

```sh
make depend    # regenerate dependencies
make           # compile
make install   # copy binary to the bin directory
```

To regenerate all Makefiles from Imakefiles at the top level:

```sh
xmkmf && make Makefiles
```

### Build variants per process

Each process `Imakefile` uses macros defined in `extra.rules` and `elle.tmpl`:

- `MakeBatch(name, ext, objs, destdir)` — headless build
- `MakeDisplayWX(name, ext, objs, destdir)` — wxWidgets build
- `MakeDisplayX(name, ext, objs, destdir)` — X11/Motif build

Libraries built into `lib/`: `libelle.a` (base), `libelle_b.a` (batch), `libelle_wx.a` (wxWidgets), `libelle_x.a` (X11).

### Adding a new process

Copy an existing process directory and edit its `Imakefile`. Set `THISAPP`, `SRCS`, and `OBJS`, then choose which `Make*` macros to include. Always include `../../elle.tmpl`. Run `xmkmf && make Makefiles` at the top level before building.

## Running a process

```sh
./bin/elle_diff -f input.elle -n 100
```

- `-f <file>` — input `.elle` or `.poly` file
- `-n <stages>` — number of simulation stages

`ParseOptions(argc, argv)` handles these flags. If no file is given, a process can call `ElleSetFile("default.elle")` as a fallback.

## Core data model

Elle represents a microstructure as three layers:

| Layer | Name | Description |
|-------|------|-------------|
| **Flynns** | `Flynn` class (`basecode/flynns.h`) | Polygonal grains/subgrains with attributes (Euler angles, viscosity, concentration, etc.) |
| **Boundary nodes (bnodes)** | `basecode/nodes.h` | Nodes on grain boundaries; carry velocity, stress, strain attributes |
| **Unit nodes (unodes)** | `Unode` (`basecode/unodes.h`) | Interior lattice points within flynns; carry scalar/vector per-point data |

Flynns form a parent/child hierarchy (grains containing subgrains). Boundary nodes form a topology of grain boundaries — they are either `TRIPLE_J` (three-grain junction) or `DOUBLE_J` (two-grain point on a boundary segment). Unodes can be arranged on hex (`HEX_GRID`), square (`SQ_GRID`), or random (`RAN_GRID`) grids, and are associated with their enclosing flynn.

Attributes on all three object types use a generic `AttributeArray` (`basecode/attribarray.h`), addressed by integer IDs defined in `basecode/file.h` (e.g., `F_ATTRIB_A`, `N_ATTRIB_A`, `U_ATTRIB_A`, `STRESS`, `CONC_A`). Node types and state are defined in `basecode/attrib.h` (e.g., `TRIPLE_J`, `DOUBLE_J`, `ACTIVE`, `INACTIVE`, `NO_NB`).

## The Elle file format (`.elle` / `.poly`)

Text-based sections delimited by keywords (`OPTIONS`, `FLYNNS`, `PARENTS`, `LOCATION`, `VELOCITY`, `STRESS`, `CONC_A`, `F_ATTRIB_A`, etc.). See `basecode/file.h` for the full set of section identifiers. The `tidy` utility in `utilities/tidy/` cleans and normalises these files.

Runtime parameters (`Stages`, `Temperature`, `Pressure`, `Timestep`, `SaveInterval`, etc.) live in the `OPTIONS` section and map to the `runtime_opts` struct in `basecode/runopts.h`. They can also be set programmatically via `ElleSet*()` functions before `StartApp()`.

## Writing a new process

Every process follows the same pattern (see `examples/simple/`):

1. Call `ElleInit()`.
2. Register callbacks: `ElleSetInitFunction()`, `ElleSetRunFunction()`, `ElleSetExitFunction()`.
3. Call `ParseOptions(argc, argv)` to handle `-f file -n stages` flags.
4. Optionally call `SetupApp(argc, argv)` for display (guard with `if (ElleDisplay())`).
5. Call `StartApp()` — this reads the input file, runs the init callback, then iterates the run callback for the requested number of stages.
6. Call `CleanUp()` after `StartApp()` returns.

The init callback must call `ElleReinit()` to clear data structures, then `ElleSetRunFunction()` to register the per-stage function, then `ElleReadData(ElleFile())` to load the input file.

The run callback loops `EllemaxStages()` times. Each iteration processes the data and ends with `ElleUpdate()` (from `basecode/update.h`) — this increments the stage counter, triggers display refresh, and handles file saving at the configured interval. Call `ElleCheckFiles()` at the top of the run loop to handle file output checks.

```c
// Typical run callback pattern
int MyRunFunction() {
    int i, k, max;
    ElleCheckFiles();
    for (i = 0; i < EllemaxStages(); i++) {
        max = ElleMaxFlynns();
        for (k = 0; k < max; k++) {
            if (ElleFlynnIsActive(k))
                MyProcess(k);
        }
        ElleUpdate();
    }
}
```

## Error handling

Use `OnError(message, err_num)` (from `basecode/error.h`) for fatal errors. Error codes are defined in `basecode/errnum.h` (e.g., `OPEN_ERR`, `READ_ERR`, `DATA_ERR`). All Elle API functions that can fail return `int` — non-zero indicates error.

## Key API headers

- `basecode/interface.h` — flynn manipulation (split, merge, area, Euler angles, attributes, unode lists)
- `basecode/nodes.h` — bnode topology and attribute I/O
- `basecode/unodes.h` — unode initialisation (`HEX_GRID`/`SQ_GRID`/`RAN_GRID`), position, and attributes
- `basecode/runopts.h` — runtime parameters (`ElleFile()`, `EllemaxStages()`, `ElleDisplay()`, `ElleTemperature()`, etc.)
- `basecode/init.h` — `ElleInit()`, `ElleReinit()`, `ElleExit()`
- `basecode/setup.h` — `StartApp()`, `SetupApp()`, `Run_App()`
- `basecode/update.h` — `ElleUpdate()` (call once per stage at end of run loop)
- `basecode/general.h` — geometry utilities (polygon centroid, intersection, angle, coordinate rotation)
- `basecode/misorient.h` — misorientation calculations between orientations
- `basecode/splitting.h` — directed and area-based flynn splitting helpers
- `basecode/error.h` / `basecode/errnum.h` — error handling

## Existing processes (`processes/`)

| Directory | Process |
|-----------|---------|
| `diffusion` | Chemical diffusion |
| `disloc_rx` | Dislocation-driven recrystallisation |
| `exchange` | Attribute exchange between neighbours |
| `gbdiff` | Grain-boundary diffusion |
| `gbm` | Grain-boundary migration |
| `growth` | Grain growth |
| `melt` / `jkb-melt` | Melting |
| `metamorphism` | Metamorphic reactions |
| `nucleation` | Nucleation |
| `phasefield` | Phase-field modelling |
| `recovery` | Recovery |
| `split` | Grain splitting |
| `viscosity` | Viscosity assignment |
| `statistics` | Per-stage statistics output |

## Platform notes

- Primary target: Linux (WSL2 or native). Cygwin/Win32 builds use `install.win` and `Make.win`.
- Fortran libraries (`-lgfortran` or `-lg2c`) are required; `gfortran` is the default FC on modern Linux.
- wxWidgets must be installed for `wx` builds; X11 and Motif for `x` builds.
- The `utilities/gpc/` directory (General Polygon Clipper library) must exist before `install.sh` can proceed.
