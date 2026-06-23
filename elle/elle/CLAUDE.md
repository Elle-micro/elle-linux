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

`bin/` is symlinked to the last selected variant. The script calls `xmkmf`, then `make Makefiles`, `make install_base`, then the appropriate install targets.

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

## Repository structure

```
basecode/     Core library: data structures, I/O, topology management
processes/    Individual simulation processes (one subdirectory each)
statscode/    Statistics utilities
plotcode/     X11 plot rendering
wxplotcode/   wxWidgets plot rendering
utilities/    Helper tools (convert, tidy, view, ranorient, etc.)
examples/     Example processes and .elle files
lib/          Compiled library archives (libelle*.a)
bin/          Symlink to active binary set (binb, binwx, or binx)
binb/         Headless executables
binwx/        wxWidgets executables
binx/         X11/Motif executables
```

## Core data model

Elle represents a microstructure as three layers:

| Layer | Name | Description |
|-------|------|-------------|
| **Flynns** | `Flynn` class (`basecode/flynns.h`) | Polygonal grains/subgrains with attributes (Euler angles, viscosity, concentration, etc.) |
| **Boundary nodes (bnodes)** | `basecode/nodes.h` | Nodes on grain boundaries; carry velocity, stress, strain attributes |
| **Unit nodes (unodes)** | `Unode` (`basecode/unodes.h`) | Interior lattice points within flynns; carry scalar/vector per-point data |

Flynns form a parent/child hierarchy (grains containing subgrains). Boundary nodes are linked into a topology of grain boundaries and triple/double junctions. Unodes can be arranged on hex, square, or random grids, and are associated with their enclosing flynn.

Attributes on all three object types use a generic `AttributeArray` (`basecode/attribarray.h`), addressed by integer IDs defined in `basecode/file.h` (e.g., `F_ATTRIB_A`, `N_ATTRIB_A`, `U_ATTRIB_A`, `STRESS`, `CONC_A`).

## The Elle file format (`.elle` / `.poly`)

Text-based sections delimited by keywords (`OPTIONS`, `FLYNNS`, `PARENTS`, `LOCATION`, `VELOCITY`, `STRESS`, `CONC_A`, `F_ATTRIB_A`, etc.). See `basecode/file.h` for the full set of section identifiers. The `tidy` utility in `utilities/tidy/` cleans and normalises these files.

## Writing a new process

Every process follows the same pattern (see `examples/simple/`):

1. Call `ElleInit()`.
2. Register callbacks: `ElleSetInitFunction()`, `ElleSetRunFunction()`, `ElleSetExitFunction()`.
3. Call `ParseOptions(argc, argv)` to handle `-f file -n stages` flags.
4. Optionally call `SetupApp(argc, argv)` for display.
5. Call `StartApp()` — this reads the input file, runs the init callback, then iterates the run callback for the requested number of stages.

The run callback typically iterates over all active flynns (`ElleMaxFlynns()` / `ElleFlynnIsActive()`), boundary nodes, or unodes, modifies attributes, and calls `ElleUpdateNodes()` to propagate topology changes.

Key API headers:
- `basecode/interface.h` — flynn manipulation (split, merge, area, Euler angles, attributes)
- `basecode/nodes.h` — bnode topology and attribute I/O
- `basecode/unodes.h` — unode initialisation, position, and attributes
- `basecode/runopts.h` — runtime parameters (`ElleFile()`, `ElleStages()`, `ElleDisplay()`, etc.)
- `basecode/init.h` — `ElleInit()`, `ElleReinit()`, `ElleExit()`
- `basecode/setup.h` — `StartApp()`, `SetupApp()`, `Run_App()`

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
