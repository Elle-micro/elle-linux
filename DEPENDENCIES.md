# Dependencies

This project can be built in two ways:

1) System packages (default, apt)
2) Conan (opt-in, reproducible pinned versions)

Both paths build the same targets (`elle_core`, `elle_cli`). GUI and MPI remain optional and off by default.

## 1) Apt-based build (default)

Install the basics on Ubuntu 24.04:

- build-essential cmake ninja-build pkg-config
- zlib1g-dev libgsl-dev

Optional:

- libwxgtk3.2-dev (GUI)
- libopenmpi-dev (MPI)
- liblapack-dev (LAPACK)

Build:

```bash
rm -rf build-apt \
  && cmake -S . -B build-apt -G Ninja \
  && cmake --build build-apt -j
```

## 2) Conan-based build (opt-in, lockfile-supported)

Install Conan 2.x:

```bash
pipx install conan
# or
pip install --user conan
```

Generate a profile and install deps (using the repo's lockfile when present):

```bash
conan profile detect --force
conan install . --lockfile=conan.lock --output-folder=build-conan --build=missing
```

Configure CMake with the Conan toolchain:

```bash
cmake -S . -B build-conan -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=build-conan/conan/conan_toolchain.cmake
cmake --build build-conan -j
```

Notes:

- Conan requirements are pinned in `conanfile.txt` (zlib 1.3, gsl 2.7).
- CI (ci-conan.yml) uses the committed `conan.lock` to guarantee reproducible dependency resolution across GCC/Clang and Debug/Release.
- Optional packages like wxWidgets, MPI, and LAPACK can be added later.
- When the Conan toolchain is in use, CMake prints: "Conan toolchain active".

### Conan CI path

The GitHub Actions workflow `.github/workflows/ci-conan.yml` builds a headless matrix using Conan v2 with the repository lockfile:

- Matrix: compiler ∈ {gcc, clang}, build_type ∈ {Debug, Release}
- Runner: ubuntu-24.04
- Cache: `~/.conan2` keyed by OS, compiler, build_type, and `hashFiles('conanfile.txt','conan.lock')`
- Configure flags: `-DELLE_HEADLESS=ON -DELLE_BUILD_GUI=OFF -DELLE_USE_MPI=OFF`

Local reproduction of the CI install step:

```bash
conan profile detect --force
conan install . --lockfile=conan.lock --output-folder=b-conan/conan --build=missing -s build_type=Release
cmake -S . -B b-conan -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=b-conan/conan/conan_toolchain.cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DELLE_HEADLESS=ON -DELLE_BUILD_GUI=OFF -DELLE_USE_MPI=OFF
cmake --build b-conan -j
```
