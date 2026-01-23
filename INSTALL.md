# Install ELLE (Conan-Based Workflow)

Primary deliverables:
- `elle_cli` (headless CLI)
- `latte` (process executable)
- Optional `elle_gui` (when GUI enabled and wxWidgets added to Conan recipe)

We standardize on a **Conan + CMake** toolchain for reproducibility and easy portability across machines.

---
## 1. Binary Packages (Release Artifacts)

GitHub Releases include:
- `.tar.gz` (portable)
- `.deb` (best effort, Ubuntu/Debian)

### Install .deb
```bash
sudo apt-get install -y ./elle-*.deb
elle_cli --version
```

### Use .tar.gz
```bash
tar -xzf elle-*.tar.gz
cd elle-*/
./bin/elle_cli --version
```

---
## 2. Build From Source (Reproducible)

Prerequisites:
```bash
sudo apt-get update
sudo apt-get install -y python3-pip cmake ninja-build g++ gcc
pip install --user conan
conan profile detect --force
```

Then install dependencies and build (Release, headless):
```bash
conan install . -s build_type=Release --build=missing -of build
cmake -S . -B build -G Ninja -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/elle_cli --version
```

Install to a staging prefix:
```bash
cmake --install build --prefix build/stage
ls build/stage/bin
```

---
## 3. Enabling Optional Features

Add (or uncomment) packages in `conanfile.txt` for GUI (wxWidgets) or MPI. Then re-run the Conan install and pass CMake options:
```bash
conan install . -s build_type=Release --build=missing -of build
cmake -S . -B build -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DELLE_BUILD_GUI=ON -DELLE_USE_MPI=ON
cmake --build build -j
```

Check resulting binaries:
```bash
ls build/ | grep elle_gui || echo "GUI not built"
```

---
## 4. Offline / Air-Gapped Environments

1. On an online machine, perform a Conan install + lock generation:
   ```bash
   conan install . -s build_type=Release --build=missing -of build
   conan lock create conanfile.txt -s build_type=Release --lockfile-out=conan.lock
   ```
2. Archive the needed Conan cache folders (minimal subset) and transfer them:
   ```bash
   tar -C ~/.conan2 -czf conan_cache_minimal.tgz p # or refine selection
   ```
3. On the offline machine, extract into `~/.conan2`, copy the repo + `conan.lock`, then run the same install command (no network required provided all packages exist in cache).

Do not modify `conan.lock` manually. Regenerate it intentionally when you change dependencies.

---
## 5. Packaging Locally

After a successful build:
```bash
(cd build && cpack -G TGZ)
(cd build && cpack -G DEB || true)
```

---
## 6. Verifying Installation

```bash
elle_cli --version
elle_cli --help || true
```

---
## 7. Troubleshooting

| Symptom | Action |
| ------- | ------ |
| Conan cannot find compiler | Run `conan profile detect --force` |
| Missing IPhreeqc | Ensure `iphreeqc/3.x` stays in `conanfile.txt`; rerun install |
| Offline failure | Confirm cache copy + `conan.lock`; avoid `--update` |
| Need debug build | Use `-s build_type=Debug` and `-DCMAKE_BUILD_TYPE=Debug` |

---
## 8. Version Info

`elle_cli --version` shows the version from `VERSION` or release tag.

---
## 9. License

See `LICENSE` (references legacy text under `elle/COPYING.txt`).

---
End of installation guide.
