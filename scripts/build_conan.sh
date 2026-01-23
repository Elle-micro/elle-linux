#!/usr/bin/env bash
set -euo pipefail

# Unified Conan + CMake build helper.
# Usage examples:
#   scripts/build_conan.sh                # Release build (default dir: build)
#   scripts/build_conan.sh debug          # Debug build
#   scripts/build_conan.sh relwithdebinfo # RelWithDebInfo build
# Options via env vars:
#   BUILD_DIR=build-gui ELLE_GUI=ON ELLE_MPI=OFF scripts/build_conan.sh
#   PACKAGE=1 scripts/build_conan.sh
#   INSTALL_PREFIX=build/stage scripts/build_conan.sh
#   WARN_AS_ERROR=ON scripts/build_conan.sh
#   CLEAN=1 scripts/build_conan.sh
#   LOCK=1 scripts/build_conan.sh          # generate/update lockfile
#   USE_LOCK=1 scripts/build_conan.sh       # use existing lockfile
# Offline mode: ensure ~/.conan2 already seeded and (optionally) USE_LOCK=1.

BUILD_TYPE_INPUT=${1:-release}
BUILD_TYPE_CANON=$(echo "$BUILD_TYPE_INPUT" | tr '[:upper:]' '[:lower:]')
case "$BUILD_TYPE_CANON" in
  release)      BUILD_TYPE_CONAN=Release; BUILD_TYPE_CMAKE=Release;;
  debug)        BUILD_TYPE_CONAN=Debug; BUILD_TYPE_CMAKE=Debug;;
  relwithdebinfo|relwithdebug|relwithdeb) BUILD_TYPE_CONAN=RelWithDebInfo; BUILD_TYPE_CMAKE=RelWithDebInfo;;
  minsizerel|minsize|min) BUILD_TYPE_CONAN=MinSizeRel; BUILD_TYPE_CMAKE=MinSizeRel;;
  *) echo "Unsupported build type: $BUILD_TYPE_INPUT" >&2; exit 2;;
esac

BUILD_DIR=${BUILD_DIR:-build}
INSTALL_PREFIX=${INSTALL_PREFIX:-$BUILD_DIR/stage}
ELLE_GUI=${ELLE_GUI:-OFF}
ELLE_MPI=${ELLE_MPI:-OFF}
WARN_AS_ERROR=${WARN_AS_ERROR:-OFF}
CLEAN=${CLEAN:-0}
PACKAGE=${PACKAGE:-0}
LOCK=${LOCK:-0}
USE_LOCK=${USE_LOCK:-0}

if [[ $CLEAN == 1 ]]; then
  echo "[clean] Removing $BUILD_DIR" >&2
  rm -rf "$BUILD_DIR"
fi

CONAN_LOCK_ARGS=()
if [[ $USE_LOCK == 1 && -f conan.lock ]]; then
  CONAN_LOCK_ARGS+=(--lockfile=conan.lock)
fi

if [[ $LOCK == 1 ]]; then
  echo "[lock] Generating lockfile" >&2
  conan lock create conanfile.txt -s build_type=$BUILD_TYPE_CONAN \
    --lockfile-out=conan.lock "${CONAN_LOCK_ARGS[@]}"
fi

# Conan install (will populate $BUILD_DIR)
conan install . -s build_type=$BUILD_TYPE_CONAN --build=missing -of "$BUILD_DIR" \
  "${CONAN_LOCK_ARGS[@]}"

# CMake configure
cmake -S . -B "$BUILD_DIR" -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE="$BUILD_DIR"/conan_toolchain.cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE_CMAKE \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DELLE_BUILD_GUI=$ELLE_GUI -DELLE_USE_MPI=$ELLE_MPI -DELLE_WARN_AS_ERROR=$WARN_AS_ERROR

# Build
cmake --build "$BUILD_DIR" -j "$(nproc)"

# Install (staging)
cmake --install "$BUILD_DIR" --prefix "$INSTALL_PREFIX"

if [[ $PACKAGE == 1 ]]; then
  echo "[package] Running CPack (TGZ + DEB)" >&2
  (cd "$BUILD_DIR" && cpack -G TGZ)
  (cd "$BUILD_DIR" && cpack -G DEB || true)
fi

echo "[done] Binaries in $BUILD_DIR; staged install in $INSTALL_PREFIX" >&2

# Print versions
"$BUILD_DIR"/elle_cli --version || true
