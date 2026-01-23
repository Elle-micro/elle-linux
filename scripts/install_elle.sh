#!/usr/bin/env bash
set -euo pipefail

# Portable installer for ELLE on Linux
# - Installs system prerequisites (cmake, ninja, gcc, g++) using distro pkg manager when available
# - Creates a local Python venv and installs Conan 2
# - Runs Conan install, CMake configure, build, and staged install
#
# Usage:
#   bash scripts/install_elle.sh [--prefix <path>] [--build-type <Release|Debug>] [--enable-gui] [--enable-mpi]
#
# Defaults:
#   --prefix       "$PWD/build/stage"
#   --build-type   Release
#
# Notes:
#   - You may be prompted for sudo password to install system packages.
#   - If a package manager is not detected or sudo is unavailable, the script will continue and
#     expect dependencies to be pre-installed (cmake, ninja or make, gcc/g++, python3 with venv).

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

PREFIX="$ROOT_DIR/build/stage"
BUILD_TYPE="Release"
ENABLE_GUI="OFF"
ENABLE_MPI="OFF"

print_msg() { printf "\033[1;32m==> %s\033[0m\n" "$*"; }
print_warn() { printf "\033[1;33m[WARN] %s\033[0m\n" "$*"; }
print_err() { printf "\033[1;31m[ERR] %s\033[0m\n" "$*"; }

usage() {
  cat <<EOF
Usage: $0 [options]
Options:
  --prefix <path>          Install prefix (default: $PREFIX)
  --build-type <type>      CMake build type: Release|Debug (default: $BUILD_TYPE)
  --enable-gui             Attempt to build GUI (requires wxWidgets in conanfile)
  --enable-mpi             Enable MPI support (requires MPI in conanfile)
  -h, --help               Show this help
EOF
}

while (( "$#" )); do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2 ;;
    --build-type) BUILD_TYPE="$2"; shift 2 ;;
    --enable-gui) ENABLE_GUI="ON"; shift ;;
    --enable-mpi) ENABLE_MPI="ON"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) print_err "Unknown option: $1"; usage; exit 2 ;;
  esac
done

# Detect package manager
PKG_MGR=""
if command -v apt-get >/dev/null 2>&1; then PKG_MGR="apt"; fi
if command -v dnf >/dev/null 2>&1; then PKG_MGR="dnf"; fi
if command -v yum >/dev/null 2>&1; then PKG_MGR="yum"; fi
if command -v pacman >/dev/null 2>&1; then PKG_MGR="pacman"; fi
if command -v zypper >/dev/null 2>&1; then PKG_MGR="zypper"; fi

need_sudo_install() {
  # use sudo if not root and sudo exists
  if [ "$(id -u)" -ne 0 ] && command -v sudo >/dev/null 2>&1; then
    echo "sudo"
  elif [ "$(id -u)" -eq 0 ]; then
    echo ""
  else
    echo "NOSUDO"
  fi
}

install_packages() {
  local pkgs=(cmake ninja-build gcc g++ python3 python3-venv python3-pip)
  local sudo_cmd
  sudo_cmd=$(need_sudo_install)
  if [ "$sudo_cmd" = "NOSUDO" ]; then
    print_warn "No sudo available and not root; skipping system package installation. Ensure dependencies are installed."
    return 0
  fi
  case "$PKG_MGR" in
    apt)
      $sudo_cmd apt-get update -y
      $sudo_cmd apt-get install -y "${pkgs[@]}"
      ;;
    dnf)
      $sudo_cmd dnf install -y cmake ninja-build gcc gcc-c++ python3 python3-pip
      ;;
    yum)
      # On some RHEL/CentOS, cmake3 may be the package name
      $sudo_cmd yum install -y cmake || $sudo_cmd yum install -y cmake3 || true
      $sudo_cmd yum install -y ninja-build gcc gcc-c++ python3 python3-pip || true
      ;;
    pacman)
      $sudo_cmd pacman -Sy --noconfirm cmake ninja gcc python python-pip
      ;;
    zypper)
      $sudo_cmd zypper refresh
      $sudo_cmd zypper install -y cmake ninja gcc gcc-c++ python3 python3-pip
      ;;
    *)
      print_warn "Unsupported or unknown package manager. Skipping automatic dependency installation."
      ;;
  esac
}

ensure_tools() {
  local missing=()
  command -v cmake >/dev/null 2>&1 || missing+=(cmake)
  # prefer ninja when present
  if ! command -v ninja >/dev/null 2>&1 && ! command -v make >/dev/null 2>&1; then
    missing+=(ninja)
  fi
  command -v gcc >/dev/null 2>&1 || missing+=(gcc)
  command -v g++ >/dev/null 2>&1 || missing+=(g++)
  command -v python3 >/dev/null 2>&1 || missing+=(python3)
  if ((${#missing[@]})); then
    print_warn "Missing system packages: ${missing[*]}"
    # Ask for confirmation to install
    local ans
    read -r -p "Install missing packages using your system package manager? [Y/n] " ans || ans=""
    ans=${ans:-Y}
    if [[ "$ans" =~ ^[Yy]$ ]]; then
      # Explicitly request sudo password upfront to avoid multiple prompts
      if [ "$(id -u)" -ne 0 ]; then
        if command -v sudo >/dev/null 2>&1; then
          print_msg "Requesting sudo credentials to install dependencies (one-time)"
          if ! sudo -v; then
            print_err "sudo authentication failed. Cannot install dependencies automatically."
            exit 1
          fi
        else
          print_err "'sudo' not found and you are not root. Please install dependencies manually and re-run."
          exit 1
        fi
      fi
      print_msg "Installing missing system packages: ${missing[*]}"
      install_packages || true
    else
      print_err "Dependencies not installed. Aborting as requested."
      exit 1
    fi
  fi
  # Re-check critical tools
  for t in cmake gcc g++; do
    if ! command -v "$t" >/dev/null 2>&1; then
      print_err "Required tool '$t' is missing. Please install it and re-run."
      exit 1
    fi
  done
}

create_venv_and_conan() {
  local venv_dir="$ROOT_DIR/.venv"
  if [ ! -d "$venv_dir" ]; then
    # Send status output to stderr so callers capturing stdout only get the path
    print_msg "Creating Python venv at $venv_dir" >&2
    python3 -m venv "$venv_dir"
  fi

  local py="$venv_dir/bin/python"
  local pip="$venv_dir/bin/pip"
  local conan_bin="$venv_dir/bin/conan"

  # Ensure pip tooling is up to date in the venv
  "$py" -m pip install -U pip wheel >/dev/null

  # Install Conan into the venv if it's not present yet
  if [ ! -x "$conan_bin" ]; then
    print_msg "Installing Conan into venv" >&2
    "$py" -m pip install -U "conan>=2.0,<3.0" >/dev/null
  fi

  # Resolve the conan binary path robustly (some environments may differ)
  if [ ! -x "$conan_bin" ]; then
    conan_bin="$($py -c 'import shutil; p=shutil.which("conan"); print(p or "")')"
  fi

  if [ -z "${conan_bin:-}" ] || [ ! -x "$conan_bin" ]; then
    print_err "Conan was not found in the virtual environment after installation." >&2
    print_err "Please ensure Python/pip can install packages and try again." >&2
    exit 1
  fi

  # Only echo the resolved conan path on stdout for clean command substitution
  echo "$conan_bin"
}

configure_build() {
  local conan_bin="$1"
  local build_dir="$ROOT_DIR/build"
  mkdir -p "$build_dir"
  print_msg "Detecting Conan profile"
  "$conan_bin" profile detect --force >/dev/null

  # Remove stale lockfile if present to avoid "Requirement not in lockfile" errors
  if [ -f "conan.lock" ]; then
    print_warn "Removing existing conan.lock to ensure fresh dependency resolution"
    rm "conan.lock"
  fi

  print_msg "Installing Conan dependencies (build_type=$BUILD_TYPE)"
  "$conan_bin" install . -s build_type="$BUILD_TYPE" --build=missing -of "$build_dir"

  # Choose generator
  local gen
  if command -v ninja >/dev/null 2>&1; then
    gen="Ninja"
  else
    gen="Unix Makefiles"
  fi
  print_msg "Configuring CMake (-G $gen)"
  cmake -S . -B "$build_dir" -G "$gen" \
    -DCMAKE_TOOLCHAIN_FILE="$build_dir/conan_toolchain.cmake" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DELLE_BUILD_GUI="$ENABLE_GUI" -DELLE_USE_MPI="$ENABLE_MPI"
}

build_targets() {
  local jobs
  jobs=$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)
  print_msg "Building (parallel jobs: $jobs)"
  cmake --build "$ROOT_DIR/build" -j "$jobs"
}

install_stage() {
  print_msg "Installing to prefix: $PREFIX"
  cmake --install "$ROOT_DIR/build" --prefix "$PREFIX"
}

verify_install() {
  print_msg "Verifying installation"
  if [ -x "$PREFIX/bin/elle_cli" ]; then
    "$PREFIX/bin/elle_cli" --version || true
  else
    print_warn "elle_cli not found under $PREFIX/bin"
  fi
  ls -l "$PREFIX/bin" || true
}

main() {
  print_msg "Starting ELLE installation"
  ensure_tools
  local conan_bin
  conan_bin=$(create_venv_and_conan)
  configure_build "$conan_bin"
  build_targets
  install_stage
  verify_install
  cat <<EONOTE

Done. Binaries are in: $PREFIX/bin
You can run:
  $PREFIX/bin/elle_cli --version

Optionally add to PATH for this session:
  export PATH="$PREFIX/bin:\$PATH"

Re-run with --build-type Debug for debug builds, or --enable-gui/--enable-mpi if your conanfile.txt includes those deps.
EONOTE
}

main "$@"
