# Elle — Installation Guide

Elle is a C/C++ microstructural modelling framework. It uses the imake/X11 build system and links against Fortran, GSL, zlib, and optionally wxWidgets or X11/Motif for display.

---

## Ubuntu / Debian

### 1. Install system packages

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    g++ \
    gfortran \
    imake \
    xutils-dev \
    libgsl-dev \
    zlib1g-dev
```

For **batch (headless) builds** that is all you need.

For **wxWidgets display builds** (`wx`):

```bash
sudo apt-get install -y \
    libwxgtk3.2-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev
```

> On Ubuntu 20.04/22.04 the package is `libwxgtk3.0-gtk3-dev`. On 24.04 it is `libwxgtk3.2-dev`. Check with `apt-cache search libwxgtk`.

For **X11/Motif display builds** (`x`):

```bash
sudo apt-get install -y \
    libx11-dev \
    libxt-dev \
    libmotif-dev
```

> `libmotif-dev` provides OpenMotif. If unavailable, try `lesstif2-dev` as a substitute.

### 2. Get the source

```bash
git clone https://github.com/Elle-micro/elle-linux elle
cd elle/elle
```

Or unpack from a tarball:

```bash
tar xzf elle.tar.gz
cd elle/elle
```

### 3. Build

```bash
./install.sh batch        # headless executables only → binb/
./install.sh wx           # wxWidgets display executables → binwx/
./install.sh x            # X11/Motif display executables → binx/
./install.sh all          # all three variants
```

`bin/` is symlinked to the last variant built.

### 4. Verify

```bash
./bin/elle_diff --help
```

---

## macOS

### 1. Install Homebrew

If you do not have Homebrew:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Follow the printed instructions to add Homebrew to your PATH. On Apple Silicon (M1/M2/M3) the prefix is `/opt/homebrew`; on Intel it is `/usr/local`.

```bash
# Apple Silicon
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"

# Intel
echo 'eval "$(/usr/local/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/usr/local/bin/brew shellenv)"
```

### 2. Install dependencies

```bash
brew install gcc          # provides gfortran (gfortran-14 or similar)
brew install gsl
brew install imake        # provides xmkmf
```

For **wxWidgets display builds** (`wx`):

```bash
brew install wxwidgets
```

> X11/Motif builds (`x`) are not supported on macOS. Use `batch` or `wx`.

### 3. Make gfortran visible as `gfortran`

Homebrew installs gfortran as a versioned binary (e.g. `gfortran-14`). Create a symlink:

```bash
# adjust the version number to match what brew installed
ln -s $(brew --prefix)/bin/gfortran-14 $(brew --prefix)/bin/gfortran
```

Verify:

```bash
gfortran --version
```

### 4. Get the source

```bash
git clone https://github.com/Elle-micro/elle-linux elle
cd elle/elle
```

### 5. Build

```bash
./install.sh batch        # headless executables only → binb/
./install.sh wx           # wxWidgets display executables → binwx/
```

### 6. Verify

```bash
./bin/elle_diff --help
```

---

## Rebuilding after code changes

From a process or utility subdirectory:

```bash
make depend
make
make install
```

To rebuild only the base library after changing `basecode/`:

```bash
cd basecode && make && make install_base
```

---

## Library summary

| Library | Purpose | Ubuntu package | macOS (Homebrew) |
|---------|---------|---------------|-----------------|
| g++ / gcc | C++ compiler | `build-essential` | Xcode CLT (`xcode-select --install`) |
| gfortran | Fortran runtime | `gfortran` | `brew install gcc` |
| GSL | Maths library | `libgsl-dev` | `brew install gsl` |
| zlib | Compression | `zlib1g-dev` | built-in |
| imake/xmkmf | Build system | `imake`, `xutils-dev` | `brew install imake` |
| wxWidgets | GUI (wx builds) | `libwxgtk3.2-dev` | `brew install wxwidgets` |
| OpenGL/GLU | Rendering (wx) | `libgl1-mesa-dev`, `libglu1-mesa-dev` | built-in (macOS) |
| X11/Motif | GUI (x builds) | `libx11-dev`, `libmotif-dev` | not supported |
