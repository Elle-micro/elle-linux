# .devcontainer/Dockerfile
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Base C++ toolchain + zsh + python
RUN apt-get update && apt-get install -y \
    build-essential g++-14 clang-18 clang++-18 cmake ninja-build \
    ccache pkg-config cppcheck doxygen graphviz ripgrep git curl ca-certificates \
    zsh fonts-powerline locales tzdata \
    python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Locale
RUN locale-gen en_US.UTF-8 && update-locale LANG=en_US.UTF-8
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# --- Conan in an isolated venv (avoids PEP 668) ---
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip conan \
    && ln -s /opt/venv/bin/conan /usr/local/bin/conan
ENV PATH="/opt/venv/bin:${PATH}"

# --- Oh My Zsh (non-interactive), theme + plugins ---
ENV ZSH="/root/.oh-my-zsh"
RUN git clone --depth=1 https://github.com/ohmyzsh/ohmyzsh.git "$ZSH" \
    && cp "$ZSH"/templates/zshrc.zsh-template /root/.zshrc \
    && git clone --depth=1 https://github.com/zsh-users/zsh-autosuggestions "$ZSH/custom/plugins/zsh-autosuggestions" \
    && git clone --depth=1 https://github.com/zsh-users/zsh-syntax-highlighting "$ZSH/custom/plugins/zsh-syntax-highlighting" \
    && sed -i 's/^ZSH_THEME=.*/ZSH_THEME="agnoster"/' /root/.zshrc \
    && sed -i 's/^plugins=(git)/plugins=(git zsh-autosuggestions zsh-syntax-highlighting)/' /root/.zshrc \
    && echo '\nexport CLICOLOR=1\nexport LSCOLORS=ExFxBxDxCxegedabagacad\n' >> /root/.zshrc \
    && chsh -s /usr/bin/zsh root

# Sanity prints
RUN cmake --version && g++-14 --version && clang++-18 --version && zsh --version && conan --version

WORKDIR /workspaces
