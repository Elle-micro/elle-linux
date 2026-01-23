#!/usr/bin/env bash
set -euo pipefail
build_dir="${1:-build}"
if [ ! -f "${build_dir}/compile_commands.json" ]; then
  echo "compile_commands.json not found in ${build_dir}. Configure with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
  exit 1
fi
clang-tidy --version
# Limit to a handful of files initially to avoid timeouts; expand later.
mapfile -t files < <(jq -r '.[].file' "${build_dir}/compile_commands.json" | grep -E '\\.(c|cc|cpp|h|hpp)$' | sort -u)
if [ "${#files[@]}" -eq 0 ]; then exit 0; fi
clang-tidy -p "${build_dir}" "${files[@]:0:200}"
