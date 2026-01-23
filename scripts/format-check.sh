#!/usr/bin/env bash
# Format check script (changed-files mode with exclusions)
# Environment variables:
#   FORMAT_BASE_REF   Git reference to diff against (default: origin/master or main fallback)
#   FORMAT_EXCLUDE_REGEX  Egrep regex of paths to exclude (default: 'wxplotcode/')
#   FORMAT_FIX        If set to 1, apply clang-format changes in-place instead of just checking

set -euo pipefail

BASE_REF=${FORMAT_BASE_REF:-}
if [ -z "${BASE_REF}" ]; then
	if git rev-parse --verify origin/master >/dev/null 2>&1; then
		BASE_REF=origin/master
	elif git rev-parse --verify origin/main >/dev/null 2>&1; then
		BASE_REF=origin/main
	else
		# Fallback: last commit (this lets initial commit still format-check staged changes)
		BASE_REF=HEAD~1
	fi
fi

EXCLUDE_REGEX=${FORMAT_EXCLUDE_REGEX:-wxplotcode/}

# Collect changed C/C++ files vs base (added, modified, renamed). Use diff-filter=ACMR
mapfile -t changed < <(git diff --name-only --diff-filter=ACMR "${BASE_REF}" -- '*.c' '*.cc' '*.cpp' '*.h' '*.hpp') || true

if [ ${#changed[@]} -eq 0 ]; then
	echo "No changed C/C++ files relative to ${BASE_REF}."; exit 0
fi

# Filter out excluded paths
filtered=()
for f in "${changed[@]}"; do
	if [[ -n "${EXCLUDE_REGEX}" ]] && echo "${f}" | egrep -q "${EXCLUDE_REGEX}"; then
		continue
	fi
	if [ -f "${f}" ]; then
		filtered+=("${f}")
	fi
done

if [ ${#filtered[@]} -eq 0 ]; then
	echo "No changed C/C++ files to check after exclusions."; exit 0
fi

clang-format --version

if [ "${FORMAT_FIX:-0}" = "1" ]; then
	echo "Applying clang-format to: ${filtered[*]}"
	clang-format -i "${filtered[@]}"
	echo "Re-running check to ensure clean state..."
fi

# Build argument list safely
fail=0
for f in "${filtered[@]}"; do
	if ! clang-format --dry-run --Werror "${f}" >/dev/null 2>&1; then
		echo "FORMAT VIOLATION: ${f}" >&2
		fail=1
	fi
done

if [ ${fail} -ne 0 ]; then
	echo "One or more files need formatting." >&2
	echo "Tip: FORMAT_FIX=1 $0 to auto-fix just the changed files." >&2
	exit 1
fi

echo "Formatting OK for changed files.";
