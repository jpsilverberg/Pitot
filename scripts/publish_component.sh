#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <ComponentName>"
  echo "Example: $0 StringUtils"
}

if [[ ${1:-} == "" ]]; then
  usage
  exit 1
fi

component_name="$1"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dstl_root="$(cd "${script_dir}/.." && pwd)"
component_dir="${dstl_root}/${component_name}"
lower_name="$(printf "%s" "${component_name}" | tr '[:upper:]' '[:lower:]')"

if [[ ! -d "${component_dir}" ]]; then
  echo "Error: ${component_dir} does not exist. Run new_component.sh first."
  exit 1
fi

cd "${component_dir}"

if [[ ! -d ".git" ]]; then
  git init
  git branch -M main
fi

if ! git rev-parse --verify HEAD >/dev/null 2>&1; then
  git add .
  git commit -m "Initial ${component_name} library"
fi

repo_fullname="jpsilverberg/${component_name}"

if gh repo view "${repo_fullname}" >/dev/null 2>&1; then
  if ! git remote get-url origin >/dev/null 2>&1; then
    git remote add origin "https://github.com/${repo_fullname}.git"
  fi
  git push -u origin main
else
  gh repo create "${component_name}" --private --source=. --remote=origin --push
fi

cd "${dstl_root}"
git submodule add -f "https://github.com/${repo_fullname}.git" "${component_name}"
git add "${dstl_root}/.gitmodules" "${dstl_root}/${component_name}"

export DSTL_ROOT="${dstl_root}"
export COMPONENT_NAME="${component_name}"
export LOWER_NAME="${lower_name}"

python3 - <<'PY'
import os

dstl_root = os.environ["DSTL_ROOT"]
component = os.environ["COMPONENT_NAME"]
lower = os.environ["LOWER_NAME"]

def insert_into_block(path, block_name, new_item, before_item=None):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    try:
        start = next(i for i, l in enumerate(lines) if l.strip().startswith(f"set({block_name}"))
        end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == ")")
    except StopIteration:
        return

    block = lines[start + 1:end]
    if any(l.strip() == new_item for l in block):
        return

    insert_at = end
    if before_item:
        for i in range(start + 1, end):
            if lines[i].strip() == before_item:
                insert_at = i
                break

    lines.insert(insert_at, f"  {new_item}")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

def ensure_include(path, header):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    include_line = f'#include "{header}"'
    if any(l.strip() == include_line for l in lines):
        return
    lines.append(include_line)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

root_cmake = os.path.join(dstl_root, "CMakeLists.txt")
dstl_header = os.path.join(dstl_root, "include", "dstl", "dstl.h")
forward_header = os.path.join(dstl_root, "include", "dstl", f"{component}.h")

insert_into_block(root_cmake, "DSTL_COMPONENT_DIRS", component, "DGSTFEM")
insert_into_block(root_cmake, "DSTL_INTERFACE_LINK_TARGETS", f"dstl::{lower}", "dstl::dg0")

os.makedirs(os.path.dirname(forward_header), exist_ok=True)
if not os.path.exists(forward_header):
    with open(forward_header, "w", encoding="utf-8", newline="\n") as f:
        f.write("#pragma once\n\n")
        f.write(f'#include "../../{component}/include/dstl/{component}.h"\n')

ensure_include(dstl_header, f"{component}.h")
PY

echo "Submodule added and DSTL integration updated. Commit in DSTL with:"
echo "  git commit -m \"Add ${component_name} submodule\""
