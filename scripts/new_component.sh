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
template_dir="${script_dir}/templates/component_base"

shopt -s nullglob globstar

if [[ -e "${component_dir}" ]]; then
  echo "Error: ${component_dir} already exists."
  exit 1
fi

lower_name="$(printf "%s" "${component_name}" | tr '[:upper:]' '[:lower:]')"
upper_name="$(printf "%s" "${component_name}" | tr '[:lower:]' '[:upper:]')"

if [[ -d "${template_dir}" ]]; then
  cp -R "${template_dir}/." "${component_dir}/"
else
  mkdir -p "${component_dir}/include/dstl" "${component_dir}/.github/workflows" "${component_dir}/.vscode"
fi

if [[ ! -f "${component_dir}/CMakeLists.txt" ]]; then
  cat > "${component_dir}/CMakeLists.txt" <<EOF
cmake_minimum_required(VERSION 3.12)
project(${component_name} LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(dstl_${lower_name} INTERFACE)
add_library(dstl::${lower_name} ALIAS dstl_${lower_name})

target_include_directories(dstl_${lower_name}
  INTERFACE
    \$<BUILD_INTERFACE:\${CMAKE_CURRENT_SOURCE_DIR}/include>
    \$<INSTALL_INTERFACE:include>
)

target_compile_features(dstl_${lower_name} INTERFACE cxx_std_17)

set_target_properties(dstl_${lower_name} PROPERTIES EXPORT_NAME ${lower_name})

install(TARGETS dstl_${lower_name}
  EXPORT DSTLTargets
)

install(DIRECTORY include/
  DESTINATION \${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
EOF
fi

template_header="${component_dir}/include/dstl/__COMPONENT_NAME__.h"
if [[ -f "${template_header}" ]]; then
  mv "${template_header}" "${component_dir}/include/dstl/${component_name}.h"
elif [[ ! -f "${component_dir}/include/dstl/${component_name}.h" ]]; then
  cat > "${component_dir}/include/dstl/${component_name}.h" <<EOF
#pragma once

namespace dstl
{
namespace ${component_name}
{
}
} // namespace dstl
EOF
fi

if [[ ! -f "${component_dir}/.clang-format" && -f "${dstl_root}/DataPools/.clang-format" ]]; then
  cp "${dstl_root}/DataPools/.clang-format" "${component_dir}/.clang-format"
fi
if [[ ! -f "${component_dir}/.gitignore" && -f "${dstl_root}/DataPools/.gitignore" ]]; then
  cp "${dstl_root}/DataPools/.gitignore" "${component_dir}/.gitignore"
fi
if [[ ! -d "${component_dir}/.github/workflows" && -d "${dstl_root}/DataPools/.github/workflows" ]]; then
  cp -R "${dstl_root}/DataPools/.github/workflows" "${component_dir}/.github/"
fi

for file in "${component_dir}"/.github/workflows/*.yml; do
  sed -i \
    -e "s/DataPools/${component_name}/g" \
    -e "s/DATAPOOLS/${upper_name}/g" \
    -e "s/datapools/${lower_name}/g" \
    "${file}"
done

for file in "${component_dir}"/**/*; do
  if [[ -f "${file}" ]]; then
    sed -i \
      -e "s/__COMPONENT_NAME__/${component_name}/g" \
      -e "s/__COMPONENT_NAME_UPPER__/${upper_name}/g" \
      -e "s/__COMPONENT_NAME_LOWER__/${lower_name}/g" \
      "${file}"
  fi
done

echo "Created ${component_name} at ${component_dir}"
echo "Next: run scripts/publish_component.sh ${component_name} to create the repo and submodule."
