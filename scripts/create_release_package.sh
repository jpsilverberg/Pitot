#!/bin/bash
# Script to create a minimal release package of DSTL
# This creates an archive with just headers and CMake config files

set -e

VERSION=${1:-"1.0.0"}
BUILD_DIR="build-release"
INSTALL_DIR="install-release"
PACKAGE_NAME="dstl-${VERSION}-headers"

echo "Building DSTL release package version ${VERSION}..."

# Clean previous build
rm -rf "${BUILD_DIR}" "${INSTALL_DIR}" "${PACKAGE_NAME}.tar.gz" "${PACKAGE_NAME}.zip"

# Configure with tests disabled
cmake -S . -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -DBUILD_TESTING=OFF \
  -DDSTL_CXX_STANDARD=17

# Build (mostly a no-op for header-only libs)
cmake --build "${BUILD_DIR}"

# Install to local directory
cmake --install "${BUILD_DIR}"

# Clean up unwanted files and directories
echo "Cleaning up build artifacts from installation..."

# Remove all build directories (recursively)
find "${INSTALL_DIR}" -type d \( \
  -name "build*" -o \
  -name "Testing" -o \
  -name "CMakeFiles" -o \
  -name "Debug" -o \
  -name "Trials" -o \
  -name "Tests" -o \
  -name "doxyout" -o \
  -name "Doc" -o \
  -name "__pycache__" -o \
  -name ".continue" -o \
  -name ".vscode" \
\) -exec rm -rf {} + 2>/dev/null || true

# Remove unwanted file types
find "${INSTALL_DIR}" -type f \( \
  -name "*.gcov" -o \
  -name "Makefile" -o \
  -name "*.xml" -o \
  -name "*.sh" -o \
  -name ".gitignore" -o \
  -name ".gdbinit" -o \
  -name ".clang-format" -o \
  -name "doxygen.conf" -o \
  -name "*.md" -o \
  -name "*.txt" -o \
  -name "*.py" -o \
  -name "*.gdb" -o \
  -name "*.natvis" \
\) -delete 2>/dev/null || true

# List what's left
echo ""
echo "Package contents:"
find "${INSTALL_DIR}" -type f | head -20
TOTAL_FILES=$(find "${INSTALL_DIR}" -type f | wc -l)
echo "... (${TOTAL_FILES} files total)"
echo ""

# Create tarball
echo "Creating archive ${PACKAGE_NAME}.tar.gz..."
tar -czf "${PACKAGE_NAME}.tar.gz" -C "${INSTALL_DIR}" .

# Create zip for Windows users  
echo "Creating archive ${PACKAGE_NAME}.zip..."
(cd "${INSTALL_DIR}" && zip -q -r "../${PACKAGE_NAME}.zip" .)

# Get sizes
TAR_SIZE=$(du -h "${PACKAGE_NAME}.tar.gz" | cut -f1)
ZIP_SIZE=$(du -h "${PACKAGE_NAME}.zip" | cut -f1)

echo ""
echo "✓ Release packages created:"
echo "  - ${PACKAGE_NAME}.tar.gz (${TAR_SIZE})"
echo "  - ${PACKAGE_NAME}.zip (${ZIP_SIZE})"
echo ""
echo "Next steps:"
echo "  1. Create a release on GitHub: https://github.com/jpsilverberg/DSTL/releases/new"
echo "  2. Tag it as v${VERSION}"
echo "  3. Upload both files to the release"
echo ""
echo "Example CMakeLists.txt for consumers:"
echo "  include(FetchContent)"
echo "  FetchContent_Declare(dstl"
echo "    URL https://github.com/jpsilverberg/DSTL/releases/download/v${VERSION}/${PACKAGE_NAME}.tar.gz"
echo "  )"
echo "  FetchContent_MakeAvailable(dstl)"
echo "  find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)"
