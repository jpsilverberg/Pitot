#!/bin/bash
# Quick test script to verify a release package works correctly

set -e

VERSION=${1}
if [ -z "$VERSION" ]; then
  echo "Usage: $0 <version>"
  echo "Example: $0 1.0.0"
  exit 1
fi

PACKAGE="dstl-${VERSION}-headers.tar.gz"

if [ ! -f "$PACKAGE" ]; then
  echo "Error: Package $PACKAGE not found"
  echo "Run: ./scripts/create_release_package.sh ${VERSION}"
  exit 1
fi

echo "Testing release package: $PACKAGE"
echo ""

# Create a temporary test directory
TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

cd "$TEST_DIR"
echo "Working in: $TEST_DIR"
echo ""

# Extract the package
echo "1. Extracting package..."
tar -xzf "$OLDPWD/$PACKAGE"
echo "   ✓ Package extracted successfully"
echo ""

# Verify key files exist
echo "2. Verifying package contents..."
REQUIRED_FILES=(
  "include/dstl/dstl.h"
  "lib/cmake/DSTL/DSTLConfig.cmake"
  "lib/cmake/DSTL/DSTLTargets.cmake"
)

for file in "${REQUIRED_FILES[@]}"; do
  if [ -f "$file" ]; then
    echo "   ✓ Found: $file"
  else
    echo "   ✗ Missing: $file"
    exit 1
  fi
done
echo ""

# Count headers
HEADER_COUNT=$(find . -name "*.h" -o -name "*.hpp" | wc -l)
echo "3. Package statistics:"
echo "   Total headers: $HEADER_COUNT"
echo "   Package size: $(du -h "$OLDPWD/$PACKAGE" | cut -f1)"
echo ""

# Check for unwanted files
echo "4. Checking for unwanted files..."
UNWANTED_COUNT=0

if find . -name "*.o" -o -name "*.a" -o -name "*.so" | grep -q .; then
  echo "   ⚠ Warning: Found compiled binaries"
  UNWANTED_COUNT=$((UNWANTED_COUNT + 1))
fi

if find . -type d -name "build*" | grep -q .; then
  echo "   ⚠ Warning: Found build directories"
  UNWANTED_COUNT=$((UNWANTED_COUNT + 1))
fi

if find . -type d -name "CMakeFiles" | grep -q .; then
  echo "   ⚠ Warning: Found CMakeFiles directories"
  UNWANTED_COUNT=$((UNWANTED_COUNT + 1))
fi

if [ $UNWANTED_COUNT -eq 0 ]; then
  echo "   ✓ No unwanted files found"
fi
echo ""

# Create a test CMake project
echo "5. Testing CMake integration..."
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.14)
project(TestDSTL)

set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex)

add_executable(test_app test.cpp)
target_link_libraries(test_app PRIVATE dstl::dmutex)
EOF

cat > test.cpp << 'EOF'
#include <dstl/dmutex.h>
int main() {
    return 0;
}
EOF

mkdir build
cd build

if cmake .. > /dev/null 2>&1; then
  echo "   ✓ CMake configuration successful"
else
  echo "   ✗ CMake configuration failed"
  exit 1
fi

if cmake --build . > /dev/null 2>&1; then
  echo "   ✓ Build successful"
else
  echo "   ✗ Build failed"
  exit 1
fi
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✓ All tests passed!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Package is ready for release:"
echo "  $PACKAGE"
echo ""
echo "Next steps:"
echo "  1. git tag v${VERSION}"
echo "  2. git push origin v${VERSION}"
echo "  3. GitHub Actions will create the release automatically"
echo ""
echo "Or upload manually to:"
echo "  https://github.com/jpsilverberg/DSTL/releases"
