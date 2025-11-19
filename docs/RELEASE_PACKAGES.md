# Creating DSTL Release Packages

This document explains how to create and publish lightweight DSTL release packages that contain only headers and CMake configuration files, without test infrastructure or build artifacts.

## Why Release Packages?

Release packages provide several benefits for consumers:

1. **Faster downloads** - Only headers and config files, not full repository history
2. **Faster CMake configuration** - No need to configure/build test suites
3. **Simpler dependency management** - Similar to how GoogleTest can be consumed
4. **Reproducible builds** - Fixed version tags ensure consistency

## Creating a Release Package

### Using the Script

```bash
./scripts/create_release_package.sh <version>
```

Example:
```bash
./scripts/create_release_package.sh 1.2.3
```

This generates:
- `dstl-1.2.3-headers.tar.gz` - For Linux/macOS users
- `dstl-1.2.3-headers.zip` - For Windows users

### What's Included

The release package contains:
- All public headers from `include/`, `Bits/include/`
- Component headers from `Mutex/include/`, `Logger/include/`, etc.
- CMake config files (`DSTLConfig.cmake`, `DSTLConfigVersion.cmake`, `DSTLTargets.cmake`)
- Minimal source files needed for compiled components

### What's Excluded

The release package excludes:
- All `build*` directories
- Test suites (`Tests/`, `Testing/`)
- Documentation build files (`doxyout/`, `Doc/`)
- IDE configurations (`.vscode/`)
- Development tools (`.gitignore`, `.clang-format`, Makefiles)
- Coverage reports and debugging helpers

## Publishing to GitHub Releases

### Manual Process

1. Go to https://github.com/jpsilverberg/DSTL/releases/new
2. Create a new tag (e.g., `v1.2.3`)
3. Set the release title (e.g., "DSTL v1.2.3")
4. Upload both `dstl-1.2.3-headers.tar.gz` and `dstl-1.2.3-headers.zip`
5. Publish the release

### Automated Process (GitHub Actions)

The repository includes a GitHub Actions workflow that automatically creates release packages when you publish a release:

1. Go to https://github.com/jpsilverberg/DSTL/releases/new
2. Create a tag and publish the release
3. GitHub Actions will automatically build and attach the release packages

You can also trigger the workflow manually from the Actions tab.

## Consuming Release Packages

Once published, users can consume DSTL without cloning the full repository:

### FetchContent with Release Archive

```cmake
cmake_minimum_required(VERSION 3.14)
project(MyProject)

include(FetchContent)
FetchContent_Declare(dstl
  URL https://github.com/jpsilverberg/DSTL/releases/download/v1.2.3/dstl-1.2.3-headers.tar.gz
  # Optional: verify integrity
  # URL_HASH SHA256=<checksum-of-tarball>
)
FetchContent_MakeAvailable(dstl)

find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)

add_executable(myapp main.cpp)
target_link_libraries(myapp PRIVATE dstl::dmutex dstl::dpools)
```

### Advantages Over Git Clone

Compare the release package approach:

```cmake
# Release package - fast, minimal
FetchContent_Declare(dstl
  URL https://github.com/jpsilverberg/DSTL/releases/download/v1.2.3/dstl-1.2.3-headers.tar.gz
)
```

To the traditional git approach:

```cmake
# Git clone - slower, includes everything
FetchContent_Declare(dstl
  GIT_REPOSITORY https://github.com/jpsilverberg/DSTL.git
  GIT_TAG v1.2.3
  OPTIONS "BUILD_TESTING OFF"
)
```

The release package is typically 10-50x smaller and configures much faster.

## Verifying Package Contents

Before publishing, verify the package contents:

```bash
# List contents of tarball
tar -tzf dstl-1.2.3-headers.tar.gz | head -20

# Extract to inspect
mkdir test-extract
tar -xzf dstl-1.2.3-headers.tar.gz -C test-extract
tree test-extract
```

Expected structure:
```
test-extract/
├── include/
│   └── dstl/
│       ├── BitFlag.h
│       ├── dstl.h
│       ├── types/
│       └── ...
├── lib/
│   └── cmake/
│       └── DSTL/
│           ├── DSTLConfig.cmake
│           ├── DSTLConfigVersion.cmake
│           └── DSTLTargets.cmake
├── DataPools/
│   ├── include/
│   └── src/
├── Logger/
│   ├── include/
│   └── src/
└── Mutex/
    ├── include/
    └── src/
```

## Troubleshooting

### Package Too Large

If the package is unexpectedly large, check for included build artifacts:

```bash
tar -tzf dstl-1.2.3-headers.tar.gz | grep -E '(build|\.o$|\.gcov|CMakeFiles)'
```

If you find build artifacts, update the cleanup patterns in `scripts/create_release_package.sh`.

### CMake Can't Find Package

If consumers report CMake can't find the package:

1. Verify `DSTLConfig.cmake` is in the tarball
2. Check the package version matches the tag
3. Ensure `find_package(DSTL CONFIG ...)` uses `CONFIG` mode

### Missing Headers

If headers are missing from the package:

1. Check that the component's CMakeLists.txt properly installs headers
2. Verify the install commands in the root CMakeLists.txt
3. Test a local install: `cmake --install build-release --prefix test-install`

## Best Practices

1. **Version Consistency** - Always use semantic versioning (MAJOR.MINOR.PATCH)
2. **Test Before Publishing** - Create a test project that consumes the package
3. **Document Changes** - Include a changelog in the GitHub release notes
4. **Provide Checksums** - Include SHA256 hashes in release notes for security
5. **Keep Archives Small** - Regularly review and update cleanup patterns

## Example: Complete Release Workflow

```bash
# 1. Ensure repo is clean and tagged
git status
git tag v1.2.3

# 2. Create release package
./scripts/create_release_package.sh 1.2.3

# 3. Generate checksums
sha256sum dstl-1.2.3-headers.tar.gz > checksums.txt
sha256sum dstl-1.2.3-headers.zip >> checksums.txt

# 4. Test the package locally
mkdir /tmp/test-dstl-release
cd /tmp/test-dstl-release
# Create a test CMakeLists.txt and verify it works

# 5. Push tag and create GitHub release
git push origin v1.2.3
# Then go to GitHub and create the release, uploading the archives

# 6. Update documentation
# Add release notes, migration guide, etc.
```
