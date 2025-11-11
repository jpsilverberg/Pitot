# DSTL

Modern CMake glue for the DSTL component libraries and headers. This repo wraps
the component submodules (`Mutex`, `Logger`, `DataPools`, `Numbers`, `LinAlg`,
`Quadrature`, `Bits`, etc.) into a single package and exports `dstl::` targets
that downstream projects can consume like any other modern CMake dependency.

## Build & Install

```bash
cmake -S . -B build -DDSTL_CXX_STANDARD=17 -DDSTL_ENABLE_WARNINGS=ON
cmake --build build
cmake --install build --prefix "$HOME/.local"
```
Or rely on the bundled VS Code / CMake presets:
```bash
cmake --preset default
cmake --build --preset default-build
cmake --build --preset install
```

`DSTL_CXX_STANDARD` may be set to `14`, `17`, `20`, or `23`. Disable the shared
warning flags by passing `-DDSTL_ENABLE_WARNINGS=OFF` if you need a quieter
build. Tests are disabled by default; add `-DBUILD_TESTING=ON` if you need them.
AddressSanitizer support can be toggled at the top level using
`-DDSTL_ENABLE_ASAN=ON`, which propagates into submodules that honor
`ENABLE_ASAN`.

## Consuming DSTL

You have three practical options when building another project on top of DSTL:

1. **Install + `find_package` (recommended for reuse)**  
   Install DSTL once (see above commands), then add the following to your
   project:
   ```cmake
   find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)

   add_executable(app main.cpp)
   target_link_libraries(app PRIVATE dstl::dmutex dstl::dpools)
   ```
   The exported `dstl::` targets bring in include paths, compile features, and
   transitive dependencies automatically.

2. **Vendor as a submodule**  
   Add this repository as a git submodule (for example under `extern/dstl`) and
   in your `CMakeLists.txt` call `add_subdirectory(extern/dstl)`. You can then
   link to the same `dstl::` targets without installing, at the cost of a longer
   configure/generate step for every consumer project.

3. **FetchContent**  
   If you prefer not to manage the source manually, wrap the repository in
   CMake's `FetchContent`:
   ```cmake
   include(FetchContent)
   FetchContent_Declare(dstl
     GIT_REPOSITORY https://github.com/jpsilverberg/DSTL.git
     GIT_TAG main)
   FetchContent_MakeAvailable(dstl)
   ```
   This is functionally equivalent to `add_subdirectory` but downloads the
   dependency automatically during configuration.

For header-only utilities there is also a convenience interface target:
```cmake
target_link_libraries(app PRIVATE dstl)
```
which links every `dstl::` component transitively so all module headers are
available without naming individual targets.

## Available Targets

| Target            | Description                                            |
|-------------------|--------------------------------------------------------|
| `dstl`            | Interface shim that re-exports every component target. |
| `dstl::dmutex`    | Spinlock/mutex primitives from the `Mutex` module.     |
| `dstl::dlog`      | Logging utilities, depends on `dstl::dmutex`.          |
| `dstl::dpools`    | Data pool containers (`DataPools` module).             |
| `dstl::numbers`   | Arbitrary-precision/ratio types (`Numbers`).           |
| `dstl::linalg`    | Linear algebra helpers (`LinAlg`).                     |
| `dstl::quadrature`| Quadrature routines (`Quadrature`).                    |
| `dstl::dg0`       | DG(0) FEM prototype (`DGSTFEM` module, if present).    |

## Repository Layout

- `Bits/include` – low-level header utilities (e.g., `dstl/BitFlag.h`)
- `Mutex`, `Logger`, `DataPools`, `Numbers`, `LinAlg`, `Quadrature`, `DGSTFEM`
  – submodules that define build targets (`dstl::dmutex`, `dstl::dlog`, etc.)
- `include` – umbrella headers shared by multiple components
- `cmake` – generated/installed package configuration helpers

Each component installs into the same export set (`DSTLTargets`), allowing you
to cherry-pick the specific `dstl::` libraries you need or depend on the
umbrella `dstl` interface for header-only helpers.
