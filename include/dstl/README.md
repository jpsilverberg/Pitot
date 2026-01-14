# DSTL Forwarding Headers

This directory contains **forwarding headers** that provide a unified include location for all DSTL components.

## What Are Forwarding Headers?

Forwarding headers are thin wrapper headers that include the actual library-specific headers using relative paths. They allow users to write simple, clean includes like:

```cpp
#include <dstl/Logger.h>
#include <dstl/Numbers.h>
#include <dstl/LinAlg.h>
```

Without needing to know the internal directory structure of the DSTL repository.

## How It Works

### During Source Builds

When building DSTL from source, CMake configures **two include directories**:

1. `/include` - This directory (forwarding headers)
2. `LibraryName/include` - Library-specific headers

Example from `LinAlg/CMakeLists.txt`:
```cmake
target_include_directories(dstl_linalg
  PUBLIC
    $<BUILD_INTERFACE:${DSTL_ROOT}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
```

### Forwarding Pattern

Each forwarding header uses a **relative path** to include the actual library header:

**Example:** [`/include/dstl/LinAlg.h`](file:///home/jon/src/DSTL/include/dstl/LinAlg.h)
```cpp
#pragma once

#include "../../LinAlg/include/dstl/LinAlg.h"
```

The relative path `../../` navigates from `/include/dstl/` up to the project root, then down to `LinAlg/include/dstl/`.

### After Installation

When DSTL is installed (e.g., to `/usr/local`), all library headers are copied to a **single flat directory**:

```
/usr/local/include/dstl/
├── Logger.h         # From Logger/include/dstl/Logger.h
├── Numbers.h        # From Numbers/include/dstl/Numbers.h
├── LinAlg.h         # From LinAlg/include/dstl/LinAlg.h
└── ...              # All other library headers
```

The forwarding headers from this directory are **not installed** because they become redundant - all headers are already in the same location.

## Why Use Relative Paths?

The forwarding headers use relative paths (`../../LibraryName/include/dstl/Header.h`) instead of angle brackets (`<dstl/Header.h>`) to:

1. **Avoid ambiguity** - Both include paths are active during builds, relative paths explicitly choose the library-specific version
2. **Prevent circular includes** - Ensures the forwarding header includes the implementation, not itself
3. **Make intent clear** - Explicitly shows this is a forwarding relationship

## Adding a New Library

When adding a new DSTL library, follow these steps to maintain consistency:

### 1. Create the library directory structure

```
NewLibrary/
├── CMakeLists.txt
├── include/dstl/
│   ├── NewLibrary.h      # Main library header
│   └── ...               # Other headers
├── src/
│   └── ...
└── Tests/
    └── ...
```

### 2. Create the forwarding header

Create `/include/dstl/NewLibrary.h`:

```cpp
#pragma once

#include "../../NewLibrary/include/dstl/NewLibrary.h"
```

### 3. Use PascalCase naming

Follow the naming convention:
- ✅ `Logger.h`, `Mutex.h`, `DataPools.h`, `LinAlg.h`
- ❌ `dlog.h`, `dmutex.h`, `dpools.h` (old style, being phased out)

### 4. In library code, use angle brackets

Within your library's implementation files, always use angle brackets:

```cpp
// In NewLibrary/include/dstl/SomeHeader.h
#include <dstl/Logger.h>     // ✅ Correct
#include <dstl/Numbers.h>    // ✅ Correct

// NOT relative paths:
#include "../../Logger/include/dstl/Logger.h"  // ❌ Wrong
```

### 5. Configure CMake include directories

In `NewLibrary/CMakeLists.txt`:

```cmake
target_include_directories(dstl_newlibrary
  PUBLIC
    $<BUILD_INTERFACE:${DSTL_ROOT}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
```

This makes both forwarding headers and library headers available.

## Current Forwarding Headers

| Forwarding Header | Library Header | Library |
|-------------------|----------------|---------|
| `Logger.h` | `Logger/include/dstl/Logger.h` | Logger |
| `Mutex.h` | `Mutex/include/dstl/Mutex.h` | Mutex |
| `DataPools.h` | `DataPools/include/dstl/DataPools.h` | DataPools |
| `Numbers.h` | `Numbers/include/dstl/Numbers.h` | Numbers |
| `LinAlg.h` | `LinAlg/include/dstl/LinAlg.h` | LinAlg |
| `Quadrature.h` | `Quadrature/include/dstl/Quadrature.h` | Quadrature |
| `StringUtils.h` | `StringUtils/include/dstl/StringUtils.h` | StringUtils |
| `Pitot.h` | `Pitot/include/dstl/Pitot.h` | Pitot |
| `dprint.h` | `StringUtils/include/dstl/dprint.h` | StringUtils (utility) |
| `BitFlag.h` | `Bits/include/dstl/BitFlag.h` | Bits |
| `Result.h` | `Bits/include/dstl/Result.h` | Bits |

## See Also

- [ARCHITECTURE.md](file:///home/jon/src/DSTL/ARCHITECTURE.md) - Complete DSTL architecture documentation
- [CMakeLists.txt](file:///home/jon/src/DSTL/CMakeLists.txt) - Root build configuration
