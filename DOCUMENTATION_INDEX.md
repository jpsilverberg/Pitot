# DSTL Documentation Index

Welcome to the DSTL (Data Structures and Tools Library) documentation hub. This index provides quick access to comprehensive documentation for all DSTL components, organized by functional category.

## Quick Navigation

- [Getting Started](#getting-started)
- [Component Selection Guide](#component-selection-guide)
- [Core Components](#core-components)
- [Improvement Notes](#improvement-notes)
- [Building and Installation](#building-and-installation)

---

## Getting Started

DSTL is a collection of C++ libraries designed for mathematical and systems programming. Each component is independently documented and can be used standalone or in combination with other DSTL components.

**New to DSTL?** Start here:
1. Read the [main README](README.md) for build and installation instructions
2. Browse the [Component Selection Guide](#component-selection-guide) below to find the right component for your needs
3. Dive into the specific component documentation linked in the [Core Components](#core-components) section

**Looking for something specific?**
- Memory management → [DataPools](DataPools/README.md)
- Thread synchronization → [Mutex](Mutex/README.md)
- Logging and diagnostics → [Logger](Logger/README.md)
- Fixed-point arithmetic → [Numbers](Numbers/README.md)
- Linear algebra → [LinAlg](LinAlg/README.md)
- Numerical integration → [Quadrature](Quadrature/README.md)
- Utility types → [Bits](Bits/README.md)

---

## Component Selection Guide

### When to Use Each Component

**DataPools** - Use when you need:
- High-performance memory pool allocation
- Reference-counted object management
- Deterministic memory usage patterns
- Custom memory allocation strategies
- Thread-safe or single-threaded object pools

**Mutex** - Use when you need:
- Thread synchronization primitives
- Low-overhead spinlocks for short critical sections
- RAII-based lock management
- Configurable thread safety policies
- Performance-critical concurrent code

**Logger** - Use when you need:
- Structured logging with multiple severity levels
- Performance tracking and function profiling
- Scoped logging contexts
- Thread-safe logging infrastructure
- Flexible output configuration (files, streams)

**Numbers** - Use when you need:
- Fixed-point arithmetic for embedded systems
- Deterministic numeric behavior
- Arbitrary precision calculations
- Ratio and interval arithmetic
- Probability and accumulator types

**LinAlg** - Use when you need:
- Vector and matrix operations
- Linear system solving
- LU decomposition
- Matrix determinants and inverses
- Compile-time sized linear algebra

**Quadrature** - Use when you need:
- Numerical integration
- Definite integral approximation
- Configurable accuracy vs. performance tradeoffs
- Integration of complex functions

**Bits** - Use when you need:
- Type-safe bit flag operations
- Result/error handling types
- Status code hierarchies
- Lightweight utility types

---

## Core Components

### Memory Management

#### [DataPools](DataPools/README.md)
Memory pool system with reference counting and multiple allocation strategies. Provides smart pointer-like accessors (Ptr, Ref, Shared, Owned, Iter) for managing pooled objects with automatic lifetime management.

**Key Features:**
- Multiple pool types: NormalPool, FixedPool, LinkedPool, HashPool, HostPool, ClientPool
- Reference-counted object lifecycle management
- Thread-safe and single-threaded policies
- Host/Client pattern for distributed object management
- Zero-overhead abstractions for pool access

**[View Improvement Notes](DataPools/IMPROVEMENTS.md)**

---

### Concurrency

#### [Mutex](Mutex/README.md)
Thread synchronization primitives including spinlocks, thread-aware locks, and RAII lock guards. Provides configurable mutex types optimized for different use cases.

**Key Features:**
- NoOpMutex for single-threaded contexts
- Spinlock for low-latency critical sections
- ThreadAwareSpinlock for recursive locking
- LockGuard RAII wrapper for exception-safe locking
- Policy-based mutex selection (SingleThreadMutex, MultiThreadMutex)

**[View Improvement Notes](Mutex/IMPROVEMENTS.md)**

---

### Utilities

#### [Logger](Logger/README.md)
Comprehensive logging system with multiple severity levels, scoped contexts, and performance tracking. Thread-safe by default with flexible output configuration.

**Key Features:**
- Six log levels: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
- DLOG singleton for global logging
- LogInstance for scoped logging contexts
- DFuncLog for automatic function performance tracking
- Convenient macros: LOG, DEBUG, ERR, WARN, VAR
- File and stream output support

**[View Improvement Notes](Logger/IMPROVEMENTS.md)**

#### [Bits](Bits/README.md)
Utility types for common programming patterns including type-safe bit flags and result/error handling.

**Key Features:**
- BitFlag template for type-safe enum-based flags
- Result class with hierarchical status codes
- Compile-time type safety
- Zero-overhead abstractions
- Clean error handling patterns

**[View Improvement Notes](Bits/IMPROVEMENTS.md)**

---

### Mathematics

#### [Numbers](Numbers/README.md)
Fixed-point and arbitrary precision arithmetic library with automatic type promotion and comprehensive mathematical operations.

**Key Features:**
- Multiple numeric types: FixedPoint, Ratio, Scaled, Interval, Accumulator, Prob
- Automatic type promotion for mixed arithmetic
- Saturation and overflow handling
- Mathematical functions (sin, cos, sqrt, etc.)
- Compile-time precision configuration

**[View Improvement Notes](Numbers/IMPROVEMENTS.md)**

#### [LinAlg](LinAlg/README.md)
Linear algebra library providing vector and matrix operations with compile-time sizing and efficient implementations.

**Key Features:**
- Vec class for vector operations
- Mat class for matrix operations
- LU decomposition
- Linear system solving
- Determinant and inverse calculations
- Template-based compile-time sizing

**[View Improvement Notes](LinAlg/IMPROVEMENTS.md)**

#### [Quadrature](Quadrature/README.md)
Numerical integration library implementing various quadrature methods for approximating definite integrals.

**Key Features:**
- Multiple quadrature methods
- Configurable accuracy and performance
- Support for various function types
- Efficient integration algorithms
- Error estimation capabilities

**[View Improvement Notes](Quadrature/IMPROVEMENTS.md)**

---

## Improvement Notes

Each component has an associated improvement notes document that identifies potential enhancements, performance optimizations, and areas for future development. These notes are valuable for maintainers and contributors.

| Component | Improvement Notes |
|-----------|-------------------|
| DataPools | [IMPROVEMENTS.md](DataPools/IMPROVEMENTS.md) |
| Mutex | [IMPROVEMENTS.md](Mutex/IMPROVEMENTS.md) |
| Logger | [IMPROVEMENTS.md](Logger/IMPROVEMENTS.md) |
| Numbers | [IMPROVEMENTS.md](Numbers/IMPROVEMENTS.md) |
| LinAlg | [IMPROVEMENTS.md](LinAlg/IMPROVEMENTS.md) |
| Quadrature | [IMPROVEMENTS.md](Quadrature/IMPROVEMENTS.md) |
| Bits | [IMPROVEMENTS.md](Bits/IMPROVEMENTS.md) |

---

## Building and Installation

For detailed build and installation instructions, see the [main README](README.md).

### Quick Start

```bash
# Build with CMake
cmake -S . -B build -DDSTL_CXX_STANDARD=17
cmake --build build
cmake --install build --prefix "$HOME/.local"
```

### Using DSTL in Your Project

```cmake
# Find and link DSTL components
find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)

add_executable(myapp main.cpp)
target_link_libraries(myapp PRIVATE dstl::dmutex dstl::dpools)
```

### Available CMake Targets

| Target | Component | Description |
|--------|-----------|-------------|
| `dstl::dmutex` | Mutex | Thread synchronization primitives |
| `dstl::dlog` | Logger | Logging utilities |
| `dstl::dpools` | DataPools | Memory pool management |
| `dstl::numbers` | Numbers | Fixed-point arithmetic |
| `dstl::linalg` | LinAlg | Linear algebra operations |
| `dstl::quadrature` | Quadrature | Numerical integration |
| `dstl` | All | Umbrella target for all components |

---

## Additional Resources

- [Release Packages Documentation](docs/RELEASE_PACKAGES.md)
- [Automated Release Workflow](docs/AUTOMATED_RELEASES.md)
- [Release Quick Reference](docs/RELEASE_QUICK_REFERENCE.md)

---

## Contributing

When contributing to DSTL components, please:
1. Review the relevant component's README for usage patterns
2. Check the IMPROVEMENTS.md file for known issues and planned enhancements
3. Ensure your changes maintain consistency with existing APIs
4. Add tests for new functionality
5. Update documentation to reflect your changes

---

*Last Updated: November 2024*
