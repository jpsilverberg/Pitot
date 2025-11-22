# Converting Geodesy to Standalone Repository

This guide walks through converting the Geodesy directory into its own git repository and then adding it back to DSTL as a submodule.

## Prerequisites

- Git installed
- GitHub account (or other git hosting)
- Write access to DSTL repository

## Step 1: Create Standalone Geodesy Repository

### 1.1: Initialize Git Repository

```bash
cd Geodesy

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: DSTL Geodesy Library v1.0.0

- WGS84 ECEF coordinate system
- Vincenty geodesic formulas
- RK4 geodesic integration
- Flat-earth and spherical approximations
- ENU/NED local frames
- English units support
- 192 tests with 100% coverage
"

# Create main branch (if not already)
git branch -M main
```

### 1.2: Create GitHub Repository

```bash
# Create repository on GitHub (via web or CLI)
gh repo create dstl-geodesy --public --source=. --remote=origin

# Or manually:
# 1. Go to GitHub and create new repository "dstl-geodesy"
# 2. Add remote:
git remote add origin https://github.com/yourusername/dstl-geodesy.git

# Push to GitHub
git push -u origin main
```

### 1.3: Tag Initial Release

```bash
# Create v1.0.0 tag
git tag -a v1.0.0 -m "Release v1.0.0

Features:
- Complete WGS84 geodetic coordinate system
- Vincenty inverse/direct formulas
- RK4 geodesic integration with altitude profiles
- Multiple coordinate systems (ECEF, Geodetic, ENU, NED)
- English units support (feet, nautical miles, knots)
- 192 tests with 100% coverage
- Production-ready quality

Accuracy:
- Sub-meter for short distances
- ~10m for medium distances (100-1000km)
- ~100m for long distances (>1000km)

Performance:
- Header-only library
- Lazy computation with caching
- All tests run in 1ms
"

# Push tag
git push origin v1.0.0
```

## Step 2: Remove Geodesy from DSTL

```bash
# Go back to DSTL root
cd ..

# Remove Geodesy directory (it's now in its own repo)
git rm -r Geodesy

# Commit removal
git commit -m "Remove Geodesy directory (moving to submodule)"
```

## Step 3: Add Geodesy as Submodule

```bash
# Add as submodule
git submodule add https://github.com/yourusername/dstl-geodesy.git Geodesy

# Commit submodule addition
git commit -m "Add Geodesy as submodule"

# Push changes
git push
```

## Step 4: Update DSTL Documentation

Update the main DSTL README.md to reference Geodesy as a submodule:

```markdown
## Components

- **Geodesy** (submodule): WGS84 geodetic coordinate system
  - Repository: https://github.com/yourusername/dstl-geodesy
  - Documentation: [Geodesy/README.md](Geodesy/README.md)
```

## Step 5: Verify Everything Works

```bash
# Clone DSTL with submodules
git clone --recursive https://github.com/yourusername/DSTL.git test-dstl
cd test-dstl

# Build
cmake -B build -DBUILD_TESTING=ON
cmake --build build

# Run tests
ctest --test-dir build -R Geodesy

# Should see: 192 tests passing
```

## Step 6: Update Submodule (Future Updates)

When you make changes to Geodesy:

```bash
# In Geodesy directory
cd Geodesy
git add .
git commit -m "Your changes"
git push

# In DSTL root
cd ..
git add Geodesy
git commit -m "Update Geodesy submodule"
git push
```

## Alternative: Keep as Subtree

If you prefer subtree over submodule:

```bash
# Remove directory
git rm -r Geodesy
git commit -m "Remove Geodesy (preparing for subtree)"

# Add as subtree
git subtree add --prefix Geodesy https://github.com/yourusername/dstl-geodesy.git main --squash

# Update subtree (future)
git subtree pull --prefix Geodesy https://github.com/yourusername/dstl-geodesy.git main --squash
```

## Benefits of Standalone Repository

### For Geodesy:
✅ Independent versioning
✅ Separate issue tracking
✅ Can be used without full DSTL
✅ Easier to contribute to
✅ Clearer dependency management
✅ Can have its own CI/CD

### For DSTL:
✅ Cleaner main repository
✅ Modular architecture
✅ Each component can evolve independently
✅ Easier to maintain
✅ Better separation of concerns

## Files Created for Standalone

- ✅ `.gitignore` - Git ignore patterns
- ✅ `LICENSE` - MIT license
- ✅ `README_STANDALONE.md` - Standalone README (rename to README.md)
- ✅ Updated `CMakeLists.txt` - Standalone + submodule support

## Checklist

- [ ] Initialize git repository in Geodesy/
- [ ] Create GitHub repository
- [ ] Push Geodesy to GitHub
- [ ] Tag v1.0.0 release
- [ ] Remove Geodesy from DSTL
- [ ] Add Geodesy as submodule
- [ ] Update DSTL documentation
- [ ] Test clean clone with submodules
- [ ] Update CI/CD if applicable

## Notes

- The CMakeLists.txt now detects if it's standalone or part of DSTL
- LinAlg dependency is handled automatically
- All tests continue to work in both modes
- No code changes required - just repository structure

## Support

For issues specific to Geodesy:
- File issues at: https://github.com/yourusername/dstl-geodesy/issues

For DSTL integration issues:
- File issues at: https://github.com/yourusername/DSTL/issues
