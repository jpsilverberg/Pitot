#!/bin/bash
set -e

echo "🚀 Converting Geodesy to Standalone Repository"
echo "================================================"
echo ""

# Configuration
GEODESY_REPO_NAME="Geodesy"
GITHUB_USER="jpsilverberg"
TEMP_DIR="/tmp/geodesy-standalone"

# Step 1: Create temporary directory and copy Geodesy
echo "📦 Step 1: Preparing standalone Geodesy repository..."
rm -rf "$TEMP_DIR"
mkdir -p "$TEMP_DIR"
cp -r Geodesy/* "$TEMP_DIR/"
cp Geodesy/.gitignore "$TEMP_DIR/"

# Step 2: Initialize git repo in temp directory
echo "🔧 Step 2: Initializing git repository..."
cd "$TEMP_DIR"
git init
git add .
git commit -m "Initial commit: DSTL Geodesy Library v1.0.0

High-precision geodesy library for Earth-based coordinate transformations.
Extracted from DSTL (Dynamic Simulation Template Library).

Features:
- WGS84 ellipsoid model
- ECEF ↔ LLA conversions
- Vincenty's formulae for geodesic calculations
- Flat-Earth approximations for local dynamics
- 100% test coverage with 192 passing tests"

# Step 3: Rename README_STANDALONE.md to README.md
echo "📝 Step 3: Setting up standalone README..."
if [ -f "README_STANDALONE.md" ]; then
    mv README_STANDALONE.md README.md
    git add README.md
    git commit -m "Use standalone README as main README"
fi

# Step 4: Create GitHub repository
echo "🌐 Step 4: Creating GitHub repository..."
echo "   Repository: $GITHUB_USER/$GEODESY_REPO_NAME"
gh repo create "$GEODESY_REPO_NAME" --public --source=. --remote=origin --description="High-precision geodesy library for Earth-based coordinate transformations"

# Step 5: Push to GitHub
echo "⬆️  Step 5: Pushing to GitHub..."
git push -u origin main

# Step 6: Create and push tag
echo "🏷️  Step 6: Creating release tag v1.0.0..."
git tag -a v1.0.0 -m "Release v1.0.0

Initial standalone release of DSTL Geodesy Library.

Features:
- WGS84 ellipsoid model with high-precision constants
- ECEF ↔ LLA coordinate transformations
- Vincenty's formulae for geodesic calculations
- Flat-Earth approximations for local dynamics
- Comprehensive test suite (192 tests, 100% coverage)
- 1ms total test runtime

Validated against reference implementations and real-world data."
git push origin v1.0.0

echo ""
echo "✅ Standalone repository created successfully!"
echo "   URL: https://github.com/$GITHUB_USER/$GEODESY_REPO_NAME"
echo ""

# Step 7: Return to DSTL and convert to submodule
echo "🔄 Step 7: Converting DSTL/Geodesy to submodule..."
cd -  # Return to DSTL root

# Save current state
echo "   Committing any pending changes..."
git add -A
git commit -m "Prepare for Geodesy submodule conversion" || echo "   (No changes to commit)"

# Remove Geodesy directory from git tracking
echo "   Removing Geodesy from git tracking..."
git rm -rf Geodesy
git commit -m "Remove Geodesy directory (converting to submodule)"

# Add as submodule
echo "   Adding Geodesy as submodule..."
git submodule add "https://github.com/$GITHUB_USER/$GEODESY_REPO_NAME.git" Geodesy
git commit -m "Add Geodesy as submodule

Geodesy is now maintained as a standalone library at:
https://github.com/$GITHUB_USER/$GEODESY_REPO_NAME

This allows independent versioning and easier reuse across projects."

echo ""
echo "🎉 Conversion Complete!"
echo "======================="
echo ""
echo "✅ Geodesy is now a standalone repository"
echo "✅ DSTL now uses Geodesy as a submodule"
echo ""
echo "Next steps:"
echo "1. Push DSTL changes: git push"
echo "2. Test the build: cd build && cmake .. && make"
echo "3. Update submodule: git submodule update --init --recursive"
echo ""
echo "Standalone repo: https://github.com/$GITHUB_USER/$GEODESY_REPO_NAME"
