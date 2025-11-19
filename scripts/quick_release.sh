#!/bin/bash
# Quick release helper script
# This guides you through creating a release with the automated workflow

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   DSTL Automated Release Helper       ║${NC}"
echo -e "${BLUE}╔════════════════════════════════════════╗${NC}"
echo ""

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ] || [ ! -d ".github" ]; then
  echo -e "${RED}Error: Must be run from DSTL root directory${NC}"
  exit 1
fi

# Check git status
echo -e "${YELLOW}Checking repository status...${NC}"
if ! git diff-index --quiet HEAD --; then
  echo -e "${RED}Error: You have uncommitted changes${NC}"
  echo "Please commit or stash your changes first."
  git status --short
  exit 1
fi

# Check if on main branch (optional warning)
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$CURRENT_BRANCH" != "main" ]; then
  echo -e "${YELLOW}⚠ Warning: You're on branch '$CURRENT_BRANCH', not 'main'${NC}"
  read -p "Continue anyway? (y/N) " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  fi
fi

echo -e "${GREEN}✓ Repository is clean${NC}"
echo ""

# Get version number
echo -e "${BLUE}Enter the version number (e.g., 1.2.3):${NC}"
read -p "Version: " VERSION

if [ -z "$VERSION" ]; then
  echo -e "${RED}Error: Version cannot be empty${NC}"
  exit 1
fi

# Validate version format
if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo -e "${YELLOW}⚠ Warning: Version doesn't match X.Y.Z format${NC}"
  read -p "Continue anyway? (y/N) " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  fi
fi

TAG="v${VERSION}"

# Check if tag already exists
if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo -e "${RED}Error: Tag $TAG already exists${NC}"
  echo "To delete it: git tag -d $TAG && git push origin :refs/tags/$TAG"
  exit 1
fi

echo ""
echo -e "${BLUE}Release Summary:${NC}"
echo "  Version: $VERSION"
echo "  Tag: $TAG"
echo "  Branch: $CURRENT_BRANCH"
echo "  Commit: $(git rev-parse --short HEAD)"
echo ""

# Optional: Test build
read -p "Test release package build locally first? (Y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Nn]$ ]]; then
  echo -e "${YELLOW}Building test package...${NC}"
  if ! ./scripts/create_release_package.sh "$VERSION"; then
    echo -e "${RED}Error: Package build failed${NC}"
    exit 1
  fi
  
  echo ""
  echo -e "${YELLOW}Testing package...${NC}"
  if ! ./scripts/test_release_package.sh "$VERSION"; then
    echo -e "${RED}Error: Package tests failed${NC}"
    exit 1
  fi
  
  # Clean up test artifacts
  rm -rf build-release install-release
  
  echo -e "${GREEN}✓ Package build and tests passed${NC}"
fi

echo ""
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Ready to create release!${NC}"
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "This will:"
echo "  1. Create and push tag $TAG"
echo "  2. Trigger GitHub Actions to:"
echo "     - Build release packages"
echo "     - Create GitHub Release"
echo "     - Upload artifacts automatically"
echo ""
read -p "Proceed with automated release? (y/N) " -n 1 -r
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  echo -e "${YELLOW}Release cancelled${NC}"
  exit 0
fi

# Create and push the tag
echo ""
echo -e "${YELLOW}Creating tag $TAG...${NC}"
git tag -a "$TAG" -m "Release $TAG"

echo -e "${YELLOW}Pushing tag to origin...${NC}"
git push origin "$TAG"

echo ""
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}✓ Tag pushed successfully!${NC}"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "GitHub Actions is now building your release."
echo ""
echo "Monitor progress at:"
echo -e "${BLUE}https://github.com/jpsilverberg/DSTL/actions${NC}"
echo ""
echo "Release will be available at:"
echo -e "${BLUE}https://github.com/jpsilverberg/DSTL/releases/tag/${TAG}${NC}"
echo ""
echo "Consumers can use:"
echo -e "${GREEN}FetchContent_Declare(dstl${NC}"
echo -e "${GREEN}  URL https://github.com/jpsilverberg/DSTL/releases/download/${TAG}/dstl-${VERSION}-headers.tar.gz${NC}"
echo -e "${GREEN})${NC}"
echo ""
