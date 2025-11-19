# Quick Release Reference

## TL;DR - Fastest Way to Release

```bash
./scripts/quick_release.sh
# Follow the prompts - it handles everything!
```

Or manually:
```bash
git tag v1.0.0
git push origin v1.0.0
# GitHub Actions does the rest
```

---

## Available Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `quick_release.sh` | Interactive release wizard | `./scripts/quick_release.sh` |
| `create_release_package.sh` | Build release package locally | `./scripts/create_release_package.sh 1.0.0` |
| `test_release_package.sh` | Test a release package | `./scripts/test_release_package.sh 1.0.0` |

---

## Automated Workflows

### Auto-Release (`.github/workflows/auto-release.yml`)
- **Trigger:** Push a tag like `v1.0.0`
- **Result:** Automatic GitHub Release with packages

### Manual Upload (`.github/workflows/release-package.yml`)
- **Trigger:** Create a release on GitHub, or run manually
- **Result:** Packages uploaded to release

---

## Complete Example

```bash
# 1. Make your changes
git add .
git commit -m "Add new feature"
git push

# 2. Create release (choose one method)

# Method A: Interactive (recommended for first-timers)
./scripts/quick_release.sh

# Method B: Direct (for experienced users)
git tag v1.2.3
git push origin v1.2.3

# 3. Wait for GitHub Actions (2-3 minutes)
# Visit: https://github.com/jpsilverberg/DSTL/actions

# 4. Verify release
# Visit: https://github.com/jpsilverberg/DSTL/releases

# Done! Users can now:
# FetchContent_Declare(dstl
#   URL .../v1.2.3/dstl-1.2.3-headers.tar.gz
# )
```

---

## Testing Before Release

```bash
# Build package locally
./scripts/create_release_package.sh 1.2.3

# Test it
./scripts/test_release_package.sh 1.2.3

# If all good, release it
git tag v1.2.3
git push origin v1.2.3
```

---

## Troubleshooting

**Q: Workflow failed?**
- Check Settings → Actions → General → Workflow permissions
- Must be "Read and write permissions"

**Q: Wrong version in package?**
- Tag must match `v*.*.*` format (e.g., `v1.2.3`)
- Delete and recreate: `git tag -d v1.2.3 && git push origin :refs/tags/v1.2.3`

**Q: Want to test without releasing?**
- Use: `./scripts/create_release_package.sh 1.2.3`
- Or trigger workflow manually from Actions tab

---

## Full Documentation

- [Automated Releases Guide](AUTOMATED_RELEASES.md) - Complete workflow documentation
- [Release Packages Guide](RELEASE_PACKAGES.md) - Package creation details
- [Main README](../README.md) - General DSTL documentation

---

## Quick Links

- [Create New Release](https://github.com/jpsilverberg/DSTL/releases/new)
- [View Actions](https://github.com/jpsilverberg/DSTL/actions)
- [View Releases](https://github.com/jpsilverberg/DSTL/releases)
