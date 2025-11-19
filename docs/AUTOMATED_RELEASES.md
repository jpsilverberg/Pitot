# Automated Release Workflows

DSTL includes two GitHub Actions workflows to automate the creation and distribution of release packages. Choose the workflow that best fits your release process.

## Quick Start

### Option 1: Fully Automatic (Recommended)
Just push a version tag:
```bash
git tag v1.2.3
git push origin v1.2.3
```
GitHub Actions will automatically:
- Build the release package
- Generate checksums
- Create the GitHub release
- Upload all artifacts
- Add usage instructions to the release notes

### Option 2: Manual Release with Auto-Upload
1. Create a release manually on GitHub
2. GitHub Actions will automatically build and upload the packages

### Option 3: Manual Test Build
Trigger the workflow manually from the Actions tab to test without creating a release.

---

## Workflow Details

### 1. Auto-Release on Tag (`.github/workflows/auto-release.yml`)

**Trigger:** Automatically runs when you push a tag matching `v*.*.*` (e.g., `v1.0.0`, `v2.3.1`)

**What it does:**
1. ✅ Checks out code with submodules
2. ✅ Builds release packages (`.tar.gz` and `.zip`)
3. ✅ Generates SHA256 checksums
4. ✅ Creates release notes with usage examples
5. ✅ Creates a GitHub Release
6. ✅ Uploads all artifacts automatically
7. ✅ Includes git changelog since last tag

**Usage:**
```bash
# Make sure your changes are committed
git add .
git commit -m "Prepare for v1.2.3 release"
git push

# Create and push a tag
git tag v1.2.3
git push origin v1.2.3

# GitHub Actions does the rest automatically!
# Watch progress at: https://github.com/jpsilverberg/DSTL/actions
```

**Advantages:**
- ✅ Fully automated - just push a tag
- ✅ Consistent release notes format
- ✅ Automatic changelog generation
- ✅ No manual GitHub UI interaction needed

**Best for:** Regular releases, CI/CD pipelines, quick releases

---

### 2. Release Package Upload (`.github/workflows/release-package.yml`)

**Triggers:**
- When you manually create a GitHub Release
- When manually triggered from the Actions tab

**What it does:**
1. ✅ Detects version from release tag or manual input
2. ✅ Builds release packages
3. ✅ Generates checksums
4. ✅ Uploads to the release (if triggered by release creation)
5. ✅ Or saves as workflow artifacts (if manually triggered)

**Usage - Automatic on Release:**
```bash
# Push your code
git push

# Go to GitHub:
# 1. Navigate to https://github.com/jpsilverberg/DSTL/releases/new
# 2. Create a new tag (e.g., v1.2.3)
# 3. Write release notes
# 4. Click "Publish release"
# GitHub Actions automatically uploads the packages
```

**Usage - Manual Trigger:**
```bash
# Go to: https://github.com/jpsilverberg/DSTL/actions/workflows/release-package.yml
# Click "Run workflow"
# Enter version (e.g., 1.2.3)
# Click "Run workflow"
# Download artifacts from the workflow run
```

**Advantages:**
- ✅ Full control over release notes
- ✅ Can test package creation without making a release
- ✅ Useful for pre-releases or drafts

**Best for:** Custom release notes, pre-releases, testing

---

## Comparison Table

| Feature | Auto-Release | Manual Release |
|---------|-------------|----------------|
| Trigger | Push tag | Create release or manual |
| Automation | Fully automatic | Semi-automatic |
| Release Notes | Auto-generated | You write them |
| Changelog | Automatic | Manual |
| Best for | Quick releases | Custom notes |
| Steps | 1 (push tag) | 2+ (GitHub UI) |

---

## Release Checklist

### Before Creating a Release

- [ ] All changes committed and pushed
- [ ] Tests passing locally
- [ ] Version bumped in relevant files (if applicable)
- [ ] Submodules updated to latest (if needed)
- [ ] CHANGELOG.md updated (if you maintain one)

### Creating the Release

**Using Auto-Release (Recommended):**
```bash
# 1. Ensure clean state
git status

# 2. Create and push tag
git tag v1.2.3 -m "Release v1.2.3"
git push origin v1.2.3

# 3. Monitor the action
# Visit: https://github.com/jpsilverberg/DSTL/actions

# 4. Verify release
# Visit: https://github.com/jpsilverberg/DSTL/releases
```

**Using Manual Release:**
```bash
# 1. Push code
git push

# 2. Create release on GitHub
# https://github.com/jpsilverberg/DSTL/releases/new

# 3. Wait for Actions to upload packages

# 4. Edit release notes if needed
```

### After Release

- [ ] Verify packages are attached to release
- [ ] Check checksums file is present
- [ ] Test download URL works
- [ ] Test consumption in a sample project:
  ```cmake
  FetchContent_Declare(dstl
    URL https://github.com/jpsilverberg/DSTL/releases/download/v1.2.3/dstl-1.2.3-headers.tar.gz
  )
  ```
- [ ] Announce release (if applicable)

---

## Troubleshooting

### Workflow Fails to Upload

**Problem:** "Resource not accessible by integration" error

**Solution:** Ensure GitHub Actions has write permissions:
1. Go to Settings → Actions → General
2. Under "Workflow permissions", select "Read and write permissions"
3. Save and re-run the workflow

### Wrong Version Number

**Problem:** Package has wrong version number

**Solution:**
- For auto-release: Ensure tag matches `v*.*.*` format (e.g., `v1.2.3`)
- For manual: Check the version input when triggering the workflow
- Delete the tag and recreate: `git tag -d v1.2.3 && git push origin :refs/tags/v1.2.3`

### Missing Submodules in Package

**Problem:** Package is missing component headers

**Solution:**
- Ensure `.gitmodules` is properly configured
- Both workflows check out with `submodules: recursive`
- If submodules are private, may need to configure deploy keys

### Package Too Large

**Problem:** Package exceeds GitHub's file size limits

**Solution:**
- Review cleanup patterns in `scripts/create_release_package.sh`
- Add more exclusions for build artifacts
- Consider splitting into multiple archives if truly necessary

### Release Already Exists

**Problem:** Can't create release because tag already exists

**Solution:**
```bash
# Delete the release on GitHub first, then:
git tag -d v1.2.3
git push origin :refs/tags/v1.2.3
# Now create it again
git tag v1.2.3
git push origin v1.2.3
```

---

## Advanced Configuration

### Customizing Auto-Release Notes

Edit `.github/workflows/auto-release.yml` in the "Generate release notes" step to customize the format.

### Adding Pre-release Support

Modify the tag pattern to distinguish stable vs pre-release:

```yaml
on:
  push:
    tags:
      - 'v*.*.*'        # Stable: v1.2.3
      - 'v*.*.*-rc*'    # Release candidate: v1.2.3-rc1
      - 'v*.*.*-beta*'  # Beta: v1.2.3-beta1
```

Then set `prerelease: true` conditionally:

```yaml
- name: Create GitHub Release
  uses: softprops/action-gh-release@v1
  with:
    prerelease: ${{ contains(github.ref, '-rc') || contains(github.ref, '-beta') }}
```

### Running Tests Before Release

Add a test step before building packages:

```yaml
- name: Run tests
  run: |
    cmake -S . -B build -DBUILD_TESTING=ON
    cmake --build build
    ctest --test-dir build --output-on-failure
```

### Multi-Platform Packages

To create packages for different platforms:

```yaml
strategy:
  matrix:
    os: [ubuntu-latest, windows-latest, macos-latest]
runs-on: ${{ matrix.os }}
```

---

## Example: Complete Automated Workflow

```bash
# Day-to-day development
git checkout main
git pull
# ... make changes ...
git add .
git commit -m "Add feature X"
git push

# Ready to release?
# Just tag it:
git tag v1.2.3 -m "Release v1.2.3: Added feature X"
git push origin v1.2.3

# GitHub Actions automatically:
# ✓ Builds packages
# ✓ Creates release
# ✓ Uploads artifacts
# ✓ Adds release notes
# ✓ Done!

# Users can immediately consume:
# FetchContent_Declare(dstl
#   URL .../v1.2.3/dstl-1.2.3-headers.tar.gz
# )
```

That's it! No manual uploading, no clicking through GitHub UI, just `git tag` and `git push`. 🚀
